#!/usr/bin/env python
# coding: utf-8
"""
Driver/worker variant of New-VAE-gridsearch-latest.py.

Why this exists
----------------
TensorFlow/CUDA does not reliably release GPU (or host) memory back to the OS
within a long-running process, even after `del`, `gc.collect()`, and
`tf.keras.backend.clear_session()`. Those calls only free Python-level
references and reset the Keras graph state, but the CUDA allocator keeps its
memory pool for the lifetime of the process. This causes memory to build up
across grid-search combos until the machine crashes.

The fix: run each hyperparameter combo in its own **subprocess**. When a
subprocess exits, the OS reclaims 100% of its memory (RAM + GPU), guaranteed.
This script is the "driver" — it never imports TensorFlow or touches the
dataset. It just loops over combos and spawns one child process per combo
via `subprocess.run(...)`, using this same file's "worker mode"
(triggered by a hidden worker argument) to do the actual training.

Usage
-----
    python New-VAE-gridsearch-subprocess.py <your_data.pkl> -f <output_folder>

This single command runs the *entire* grid search — no need for a PowerShell
restart loop. Progress is checkpointed to `<output_folder>/.checkpoint.json`,
so if the machine crashes or you Ctrl+C, re-running the same command resumes
from where it left off.
"""
import os, sys, json, hashlib, argparse, subprocess
from collections import Counter
from datetime import date
from itertools import product


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('train', metavar='training data', type=str, help="Path to MoD")
    p.add_argument('-f', '--savefld', dest='savefld', type=str, default=date.today().isoformat(), help="Output folder (default: today's date)")
    p.add_argument('-p', '--params', '--txt-params', dest='params_file', type=str, default=None,
                   help="Path to a JSON file containing parameter combinations. If omitted, the built-in grid is used and written to grid-search-params.json.")
    p.add_argument('-e', '--patience', dest='patience', type=int, default=10, help="Early stopping patience (default: 10 epochs)")
    p.add_argument('-d', '--min-delta', dest='min_delta', type=float, default=0.0,
                   help="Early stopping min_delta (default: 0.0). Recorded in provenance.json.")
    p.add_argument('-b', '--batch-size', dest='batch_size', type=int, default=256, help="Batch size for training (default: 256)")
    p.add_argument('-l', '--learning-rate', dest='learning_rate', type=float, default=1e-3,
                   help="Adam learning rate (default: 1e-3). Recorded in provenance.json.")
    p.add_argument('--gpu', dest='gpu', type=str, default='0', help="Which GPU to use (default: 0)")
    p.add_argument('--nondeterministic-ops', dest='nondeterministic_ops', action='store_true',
                   help="Skip tf.config.experimental.enable_op_determinism(). Faster, but runs are "
                        "no longer bitwise reproducible. The choice is recorded in provenance.json.")
    p.add_argument('--free-bits', dest='free_bits', type=float, default=0.001,
                   help="Nats per latent dimension below which a dimension is not penalised "
                        "(default: 0.001, safety valve). Only meaningful once the KL term is strong "
                        "enough relative to reconstruction for posterior collapse to be possible, "
                        "that is once rec_weight is small enough.")
    p.add_argument('--restart-from', dest='restart_from', type=int, default=None,
                   help="Restart the grid search from this displayed iteration number (1-based) and rerun every iteration from there onward. If omitted, resume normally from the checkpoint.")
    p.add_argument('--worker-idx', dest='worker_idx', type=int, default=None,
                   help=argparse.SUPPRESS)
    return p.parse_args()


def get_checkpoint_path(savefld):
    return os.path.join(savefld, '.checkpoint.json')


CHECKPOINT_VERSION = 2
PARAMETER_COUNT = 6


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

def combo_seed(combo) -> int:
    """Derive a reproducible seed from the parameter tuple itself.

    Keying the seed on the tuple rather than on the grid index means the seed
    survives an edit to the grid, for the same reason the checkpoint stores
    tuples rather than indices. Two runs of the same combination therefore
    reproduce each other, while different combinations, including different
    replicate indices of the same architecture, get unrelated seeds.

    Note that the tuple is serialised before hashing, so 3072 and 3072.0 are
    different combinations and will produce different seeds.
    """
    payload = json.dumps(list(combo), sort_keys=True).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "big") % (2 ** 31 - 1)


def set_seed(seed: int, deterministic_ops: bool = True) -> None:
    """Make a training run bitwise reproducible on the same hardware.

    Must be called before the train/test split and before the model is
    constructed, since both consume random state. Imports are local so that
    the driver process never loads TensorFlow.

    `enable_op_determinism()` costs some training speed and requires that any
    tf.data shuffling carry an explicit `seed=` argument.
    """
    import random
    import numpy as np
    import tensorflow as tf

    random.seed(seed)
    np.random.seed(seed)
    tf.keras.utils.set_random_seed(seed)  # python, numpy and TF in one call
    if deterministic_ops:
        try:
            tf.config.experimental.enable_op_determinism()
        except Exception as exc:  # pragma: no cover - depends on TF build
            print(f"WARNING: enable_op_determinism() failed ({exc}). "
                  "Runs will not be bitwise reproducible.")


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def file_sha256(path: str, chunk_size: int = 1 << 20) -> str:
    """SHA-256 of a file, read in chunks so large inputs are not slurped."""
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b''):
            digest.update(chunk)
    return digest.hexdigest()


def git_provenance() -> dict:
    """Record the commit and whether the working tree was clean.

    A commit hash alone is misleading when the tree has uncommitted changes,
    so the dirty flag and the modified paths are recorded alongside it.
    """
    repo_dir = os.path.dirname(os.path.abspath(__file__))

    def _run(cmd):
        try:
            out = subprocess.run(cmd, cwd=repo_dir, capture_output=True,
                                 text=True, timeout=10)
        except (OSError, subprocess.SubprocessError):
            return None
        return out.stdout.strip() if out.returncode == 0 else None

    commit = _run(['git', 'rev-parse', 'HEAD'])
    status = _run(['git', 'status', '--porcelain'])
    return {
        'git_commit': commit,
        'git_dirty': None if status is None else bool(status),
        'git_uncommitted': status.splitlines() if status else [],
    }


def weights_fingerprint(model) -> str:
    """SHA-256 over the model's weight values.

    Two bitwise reproducible runs produce the same digest, which makes
    reproducibility checkable by comparing one field of two provenance files.
    Hashing the .h5 file instead would be unreliable, because HDF5 byte
    layout can differ even when the stored arrays are identical.
    """
    import numpy as np
    digest = hashlib.sha256()
    for weight in model.get_weights():
        arr = np.ascontiguousarray(weight)
        digest.update(str(arr.dtype).encode())
        digest.update(str(arr.shape).encode())
        digest.update(arr.tobytes())
    return digest.hexdigest()


def unique_model_dir(parent: str, stem: str) -> str:
    """Create a model folder, refusing to merge into an existing one.

    VAE_bilayer.save() calls os.makedirs(exist_ok=True), so a name collision
    silently merges two runs: config.json and the weights are replaced, while
    stale files from the earlier run can survive in saved_model/. The result
    loads without error and contains parts of two models.

    Creating the directory here with exist_ok=False makes the collision
    visible and resolves it by suffixing. os.mkdir is atomic on both POSIX
    and NTFS, so this is also correct when several workers run concurrently.
    """
    os.makedirs(parent, exist_ok=True)
    candidate, suffix = os.path.join(parent, stem), 1
    while True:
        try:
            os.makedirs(candidate, exist_ok=False)
            if suffix > 1:
                print(f"NOTE: {stem} already existed; this run was written to "
                      f"{os.path.basename(candidate)} rather than merged into it.")
            return candidate
        except FileExistsError:
            suffix += 1
            candidate = os.path.join(parent, f"{stem}-{suffix}")


class LegacyCheckpointError(RuntimeError):
    """Raised when a checkpoint still identifies completed runs by index."""


def _normalise_combo(combo):
    """Return one JSON-loaded parameter combination as an immutable tuple."""
    if not isinstance(combo, (list, tuple)) or len(combo) != PARAMETER_COUNT:
        raise ValueError(
            f"Each parameter combination must contain exactly {PARAMETER_COUNT} values; got {combo!r}."
        )
    return tuple(combo)


def load_checkpoint(savefld):
    """Load completed parameter tuples from the checkpoint.

    Version-1 checkpoints stored integer grid indices. They are deliberately
    rejected because an index cannot be mapped safely after the grid changes.
    """
    checkpoint_path = get_checkpoint_path(savefld)
    if not os.path.exists(checkpoint_path):
        return {'version': CHECKPOINT_VERSION, 'completed': []}

    with open(checkpoint_path, 'r') as f:
        checkpoint = json.load(f)

    completed = checkpoint.get('completed', [])
    if any(isinstance(item, int) for item in completed):
        raise LegacyCheckpointError(
            "Legacy checkpoint detected: completed runs are stored as grid indices. "
            "Those indices cannot be converted safely if the parameter grid has changed. "
            "Back up or remove .checkpoint.json and restart the search (use --restart-from 1 "
            "for a full rerun)."
        )

    return {
        'version': CHECKPOINT_VERSION,
        'completed': [_normalise_combo(combo) for combo in completed],
    }


def save_checkpoint(savefld, completed_list):
    """Persist completed parameter tuples as JSON arrays."""
    checkpoint_path = get_checkpoint_path(savefld)
    payload = {
        'version': CHECKPOINT_VERSION,
        'completed': [list(_normalise_combo(combo)) for combo in completed_list],
    }
    with open(checkpoint_path, 'w') as f:
        json.dump(payload, f, indent=2)


def _completion_mask(multi_list, completed_combos):
    """Match completed tuples to the current grid one occurrence at a time.

    This is intentionally multiset-based: if the grid contains the same
    parameter tuple more than once, one checkpoint entry completes only one
    occurrence.
    """
    remaining_counts = Counter(completed_combos)
    mask = []
    for combo in multi_list:
        combo = tuple(combo)
        if remaining_counts[combo] > 0:
            mask.append(True)
            remaining_counts[combo] -= 1
        else:
            mask.append(False)
    return mask


def build_param_grid(args):
    """Build or load the hyperparameter grid using JSON only.

    The final tuple element is an explicit replicate index. By default each
    latent dimension is trained once, while selected latent dimensions can be
    assigned additional stochastic replicates through ``replicates_by_latent_dim``.
    """
    if args.params_file is None:
        lossfunc = ["cos+KL"]
        # alternatives: "mean_squared_error", "cosine_similarity", "MSE+KL", "cos+KL"
        hidden_layer = [1024*3]
        hidden_2 = [512]
        latent_dim = [96]
        reconstruct_loss_weight_exp = [0,1,2,3,4,5,6]

        # Number of independent stochastic training runs for selected latent dimensions.
        # Any latent dimension not listed here is trained once.
        replicates_by_latent_dim = {
            _: 5 for _ in latent_dim
        }

        multi_list = []
        for L, h1, h2, ld, rwexp in product(
            lossfunc,
            hidden_layer,
            hidden_2,
            latent_dim,
            reconstruct_loss_weight_exp,
        ):
            n_replicates = replicates_by_latent_dim.get(ld, 1)
            for replicate_idx in range(n_replicates):
                multi_list.append((L, h1, h2, ld, rwexp, replicate_idx))
        return multi_list

    with open(args.params_file, 'r') as f:
        data = json.load(f)

    if not isinstance(data, list):
        raise ValueError("Parameter JSON must contain a list of parameter combinations.")
    return [_normalise_combo(combo) for combo in data]


def run_driver(args):
    """Top-level orchestrator: no TensorFlow import here on purpose."""
    if not os.path.exists(args.savefld):
        os.mkdir(args.savefld)

    multi_list = build_param_grid(args)
    with open(os.path.join(args.savefld, 'grid-search-params.json'), 'w') as f:
        json.dump([list(combo) for combo in multi_list], f, indent=2)

    try:
        checkpoint = load_checkpoint(args.savefld)
    except LegacyCheckpointError as exc:
        if args.restart_from == 1:
            print(
                "Legacy index-based checkpoint detected. Because --restart-from 1 was requested, "
                "the checkpoint will be reset and the full grid will be rerun using tuple-based tracking."
            )
            save_checkpoint(args.savefld, [])
            checkpoint = {'version': CHECKPOINT_VERSION, 'completed': []}
        else:
            print(f"Checkpoint error: {exc}")
            sys.exit(1)
    except (ValueError, json.JSONDecodeError) as exc:
        print(f"Checkpoint error: {exc}")
        sys.exit(1)

    completed_combos = checkpoint['completed']
    completion_mask = _completion_mask(multi_list, completed_combos)

    if args.restart_from is not None:
        if args.restart_from < 1 or args.restart_from > len(multi_list):
            print(f"--restart-from must be between 1 and {len(multi_list)} (inclusive).")
            sys.exit(1)

        start_idx = args.restart_from - 1

        # Preserve only genuinely completed occurrences before the requested
        # restart point. Every current-grid occurrence from start_idx onward is
        # deliberately removed from the checkpoint and rerun.
        completed_combos = [
            multi_list[i]
            for i in range(start_idx)
            if completion_mask[i]
        ]
        save_checkpoint(args.savefld, completed_combos)
        remaining = list(range(start_idx, len(multi_list)))

        incomplete_prefix = [i + 1 for i in range(start_idx) if not completion_mask[i]]
        if incomplete_prefix:
            print(
                "WARNING: the checkpoint does not contain completed parameter tuples for "
                f"{len(incomplete_prefix)} iteration(s) before --restart-from: {incomplete_prefix}. "
                "They will not be run in this invocation."
            )

        print(
            f"Grid search: {len(multi_list)} total combos. "
            f"Restarting from iteration {args.restart_from}; "
            f"{len(remaining)} iterations will be run."
        )
    else:
        remaining = [i for i, is_done in enumerate(completion_mask) if not is_done]
        matched_completed = sum(completion_mask)
        print(
            f"Grid search: {len(multi_list)} total combos, "
            f"{matched_completed} already done, {len(remaining)} remaining."
        )

    this_script = os.path.abspath(__file__)

    # Environment for every child. PYTHONHASHSEED is read by the interpreter at
    # startup, so it cannot be fixed from inside a running process and has to be
    # passed to the subprocess here. CUDA_VISIBLE_DEVICES is read when
    # TensorFlow initialises its device list, so setting it in the child's
    # environment is more reliable than assigning to os.environ after import.
    child_env = os.environ.copy()
    child_env['PYTHONHASHSEED'] = '0'
    child_env['CUDA_VISIBLE_DEVICES'] = args.gpu
    child_env['TF_DETERMINISTIC_OPS'] = '1'      # only needed for TensorFlow < 2.9
    child_env['TF_CUDNN_DETERMINISTIC'] = '1'

    for combo_idx in remaining:
        combo_preview = tuple(multi_list[combo_idx])
        print(f"\n=== Launching subprocess for combo {combo_idx+1} of {len(multi_list)} "
              f"(seed {combo_seed(combo_preview)}) ===")
        cmd = [
            sys.executable, this_script,
            args.train,
            '-f', args.savefld,
            '-e', str(args.patience),
            '-b', str(args.batch_size),
            '-l', str(args.learning_rate),
            '-d', str(args.min_delta),
            '--gpu', args.gpu,
            '--worker-idx', str(combo_idx),
        ]
        if args.nondeterministic_ops:
            cmd.append('--nondeterministic-ops')
        cmd.extend(['--free-bits', str(args.free_bits)])
        if args.params_file is not None:
            cmd.extend(['-p', args.params_file])

        combo = tuple(multi_list[combo_idx])
        before_count = Counter(load_checkpoint(args.savefld)['completed'])[combo]
        result = subprocess.run(cmd, env=child_env)

        if result.returncode != 0:
            print(f"\n✗ Subprocess for combo {combo_idx+1} exited with code {result.returncode}. "
                  f"Stopping grid search. Re-run this same command to resume.")
            sys.exit(result.returncode)

        # Re-read checkpoint written by the child (belt-and-braces in case
        # combos are re-ordered/skipped externally).
        checkpoint = load_checkpoint(args.savefld)
        completed_combos = checkpoint['completed']
        after_count = Counter(completed_combos)[combo]
        if after_count != before_count + 1:
            print(
                f"\n✗ Worker exited cleanly but the completed count for combo {combo_idx+1} "
                "did not increase by exactly one. Stopping."
            )
            sys.exit(1)

    final_checkpoint = load_checkpoint(args.savefld)
    final_mask = _completion_mask(multi_list, final_checkpoint['completed'])
    still_missing = [i + 1 for i, is_done in enumerate(final_mask) if not is_done]
    if still_missing:
        print(
            f"\nRequested iterations completed, but {len(still_missing)} current-grid "
            f"iteration(s) remain incomplete: {still_missing}"
        )
    else:
        print("\nAll parameter combinations have been completed!")


def run_worker(args):
    """Train exactly ONE combo (args.worker_idx) in this process, then exit.

    All TensorFlow / heavy imports are deferred to here so the driver process
    never loads them.
    """
    # --- Environment first: these must be set BEFORE TensorFlow is imported.
    # TensorFlow reads CUDA_VISIBLE_DEVICES when it builds its device list, so
    # assigning it after the import (as this script previously did) has no
    # effect. The driver also sets these in the child environment; the
    # assignments here cover direct invocation of worker mode.
    os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu
    os.environ.setdefault('TF_DETERMINISTIC_OPS', '1')
    os.environ.setdefault('TF_CUDNN_DETERMINISTIC', '1')

    if os.environ.get('PYTHONHASHSEED') != '0':
        print("WARNING: PYTHONHASHSEED is not '0'. Python reads it before the "
              "interpreter starts, so it cannot be set from inside this process. "
              "Launch through the driver, which sets it for every subprocess.")

    import gc
    import socket
    import inspect
    import names_generator
    import numpy as np
    import pandas as pd
    import tensorflow as tf
    from tensorflow.keras import optimizers
    from vae import VAE_bilayer
    from datetime import datetime

    # Memory growth must be configured before any GPU is initialised.
    for gpu in tf.config.list_physical_devices('GPU'):
        tf.config.experimental.set_memory_growth(gpu, True)

    patience = args.patience
    multi_list = build_param_grid(args)
    combo_idx = args.worker_idx
    if combo_idx < 0 or combo_idx >= len(multi_list):
        print(f"Internal worker index {combo_idx} out of range (0..{len(multi_list)-1})")
        sys.exit(1)

    # --- Seed before anything consumes random state: the split, the weight
    # initialisation, the dropout masks and the reparameterisation noise.
    combo = _normalise_combo(multi_list[combo_idx])
    seed = combo_seed(combo)
    set_seed(seed, deterministic_ops=not args.nondeterministic_ops)
    print(f"Seed {seed}, derived from parameter tuple {list(combo)}")

    # 1. Load data
    Input = pd.read_pickle(args.train)

    n_missing = Input.isna().sum().sum()
    if n_missing:
        msg = f"+  Dataset has {n_missing} missing entries ({n_missing/Input.size:.1%} of all values).  +"
        print("+" * len(msg))
        print(msg)
        print("+" * len(msg))

    print("Input shape:", Input.shape)
    Train, Test = Input.copy(deep=True), Input.copy(deep=True)

    original_path = os.path.join(args.savefld, 'original.pkl.gz')
    if not os.path.exists(original_path):
        Input.to_pickle(original_path, compression='gzip')

    L, h1, h2, ld, rwexp, replicate_idx = combo
    print(
        f"\n\nTraining VAE combo {combo_idx+1} of {len(multi_list)}: "
        f"loss={L}, hidden1={h1}, hidden2={h2}, latent_dim={ld}, "
        f"rec_weight_exp={rwexp}, replicate={replicate_idx+1}\n"
    )

    print(f"free_bits = {args.free_bits}")

    vae = VAE_bilayer(
        original_dim=Train.shape[1],
        hidden_dim1=h1,
        hidden_dim2=h2,
        latent_dim=ld,
        loss_type=L,
        rec_weight=10 ** rwexp,
        free_bits=args.free_bits,
    )

    vae.compile(optimizer=optimizers.Adam(learning_rate=args.learning_rate))

    batch_size = args.batch_size
    train_np = Train.fillna(0).values
    test_np = Test.fillna(0).values
    # train_ds = tf.data.Dataset.from_tensor_slices(train_np).batch(batch_size)
    train_ds = tf.data.Dataset.from_tensor_slices(train_np).shuffle(
            buffer_size=len(train_np),
            seed=seed,
            reshuffle_each_iteration=True,
        ).batch(batch_size)
    val_ds = tf.data.Dataset.from_tensor_slices(test_np).batch(batch_size)
    del train_np, test_np

    early_stopping = tf.keras.callbacks.EarlyStopping(
        monitor='val_loss', mode='min',
        patience=patience, verbose=1, min_delta=args.min_delta,
        restore_best_weights=True
    )
    history = vae.fit(
        train_ds,
        epochs=999,
        batch_size=None,
        validation_data=val_ds,
        callbacks=[early_stopping]
    )

    results = vae.evaluate(val_ds, return_dict=True)
    print("\nTest results:\n" + "\n".join(f"  {k:<22}{v:.6g}" for k, v in results.items()) + "\n")
    print(f"Parameters: {L}, {h1}, {h2}, {ld}, {rwexp}, replicate={replicate_idx+1}")

    vae.build(input_shape=(None, Train.shape[1]))
    # The generated name is now a function of the combo, since names_generator
    # draws from the seeded random state. That makes it an identifier rather
    # than decoration, so the seed is included alongside it and the directory
    # is created with a guard instead of a merge.
    stem = '{ts}-s{seed}-{name}'.format(
        ts=datetime.now().strftime('%Y%m%d-%H%M'),
        seed=seed,
        name=names_generator.generate_name(),
    )
    savefld = unique_model_dir(args.savefld, stem)
    vae.save(savefld)

    # --- Provenance: everything needed to identify why two runs differ ---
    try:
        gpu_details = [
            tf.config.experimental.get_device_details(d).get('device_name', d.name)
            for d in tf.config.list_physical_devices('GPU')
        ]
    except Exception:
        gpu_details = [d.name for d in tf.config.list_physical_devices('GPU')]

    provenance = {
        'seed': seed,
        'seed_source': 'sha256 of the parameter tuple',
        'model_folder': os.path.basename(savefld),
        'parameter_tuple': list(combo),
        'parameters': {
            'loss_type': L,
            'hidden_dim1': h1,
            'hidden_dim2': h2,
            'latent_dim': ld,
            'rec_weight_exp': rwexp,
            'rec_weight': 10 ** rwexp,
            'replicate_idx': replicate_idx,
        },
        'free_bits': args.free_bits,
        'input_file': os.path.abspath(args.train),
        'input_sha256': file_sha256(args.train),
        'input_shape': list(Input.shape),
        'train_shape': list(Train.shape),
        'n_missing': int(n_missing),
        'learning_rate': args.learning_rate,
        'batch_size': args.batch_size,
        'patience': args.patience,
        'min_delta': args.min_delta,
        'epochs_run': len(history.history['loss']),
        'best_epoch': getattr(early_stopping, 'best_epoch', None),
        'stopped_epoch': int(early_stopping.stopped_epoch),
        'restore_best_weights': True,
        'final_metrics': {k: float(v) for k, v in results.items()},
        'weights_sha256': weights_fingerprint(vae),
        'deterministic_ops': not args.nondeterministic_ops,
        'environment': {
            'PYTHONHASHSEED': os.environ.get('PYTHONHASHSEED'),
            'CUDA_VISIBLE_DEVICES': os.environ.get('CUDA_VISIBLE_DEVICES'),
            'TF_DETERMINISTIC_OPS': os.environ.get('TF_DETERMINISTIC_OPS'),
            'TF_CUDNN_DETERMINISTIC': os.environ.get('TF_CUDNN_DETERMINISTIC'),
        },
        'hardware': {
            'hostname': socket.gethostname(),
            'gpu_devices': gpu_details,
        },
        'code': {
            **git_provenance(),
            'script': os.path.abspath(__file__),
            'script_sha256': file_sha256(os.path.abspath(__file__)),
            'vae_module': os.path.abspath(inspect.getfile(VAE_bilayer)),
            'vae_module_sha256': file_sha256(inspect.getfile(VAE_bilayer)),
        },
        'versions': {
            'python': sys.version.split()[0],
            'tensorflow': tf.__version__,
            'keras': tf.keras.__version__,
            'numpy': np.__version__,
            'pandas': pd.__version__,
        },
        'command': ' '.join(sys.argv),
        'timestamp': datetime.now().isoformat(timespec='seconds'),
    }
    with open(os.path.join(savefld, 'provenance.json'), 'w') as fh:
        json.dump(provenance, fh, indent=2)
    print(f"Weights fingerprint: {provenance['weights_sha256']}")

    del vae, train_ds, val_ds, results
    tf.keras.backend.clear_session()
    _ = gc.collect()

    checkpoint = load_checkpoint(args.savefld)
    completed_combos = checkpoint['completed']
    completed_combos.append(tuple(multi_list[combo_idx]))
    save_checkpoint(args.savefld, completed_combos)
    print(f"\n✓ Completed combo {combo_idx+1}. Parameter tuple checkpointed.")
    sys.exit(0)


if __name__ == '__main__':
    args = parse_cli()
    if args.worker_idx is None:
        run_driver(args)
    else:
        run_worker(args)
