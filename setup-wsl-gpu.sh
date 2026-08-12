#!/usr/bin/env bash
# ==============================================================================
# setup-wsl-gpu.sh — Set up and verify the MoDPA GPU environment on WSL
#
# Usage:
#   bash setup-wsl-gpu.sh [setup|check]
#
#   setup  — Create the conda environment from env-wsl-gpu.yml + install PTMmap
#   check  — Verify GPU is visible to TensorFlow
#
# Example:
#   bash setup-wsl-gpu.sh setup
#   bash setup-wsl-gpu.sh check
# ==============================================================================
set -eo pipefail

# ── Resolve the project root (where this script lives) ──────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Colours ──────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
info()  { echo -e "${GREEN}[INFO]${NC}  $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*"; }

# ── Conda helpers ─────────────────────────────────────────────────────────────
setup_conda() {
    if ! command -v conda &>/dev/null; then
        # Non-interactive shells skip .bashrc; try sourcing conda init directly
        local conda_sh=""
        for candidate in \
            "$HOME/miniconda3/etc/profile.d/conda.sh" \
            "$HOME/anaconda3/etc/profile.d/conda.sh" \
            "$HOME/miniforge3/etc/profile.d/conda.sh" \
            "$HOME/mambaforge/etc/profile.d/conda.sh" \
            "/opt/conda/etc/profile.d/conda.sh"; do
            if [ -f "$candidate" ]; then
                conda_sh="$candidate"
                break
            fi
        done
        if [ -z "$conda_sh" ]; then
            error "conda not found in PATH and no conda.sh init script detected."
            error "Install Miniconda/Anaconda in WSL first, or run this script from an interactive WSL shell."
            exit 1
        fi
        # shellcheck source=/dev/null
        source "$conda_sh"
    fi
    eval "$(conda shell.bash hook)"
}

env_exists() {
    conda env list | awk 'NR > 2 {print $1}' | grep -Fxq "MoDPA-env"
}

activate_modpa_env() {
    setup_conda
    if ! env_exists; then
        error "Environment 'MoDPA-env' does not exist. Run: bash setup-wsl-gpu.sh setup"
        exit 1
    fi
    if [ "${CONDA_DEFAULT_ENV:-}" != "MoDPA-env" ]; then
        conda activate MoDPA-env
    fi
}

# ── Helper: install conda activation hook for NVIDIA LD_LIBRARY_PATH ─────────
install_conda_hook() {
    local hook_dir="${CONDA_PREFIX}/etc/conda/activate.d"
    mkdir -p "$hook_dir"
    cat > "$hook_dir/nvidia-ld-path.sh" << 'HOOKEOF'
#!/usr/bin/env bash
# Auto-set LD_LIBRARY_PATH for tensorflow[and-cuda] NVIDIA pip packages.
_site="$(python -c 'import site; print(site.getsitepackages()[0])')"
_nvidia_dir="${_site}/nvidia"
if [ -d "$_nvidia_dir" ]; then
    for _d in "$_nvidia_dir"/*/lib; do
        [ -d "$_d" ] || continue
        case ":${LD_LIBRARY_PATH:-}:" in
            *":${_d}:"*) ;;
            *) export LD_LIBRARY_PATH="${_d}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" ;;
        esac
    done
fi
if [ -d "/usr/lib/wsl/lib" ]; then
    case ":${LD_LIBRARY_PATH:-}:" in
        *":/usr/lib/wsl/lib:"*) ;;
        *) export LD_LIBRARY_PATH="/usr/lib/wsl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" ;;
    esac
fi
unset _site _nvidia_dir _d
HOOKEOF
    chmod +x "$hook_dir/nvidia-ld-path.sh"
    info "Conda GPU hook installed to: $hook_dir/nvidia-ld-path.sh"
}

# ── Helper: prepend a dir to LD_LIBRARY_PATH only once ───────────────────────
prepend_ld_library_path_once() {
    local dir="$1"
    [ -d "$dir" ] || return 0
    case ":${LD_LIBRARY_PATH:-}:" in
        *":${dir}:"*) ;;
        *) export LD_LIBRARY_PATH="${dir}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" ;;
    esac
}

# ── Helper: set LD_LIBRARY_PATH for pip-installed NVIDIA libs ────────────────
setup_nvidia_ld_path() {
    # tensorflow[and-cuda] installs nvidia-* packages in site-packages.
    # Their shared libraries must be on LD_LIBRARY_PATH.
    local site_pkgs
    site_pkgs="$(python -c 'import site; print(site.getsitepackages()[0])')"
    local nvidia_dir="${site_pkgs}/nvidia"

    if [ -d "$nvidia_dir" ]; then
        local lib_dirs=""
        for d in "$nvidia_dir"/*/lib; do
            [ -d "$d" ] && lib_dirs="${lib_dirs:+$lib_dirs:}$d"
        done
        if [ -n "$lib_dirs" ]; then
            local old_ifs="$IFS"
            IFS=':'
            for d in $lib_dirs; do
                prepend_ld_library_path_once "$d"
            done
            IFS="$old_ifs"
            info "LD_LIBRARY_PATH updated with $(echo "$lib_dirs" | tr ':' '\n' | wc -l) NVIDIA lib paths"
        fi
    fi

    # Also add the CUDA stubs directory if present (needed for libcuda.so fallback)
    local cuda_compat="/usr/lib/wsl/lib"
    if [ -d "$cuda_compat" ]; then
        prepend_ld_library_path_once "$cuda_compat"
    fi
}

# ── Command: setup ───────────────────────────────────────────────────────────
do_setup() {
    info "=== Step 1: Checking NVIDIA driver in WSL ==="
    if ! command -v nvidia-smi &>/dev/null; then
        error "nvidia-smi not found!"
        echo "  Make sure you have the latest NVIDIA GPU driver installed on Windows."
        echo "  The driver passes through to WSL automatically. No CUDA toolkit"
        echo "  install is needed inside WSL — tensorflow[and-cuda] bundles it."
        echo "  Download from: https://www.nvidia.com/Download/index.aspx"
        exit 1
    fi
    nvidia-smi
    echo ""

    info "=== Step 2: Creating conda environment from env-wsl-gpu.yml ==="
    setup_conda
    if env_exists; then
        warn "Environment 'MoDPA-env' already exists. Updating it..."
        conda env update -f env-wsl-gpu.yml --prune
    else
        conda env create -f env-wsl-gpu.yml
    fi

    info "=== Step 3: Installing PTMmap wheel ==="
    activate_modpa_env

    if [ -f "PTMmap-0.1.3-py3-none-any.whl" ]; then
        pip install PTMmap-0.1.3-py3-none-any.whl
        info "PTMmap installed."
    else
        warn "PTMmap wheel not found at ${SCRIPT_DIR}/PTMmap-0.1.3-py3-none-any.whl"
        warn "You may need to install it manually."
    fi

    info "=== Step 4: Installing conda GPU activation hook ==="
    install_conda_hook

    info "=== Step 5: Verifying GPU visibility ==="
    setup_nvidia_ld_path
    do_check

    echo ""
    info "=== Setup complete! ==="
    echo ""
    echo "To activate the environment in future sessions:"
    echo "  conda activate MoDPA-env"
}

# ── Command: check ───────────────────────────────────────────────────────────
do_check() {
    info "Checking TensorFlow GPU support..."
    activate_modpa_env
    setup_nvidia_ld_path

    python - <<'PYEOF'
import sys

print(f"Python: {sys.executable}")
print(f"Python version: {sys.version}")

import tensorflow as tf
print(f"TensorFlow version: {tf.__version__}")

# Check build info
print(f"Built with CUDA: {tf.test.is_built_with_cuda()}")
print(f"Built with GPU:  {tf.test.is_built_with_gpu_support()}")

# List physical devices
gpus = tf.config.list_physical_devices('GPU')
cpus = tf.config.list_physical_devices('CPU')
print(f"CPUs found: {len(cpus)}")
print(f"GPUs found: {len(gpus)}")

if gpus:
    for i, gpu in enumerate(gpus):
        details = tf.config.experimental.get_device_details(gpu)
        cc = details.get("compute_capability")
        if cc:
            print(f"  GPU {i}: {gpu} (compute capability: {cc[0]}.{cc[1]})")
        else:
            print(f"  GPU {i}: {gpu}")

    # Quick computation test
    try:
        with tf.device('/GPU:0'):
            a = tf.random.normal([1000, 1000])
            b = tf.random.normal([1000, 1000])
            c = tf.matmul(a, b)
        print(f"\n✅ GPU computation successful! Result shape: {c.shape}")
    except Exception as e:
        msg = str(e)
        print("\n❌ GPU is visible, but TensorFlow failed to run a GPU kernel.")
        if "CUDA_ERROR_INVALID_PTX" in msg or "CUDA_ERROR_INVALID_HANDLE" in msg:
            print("\nUNSUPPORTED GPU IN THIS TF BUILD")
            print("TensorFlow can see your GPU, but this TensorFlow build does not include")
            print("compatible CUDA kernels/PTX for this GPU architecture.")
            print("Try updating NVIDIA drivers and/or installing a newer TensorFlow build.")
        print(f"\nOriginal error: {msg}")
        raise SystemExit(2)
else:
    print("\n❌ No GPU detected by TensorFlow!")
    print("\nTroubleshooting:")
    print("  1. Run 'nvidia-smi' to confirm the driver works in WSL")
    print("  2. Check LD_LIBRARY_PATH includes the nvidia pip package libs")
    print("  3. Ensure tensorflow[and-cuda] was installed (not plain tensorflow)")

    # Check if nvidia packages are installed
    import importlib
    for pkg in ['nvidia.cudnn', 'nvidia.cublas', 'nvidia.cuda_runtime']:
        try:
            importlib.import_module(pkg)
            print(f"  ✓ {pkg} is installed")
        except ImportError:
            print(f"  ✗ {pkg} is NOT installed — run: pip install tensorflow[and-cuda]==2.21.0")
            break
PYEOF
}

# ── Main ─────────────────────────────────────────────────────────────────────
case "${1:-help}" in
    setup) do_setup ;;
    check) do_check ;;
    *)
        echo "Usage: bash setup-wsl-gpu.sh [setup|check]"
        echo ""
        echo "  setup  — Create conda env + install dependencies + verify GPU"
        echo "  check  — Verify TensorFlow can see the GPU"
        echo ""
        echo "Examples:"
        echo "  bash setup-wsl-gpu.sh setup"
        echo "  bash setup-wsl-gpu.sh check"
        ;;
esac
