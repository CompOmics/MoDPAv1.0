import os
import json

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

tf.keras.backend.set_floatx('float64')

# In TensorFlow 2.17 the functional alias keras.losses.mean_squared_error(y_true, y_pred)
# has been removed. Use keras.losses.mse instead.


class VAE_bilayer(keras.Model):
    """
    Two-layer Variational Autoencoder (VAE) in Keras.

    Args:
        original_dim (int): Dimension of the input vector.
        hidden_dim1  (int): Units in the first hidden layer.
        hidden_dim2  (int): Units in the second hidden layer.
        latent_dim   (int): Dimensionality of the latent code z.
        loss_type    (str): One of
                           - "mean_squared_error" (MSE only)
                           - "cosine_similarity"  (cos only)
                           - "MSE+KL"             (MSE + KL)
                           - "cos+KL"             (cos + KL)
                           Any other value raises at construction.
        rec_weight   (float): Weight on the reconstruction term. This is the
                      only knob on the reconstruction/KL balance; the KL term
                      enters the objective with weight 1.
        free_bits    (float): Nats per latent dimension below which a dimension
                      is not penalised. 0.0 disables. Prevents posterior
                      collapse once the KL term is strong enough relative to
                      reconstruction for collapse to occur, that is once
                      rec_weight is small enough.
        dropout_rate (float): Dropout fraction in the encoder.

    Note on rec_weight
    ------------------
    The objective is rec_weight * L_rec + KL. A separate KL weight (beta) would
    be redundant with rec_weight, because rec_weight * L_rec + beta * KL equals
    rec_weight * (L_rec + (beta / rec_weight) * KL): only the ratio moves the
    optimum, and Adam is scale-invariant up to its epsilon, so the overall
    factor has essentially no effect on training. rec_weight therefore carries
    the whole balance and is the only axis to sweep. A grid cell with
    rec_weight = w reproduces the old (rec_weight, beta) pair for any
    rec_weight / beta = w.

    One consequence of sweeping rec_weight rather than a KL weight: the total
    loss moves by orders of magnitude between grid cells, so any absolute
    threshold that reads the loss value is not comparable across the sweep. In
    particular, set EarlyStopping to monitor a quantity that does not carry
    rec_weight, such as rec_raw or a held-out reconstruction metric, rather than
    loss with a fixed min_delta.

    Note on the reconstruction and KL scales
    ----------------------------------------
    The KL term is the sum over latent dimensions of the per-dimension KL,
    averaged over the batch. This is what the ELBO requires for a factorised
    posterior. Averaging over both axes instead would give the ELBO's KL divided
    by latent_dim, so models at different latent_dim would be regularised by
    different amounts and a sensitivity analysis would confound dimensionality
    with regularisation strength.

    The MSE branch multiplies by original_dim, converting the per-feature mean
    returned by keras.losses.mse into the per-feature sum a likelihood requires.
    The cosine branch does not, because cosine similarity is already a
    per-sample scalar and there is nothing to convert.

    Note on likelihood interpretation
    ---------------------------------
    Both fixed-variance branches are Gaussian ELBOs, so rec_weight is an inverse
    variance rather than an arbitrary constant:
      MSE+KL  L_rec = E[sum_d (x_d - xhat_d)^2],      rec_weight = 1 / (2 sigma^2)
      cos+KL  L_rec = E[||x_n - xhat_n||^2] / 2,      rec_weight = 1 / sigma^2
    (the second uses ||x_n - xhat_n||^2 = 2 (1 - cos) on L2-normalised rows).
    Report sigma rather than rec_weight where possible.

    Because the KL term has weight 1, rec_weight is exactly that inverse
    variance and not an arbitrary constant scaled by a second parameter.

    """

    _ALLOWED_LOSSES = (
        "mean_squared_error",
        "cosine_similarity",
        "MSE+KL",
        "cos+KL",
    )

    def __init__(self,
                 original_dim: int,
                 hidden_dim1:  int,
                 hidden_dim2:  int,
                 latent_dim:   int,
                 loss_type:    str = "MSE+KL",
                 rec_weight:   float = 1.0,
                 free_bits:    float = 0.0,
                 dropout_rate: float = 0.25,
                 **kwargs):
        super().__init__(**kwargs)

        if loss_type not in self._ALLOWED_LOSSES:
            raise ValueError(
                f"Unsupported loss_type {loss_type!r}. Choose one of "
                f"{', '.join(self._ALLOWED_LOSSES)}. Unknown values are rejected "
                "here rather than dispatched to keras.losses.get(), which would "
                "silently return a per-feature mean and put the reconstruction "
                "term on a different scale from the branches above."
            )

        # Store config
        self.original_dim = original_dim
        self.loss_type    = loss_type
        self.rec_weight   = float(rec_weight)
        self.free_bits    = float(free_bits)
        self.dropout_rate = dropout_rate


        # ===== Encoder layers =====
        self.dense_e1   = layers.Dense(hidden_dim1, activation="relu", name="enc_dense1")
        self.dropout    = layers.Dropout(dropout_rate, name="enc_dropout")
        self.dense_e2   = layers.Dense(hidden_dim2, activation="relu", name="enc_dense2")
        self.dense_mu   = layers.Dense(latent_dim,  name="enc_mu")
        self.dense_logv = layers.Dense(latent_dim,  name="enc_logvar")

        # ===== Decoder layers =====
        # Sigmoid assumes inputs scaled to [0, 1]. The cosine branches are
        # insensitive to output magnitude, but the MSE branches
        # identify it, so an input range outside [0, 1] biases those two and no
        # value of rec_weight corrects it.
        self.dense_d1   = layers.Dense(hidden_dim2, activation="relu", name="dec_dense1")
        self.dense_d2   = layers.Dense(hidden_dim1, activation="relu", name="dec_dense2")
        self.output_dec = layers.Dense(original_dim, activation="sigmoid", name="dec_output")

        # ===== Loss & metric trackers =====
        self.total_loss_tracker = keras.metrics.Mean(name="loss")
        self.rec_loss_tracker   = keras.metrics.Mean(name="loss_reconstruction")
        # Unweighted reconstruction term. loss_reconstruction is scaled by
        # rec_weight, so it is not comparable between grid cells that use
        # different rec_weight. Select models on rec_raw, not on loss.
        self.rec_raw_tracker    = keras.metrics.Mean(name="rec_raw")
        self.kl_loss_tracker    = keras.metrics.Mean(name="loss_kl")
        # Unweighted, un-floored KL, for monitoring only. The weighted term can
        # look stable while the underlying KL collapses, so track both.
        self.kl_raw_tracker     = keras.metrics.Mean(name="kl_raw")

    @property
    def uses_kl(self) -> bool:
        return "KL" in self.loss_type

    @property
    def metrics(self):
        # Keras resets these trackers at the start of each epoch. Order matters:
        # evaluate() returns values in this order, so append new trackers rather
        # than inserting them.
        return [
            self.total_loss_tracker,
            # self.rec_loss_tracker,
            self.rec_raw_tracker,
            # self.kl_loss_tracker,
            self.kl_raw_tracker
            ]

    def encode(self, x, training=False):
        h = self.dense_e1(x)
        h = self.dropout(h, training=training)
        h = self.dense_e2(h)
        mu = self.dense_mu(h)
        logvar = self.dense_logv(h)
        return mu, logvar

    def reparameterize(self, mu, logvar):
        eps = tf.random.normal(shape=tf.shape(mu), dtype=tf.float64)
        return mu + tf.exp(0.5 * logvar) * eps

    def decode(self, z):
        h = self.dense_d1(z)
        h = self.dense_d2(h)
        return self.output_dec(h)

    def call(self, inputs, training=False):
        mu, logvar = self.encode(inputs, training=training)
        z = self.reparameterize(mu, logvar)
        return self.decode(z)

    # ---------- loss components ----------
    def _reconstruction_loss(self, x, reconstruction):
        """Return (weighted reconstruction term, unweighted reconstruction term)."""
        if self.loss_type in ("mean_squared_error", "MSE+KL"):
            # keras.losses.mse averages over features. Multiplying by
            # original_dim converts that mean into the sum over features
            # that a per-feature likelihood requires.
            rec = (tf.reduce_mean(keras.losses.mse(x, reconstruction))
                   * tf.cast(self.original_dim, x.dtype))

        elif self.loss_type in ("cosine_similarity", "cos+KL"):
            # keras.losses.cosine_similarity returns the NEGATIVE cosine, so
            # 1 + it is 0 for an identical reconstruction and 2 for an opposite
            # one. No original_dim factor: the cosine is already a per-sample
            # scalar, not a per-feature mean, so there is nothing to convert.
            rec = tf.reduce_mean(1.0 + keras.losses.cosine_similarity(x, reconstruction))

        return rec * self.rec_weight, rec

    def _kl_loss(self, mu, logvar):
        """Return (KL entering the loss, unweighted raw KL).

        The KL enters the objective with weight 1: the reconstruction/KL balance
        is set entirely by rec_weight. The two returned values differ only when
        free_bits is active.

        The per-dimension KL is averaged over the batch and summed over the
        latent dimension. Summing is what the ELBO requires for a factorised
        posterior; a mean over both axes would give the KL divided by
        latent_dim, making regularisation strength a function of the latent
        dimensionality being tested.
        """
        kl_per_dim = tf.reduce_mean(
            -0.5 * (1.0 + logvar - tf.square(mu) - tf.exp(logvar)), axis=0
        )                                                   # [latent_dim]
        kl_raw = tf.reduce_sum(kl_per_dim)

        if self.free_bits > 0.0:
            # Free bits: stop penalising a dimension once it is below the floor,
            # which prevents it collapsing to the prior and keeps the number of
            # active dimensions stable between runs. This puts a floor of
            # latent_dim * free_bits on the reported KL term, so compare
            # kl_raw rather than loss_kl when checking for collapse.
            kl_per_dim = tf.maximum(kl_per_dim, tf.cast(self.free_bits, kl_per_dim.dtype))

        return tf.reduce_sum(kl_per_dim), kl_raw

    # ---------- shared loss computation ----------
    def _compute_losses(self, data, training):
        """Single implementation shared by train_step and test_step.

        Keeping one copy matters because EarlyStopping selects on the test_step
        value, so a drift between two copies would silently change model
        selection.
        """
        x = data
        mu, logvar     = self.encode(x, training=training)
        z              = self.reparameterize(mu, logvar)
        reconstruction = self.decode(z)

        rec_loss, rec_raw = self._reconstruction_loss(x, reconstruction)
        kl_loss,  kl_raw  = self._kl_loss(mu, logvar)

        if self.uses_kl:
            total_loss = rec_loss + kl_loss
        else:
            # kl_loss and kl_raw are still tracked, for monitoring only.
            total_loss = rec_loss

        return total_loss, rec_loss, rec_raw, kl_loss, kl_raw

    def _update_trackers(self, total_loss, rec_loss, rec_raw, kl_loss, kl_raw):
        self.total_loss_tracker.update_state(total_loss)
        # self.rec_loss_tracker.update_state(rec_loss)
        self.rec_raw_tracker.update_state(rec_raw)
        # self.kl_loss_tracker.update_state(kl_loss)
        self.kl_raw_tracker.update_state(kl_raw)
        return {
            "loss":                self.total_loss_tracker.result(),
            # "loss_reconstruction": self.rec_loss_tracker.result(),
            "rec_raw":             self.rec_raw_tracker.result(),
            # "loss_kl":             self.kl_loss_tracker.result(),
            "kl_raw":              self.kl_raw_tracker.result(),
        }

    def train_step(self, data):
        with tf.GradientTape() as tape:
            losses = self._compute_losses(data, training=True)
        total_loss = losses[0]

        # Backpropagate. The None-gradient filter is defensive.
        grads = tape.gradient(total_loss, self.trainable_variables)
        pairs = [(g, v) for g, v in zip(grads, self.trainable_variables) if g is not None]
        self.optimizer.apply_gradients(pairs)

        return self._update_trackers(*losses)

    def test_step(self, data):
        # Same as train_step but without gradients and with dropout disabled.
        return self._update_trackers(*self._compute_losses(data, training=False))

    def kl_per_dimension(self, x, batch_size: int = 1024, threshold: float = 0.01):
        """Average KL from the prior for each latent dimension, in nats.

        Use this to detect posterior collapse and to choose free_bits. A
        dimension whose KL is near zero has mu near 0 and logvar near 0 for
        every sample: it carries no information about the input and contributes
        nothing to the latent distances used downstream.

        free_bits must sit below the KL that healthy dimensions reach, or every
        dimension is floored, the KL term becomes a constant and the
        regularisation is switched off entirely. Read the quantiles here before
        choosing it.

        Args:
            x: input matrix, shape [n_samples, original_dim].
            batch_size: encode in chunks to bound memory.
            threshold: nats above which a dimension counts as active.

        Returns:
            dict with the per-dimension KL, its quantiles, the total, and two
            collapse counts. Dropout is disabled and mu is used directly, so the
            result is deterministic.
        """
        import numpy as np
        x = np.asarray(x)
        n_samples = x.shape[0]
        kl_sum, mu_chunks = None, []

        for start in range(0, n_samples, batch_size):
            chunk = tf.convert_to_tensor(x[start:start + batch_size], dtype=tf.float64)
            mu, logvar = self.encode(chunk, training=False)
            kl = -0.5 * (1.0 + logvar - tf.square(mu) - tf.exp(logvar))
            batch_sum = tf.reduce_sum(kl, axis=0).numpy()
            kl_sum = batch_sum if kl_sum is None else kl_sum + batch_sum
            mu_chunks.append(mu.numpy())

        kl_per_dim = kl_sum / n_samples
        mu_all = np.concatenate(mu_chunks, axis=0)
        return {
            "latent_dim":    int(kl_per_dim.shape[0]),
            "kl_per_dim":    kl_per_dim,
            "kl_total":      float(kl_per_dim.sum()),
            "n_active":      int((kl_per_dim > threshold).sum()),
            # variance of mu across samples: an independent collapse indicator
            "n_informative": int((mu_all.var(axis=0) > 1e-3).sum()),
            "quantiles":     {q: float(np.quantile(kl_per_dim, q))
                              for q in (0.05, 0.10, 0.25, 0.50, 0.75, 0.95)},
        }

    def get_config(self):
        """Arguments needed to reconstruct the model.

        Keras uses this for serialization; save() also writes it to JSON.
        """
        return {
            "original_dim": self.original_dim,
            "hidden_dim1":  self.dense_e1.units,
            "hidden_dim2":  self.dense_e2.units,
            "latent_dim":   self.dense_mu.units,
            "loss_type":    self.loss_type,
            "rec_weight":   self.rec_weight,
            "free_bits":    self.free_bits,
            "dropout_rate": self.dropout_rate,
        }

    def save(self, folder, save_weights_only: bool = False):
        """Save config and weights. save_weights_only skips the SavedModel."""
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "config.json"), "w") as f:
            json.dump(self.get_config(), f, indent=2)
        self.save_weights(os.path.join(folder, "vae.weights.h5"))
        if not save_weights_only:
            tf.saved_model.save(self, os.path.join(folder, "saved_model"))

    @classmethod
    def load_vae(cls, folder, compile_kwargs=None):
        """Rebuild from saved config.json + weights.

        Args:
            folder: directory containing config.json and vae.weights.h5
            compile_kwargs: dict passed to model.compile()
        """
        cfg_path = os.path.join(folder, "config.json")
        with open(cfg_path, "r") as f:
            config = json.load(f)

        # Configs written before censored+KL was removed record an `lod` and a
        # loss_type this class no longer implements. Fail with a clear message
        # rather than a TypeError from **kwargs.
        if "lod" in config or config.get("loss_type") == "censored+KL":
            raise ValueError(
                f"{cfg_path} was written by a version that supported "
                'loss_type="censored+KL", which this class no longer implements. '
                "Retrain with a supported loss_type."
            )

        # Configs written before the KL weight was removed record a `beta`.
        # Only the ratio rec_weight / beta affects the optimum, so folding beta
        # into rec_weight reproduces the objective the weights were trained
        # under. The absolute loss scale changes by 1 / beta, so reported loss
        # values are not comparable with the original run's logs; rec_raw and
        # kl_raw are unaffected.
        if "beta" in config:
            beta = float(config.pop("beta"))
            uses_kl = "KL" in config.get("loss_type", "")
            if not uses_kl or beta == 1.0:
                pass  # beta had no effect on this model
            elif beta == 0.0:
                raise ValueError(
                    f"{cfg_path} records beta=0.0 with loss_type="
                    f"{config.get('loss_type')!r}, i.e. a KL term switched off "
                    "entirely. This class weights the KL by 1 and cannot express "
                    "that. Use a loss_type without KL instead."
                )
            else:
                rec_weight = float(config.get("rec_weight", 1.0)) / beta
                print(f"NOTE: {cfg_path} records beta={beta:g}; folding it into "
                      f"rec_weight, which becomes {rec_weight:g}. The objective is "
                      "unchanged up to an overall factor, but absolute loss values "
                      "will differ from the original run by 1/beta.")
                config["rec_weight"] = rec_weight

        config.setdefault("free_bits", 0.0)
        model = cls(**config)

        # Build variables with a real forward pass instead of model.build(), which
        # avoids Keras warnings on subclassed models without a custom build().
        dummy = tf.zeros((1, config["original_dim"]), dtype=tf.float64)
        _ = model(dummy, training=False)

        model.load_weights(os.path.join(folder, "vae.weights.h5"))

        if compile_kwargs is None:
            compile_kwargs = {"optimizer": "adam"}
        if compile_kwargs:
            model.compile(**compile_kwargs)

        return model
