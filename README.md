# MoDPA: Modification-Dependent Protein Associations

Post-translational modifications (PTMs) are key regulators of protein function and cellular processes; however, the overall principles of PTM co-regulation and crosstalk remain to be fully understood. 

To overcome the extreme sparsity and heterogeneity of PTM calls across experiments, MoDPA utilizes a variational autoencoder (VAE) to embed per-site detection profiles into a low-dimensional latent space that preserves covariation while denoising missing data. A PTM association network is constructed by correlating latent representations across experiments.

