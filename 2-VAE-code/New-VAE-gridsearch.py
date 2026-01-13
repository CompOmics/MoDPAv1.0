#!/usr/bin/env python
# coding: utf-8
import os, gc, names_generator, re, argparse
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow.keras import optimizers
import matplotlib.pyplot as plt
import seaborn as sns
from vae import VAE_bilayer
from datetime import datetime, date

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('train', metavar='training data', type=str, help="Path to MoD")
    p.add_argument('--savefld', dest='savefld', type=str, default=date.today().isoformat(), help="Output folder (default: today's date)")
    p.add_argument('--gpu', dest='gpu', type=str, default='0', help="Which GPU to use (default: 0)")
    p.add_argument('--txt-params', dest='txt', type=str, default='0', help="Are there parameters saved in a txt file? If yes, provide the path (default: 0, meaning no). ")
    return p.parse_args() 

args = parse_cli()
os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu
os.mkdir(args.savefld)

if args.txt == '0':
    lossfunc     = ["cos+KL"]
    hidden_layer = [128, 256]
    hidden_2     = [ 32,  64, 128] 
    latent_dim   = [ 16,  32, 64]
    reconstruct_loss_weight_exp = [0] #list(range(4))
    iterations = range(1)
    multi = pd.MultiIndex.from_product([lossfunc, hidden_layer, hidden_2, latent_dim, reconstruct_loss_weight_exp, iterations])
    with open('grid-search-params.txt', 'w') as f:
        f.write(str(multi.tolist()))
else:
    with open(args.txt, 'r') as f:
        txt = f.read()
    multi = eval(txt)

# 1. Load data
Input = pd.read_pickle(args.train)
Input.fillna(0, inplace=True)
print(Input.values.max(), "|", Input.values.min(), "|", Input.shape)
Train, Test = Input.copy(deep=True), Input.copy(deep=True)
Input.to_pickle(os.path.join(args.savefld,'original.pkl.gz'), compression='gzip')

for L,h1,h2,ld,rwexp,_ in multi:
    # 2. Instantiate the VAE
    vae = VAE_bilayer(
        original_dim = Train.shape[1],
        hidden_dim1  = h1,
        hidden_dim2  = h2,
        latent_dim   = ld,
        loss_type    = L,   #alternatives: "mean_squared_error", "cosine_similarity", "RMSE+KL", "cos+KL",
        rec_weight   = 10**rwexp,
    )

    # 3. Compile with Adam
    vae.compile(optimizer = optimizers.Adam(learning_rate=1e-3))

    vae.fit(
        Train,
        epochs          = 999,
        batch_size      = 128,
        validation_data = (Test, None),  # targets are implicit
        callbacks=[tf.keras.callbacks.EarlyStopping(patience=8, verbose=1, min_delta=.25, restore_best_weights=True)]
    )

    results = vae.evaluate(
        Test,
        batch_size=256,
        # verbose=2
    )
    print(f"\nTest results:\n  Total loss = {results[0]:.4f}\n  Rec loss   = {results[1]:.4f}\n  KL loss    = {results[2]:.4f}\n")


    # # — After training —
    vae.build(input_shape=(None, Train.shape[1]))
    savefld = os.path.join(args.savefld, datetime.now().strftime('%Y%m%d-%H%M-')+names_generator.generate_name())
    vae.save(savefld) 


    # # — Get encodings —
    print('Get encodings...')
    ptms_list = [''.join(re.split(r'[][]',x)[:2]) for x in Input.index]
    exps_list = list(Input.columns)

    encoded       = np.array(vae.encode(Input.values))
    reconstructed = np.array(vae.decode(encoded[0]))

    print('Plotting...')
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(16,9))
    sns.heatmap(Input, vmin=0, vmax=1, center=0, 
                cmap='vlag', cbar=False, ax=ax[0])
    sns.heatmap(reconstructed, vmin=0, vmax=1, center=0, 
                cmap='vlag', cbar=False, ax=ax[1])
    fig.savefig(os.path.join(savefld,'reconstruction.png'), dpi=300, bbox_inches='tight')

    # Save x_encoded, x_reconstructed
    print('Saving latent space...')
    encoded = pd.DataFrame(encoded[0,:,:], index=ptms_list)
    encoded.to_pickle(os.path.join(savefld,'Latent-space.pkl.gz'), compression='gzip')
    print('Saving reconstruction...')
    reconstructed = pd.DataFrame(reconstructed, index=ptms_list)
    reconstructed.to_pickle(os.path.join(savefld,'reconstruction.pkl.gz'), compression='gzip')


    # — garbage collector —
    del vae
    tf.keras.backend.clear_session()
    _ = gc.collect()