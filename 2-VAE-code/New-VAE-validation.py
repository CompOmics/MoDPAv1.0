#!/usr/bin/env python
# coding: utf-8
import pandas as pd
from tensorflow import keras
import matplotlib.pyplot as plt
import seaborn as sns
import os, json, argparse
from tqdm import tqdm

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('folder', type=str, help="Path to models folder")
    return p.parse_args() 

# In[2]:
def evaluate_reconstruction(fld):
    x = pd.read_pickle(os.path.join(fld,'..','original.pkl.gz'))
    y = pd.read_pickle(os.path.join(fld,'reconstruction.pkl.gz'))
    cs = keras.metrics.CosineSimilarity(name="cosine_similarity", dtype=None, axis=1)
    cs.update_state(x,y)
    return cs.result().numpy()

def read_model_params(fld):
    cfg_path = os.path.join(fld, "config.json")
    with open(cfg_path, "r") as f:
        config = json.load(f)
    return config


# In[3]:
args = parse_cli()
models = []
folders = [_ for _ in os.scandir(args.folder) if os.path.isdir(_.path)]
for m in tqdm(folders):
    x = read_model_params(m)
    x = pd.DataFrame.from_dict(x, orient='index').T
    x['model_path'] = m.path
    x['model'] = m.name
    x['cosine_similarity'] = evaluate_reconstruction(m.path)
    x.set_index('model', inplace=True)
    models.append(x)
models = pd.concat(models)


# In[3]:
sns.catplot(data=models, kind='box', 
            y='cosine_similarity', x='hidden_dim2', hue='latent_dim', col='hidden_dim1',
            # notch=True,
            dodge=True, palette='pastel'
           )
plt.savefig(os.path.join(args.folder,'boxplot.png'), dpi=300, bbox_inches='tight')


# In[8]:
sns.catplot(data=models, kind='point', 
            y='cosine_similarity', x='hidden_dim2', hue='latent_dim', col='hidden_dim1',
            # linestyle='--',
            dodge=True, palette='pastel'
           )
plt.savefig(os.path.join(args.folder,'pointplot.png'), dpi=300, bbox_inches='tight')


# In[17]:
models.sort_values('cosine_similarity', ascending=False, inplace=True)
models.to_csv(os.path.join(args.folder,'models_summary.csv.gz'), compression='gzip')
print(models.iloc[0,:])
