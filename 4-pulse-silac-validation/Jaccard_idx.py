import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from more_itertools import distinct_combinations

def pseudo_jaccard(vec1,vec2):
    z = []
    for x,y in zip(vec1,vec2):
        if max(x,y) > 0:
            z.append(min(x,y)/max(x,y))
    return np.mean(z)

def jaccard_2(vec1, vec2):
    intersection = np.sum(vec1 & vec2)
    union = np.sum(vec1 | vec2)
    size1 = np.sum(vec1)
    size2 = np.sum(vec2)
    universe_size = len(vec1)

    jaccard_similarity = intersection / union if union > 0 else 0
    pvalue = hypergeom.sf(intersection - 1, universe_size, size1, size2)
    return jaccard_similarity,pvalue

    
    

def jaccard_similarity_binary(vec1, vec2):
    """Calculate Jaccard similarity between two binary vectors."""
    intersection = np.sum(vec1 & vec2)
    union = np.sum(vec1 | vec2)
    return intersection / union if union > 0 else 0

def jaccard_pvalue(vec1, vec2, universe_size):
    """
    Calculate p-value for Jaccard similarity using hypergeometric test.
    
    Parameters:
    -----------
    vec1, vec2 : array-like
        Binary vectors (0s and 1s)
    universe_size : int
        Total number of features (length of vectors)
    
    Returns:
    --------
    float : p-value
    """
    intersection = np.sum(vec1 & vec2)
    size1 = np.sum(vec1)
    size2 = np.sum(vec2)
    
    # Hypergeometric test
    pvalue = hypergeom.sf(intersection - 1, universe_size, size1, size2)
    
    return pvalue

def calculate_jaccard_matrix(data, pvalue_threshold=0.05, return_dataframe=True):
    """
    Calculate pairwise Jaccard similarities and p-values for a binary matrix.
    
    Parameters:
    -----------
    data : numpy array or pandas DataFrame
        Binary matrix where rows are observations and columns are features
    pvalue_threshold : float
        Threshold for filtering significant pairs
    return_dataframe : bool
        If True, returns pandas DataFrame; otherwise returns list of dicts
    
    Returns:
    --------
    DataFrame or list : Pairwise comparisons with Jaccard index and p-value
    """
    # Convert to numpy array if DataFrame
    if isinstance(data, pd.DataFrame):
        obs_names = data.index.tolist()
        matrix = data.values.astype(bool)
    else:
        matrix = np.array(data).astype(bool)
        obs_names = list(range(len(matrix)))
    
    n_obs, universe_size = matrix.shape
    
    # Calculate for all pairs
    results = []
    for i, j in distinct_combinations(range(n_obs), 2):
        jaccard = jaccard_similarity_binary(matrix[i], matrix[j])
        pval = jaccard_pvalue(matrix[i], matrix[j], universe_size)
        
        results.append({
            'obs1': obs_names[i],
            'obs2': obs_names[j],
            'jaccard': jaccard,
            'pvalue': pval,
            'significant': pval < pvalue_threshold
        })
    
    if return_dataframe:
        return pd.DataFrame(results)
    return results
    