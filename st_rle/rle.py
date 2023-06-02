import time
import numpy as np
from scipy.sparse import issparse
from scipy.sparse import csc_matrix, csr_matrix
from joblib import Parallel, delayed

def gm_mean(x):
    """
    Calculates the geometric mean of a numpy array after removing non-positive entries.
    
    Parameters
    ----------
    x : numpy.ndarray
        Input array.

    Returns
    -------
    float
        Geometric mean of positive elements in x. Returns 0 if no positive elements exist.
    """
    x = x[x > 0]
    if len(x) == 0:
        return 0
    elif len(x) == 1:
        return x[0]
    else:
        return np.exp(np.sum(np.log(x)) / len(x))

def calc_size_factor(data, i, gm, keep_rows, sparse_convert):
    """
    Calculates the size factor for a given column in the data.
    
    Parameters
    ----------
    data : scipy.sparse.csc_matrix or numpy.ndarray
        Input data.
    i : int
        Index of the column to calculate size factor.
    gm : numpy.ndarray
        Array of geometric means.
    keep_rows : numpy.ndarray
        Boolean array indicating which rows to keep.
    sparse_convert : bool
        Flag indicating if the input data is sparse.
        
    Returns
    -------
    float
        Calculated size factor.
    """
    if sparse_convert:
        data_subset = data[:,i].toarray().flatten()
    else:
        data_subset = data[:,i]
    data_subset = data_subset[keep_rows]
    active_subset_idxs = data_subset > 0
    data_subset = data_subset[active_subset_idxs]
    if np.sum(active_subset_idxs) == 0:
        return 0.0
    else:
        return np.median(data_subset / gm[active_subset_idxs])


def calcFactorRLE(data):
    """
    Calculates size factors for all columns in the data using RLE (Relative Log Expression) method.
    
    Parameters
    ----------
    data : scipy.sparse.csc_matrix or numpy.ndarray
        Input data.
        
    Returns
    -------
    numpy.ndarray
        Array of calculated size factors.
    """
    sparse_convert = issparse(data)
    start_gm = time.time()
    if sparse_convert:
        # faster to convert between sparse formats
        data = csr_matrix(data)
    gm = np.array([gm_mean(data[i,:].toarray().flatten() if sparse_convert else data[i,:]) for i in range(data.shape[0])])
    print("gm time:",time.time()-start_gm)
    keep_rows = gm > 0
    gm = gm[keep_rows]
    start_size=time.time()
    if sparse_convert:
        data = csc_matrix(data)
    size_factor = Parallel(n_jobs=-1)(delayed(calc_size_factor)(data, i, gm, keep_rows, sparse_convert) for i in range(data.shape[1]))
    print("size factor time:",time.time()-start_size)
    return np.array(size_factor)



def do_depth_normalization(data, size_factors = None, bin_size=7000):
    """
    Performs depth normalization on the input data.
    
    Parameters
    ----------
    data : scipy.sparse.csc_matrix or numpy.ndarray
        Input data.
    size_factors : numpy.ndarray, optional
        Precomputed size factors. If None, size factors are calculated (this is recommended).
    bin_size : int, optional
        Number of columns to be processed at a time during normalization.
        
    Returns
    -------
    scipy.sparse.csc_matrix or numpy.ndarray
        Depth-normalized data.
    """
    if issparse(data):
        # log the original type
        orig_format = type(data)
    if size_factors is None:
        size_factors = calcFactorRLE(data)
    norm_start = time.time()
    ## make sure we don't explode anything with zeros!
    size_factors[size_factors==0]=1
    if issparse(data):
        data=csc_matrix(data)
    num_bins = int(np.ceil(size_factors.shape[0] / bin_size))
    for b in range(num_bins):
        start_idx = b*bin_size
        end_idx = min(start_idx + bin_size,data.shape[1])
        print("\tnormalizing cols:", start_idx,"to",end_idx)
        if issparse(data):
            sub_matrix = data[:, start_idx:end_idx].toarray()
        else:
            sub_matrix = data[:, start_idx:end_idx]
        sub_size_factors = size_factors[start_idx:end_idx]
        # Normalize the submatrix
        sub_matrix /= sub_size_factors
        if issparse(data):
            data[:, start_idx:end_idx] = csc_matrix(sub_matrix)
        else:
            data[:, start_idx:end_idx] = sub_matrix
    if issparse(data):
        data = orig_format(data)
    print("normalization time:",time.time()-norm_start)
    return data







