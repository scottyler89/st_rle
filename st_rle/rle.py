import os
import time
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.sparse import issparse

obj_save_path = "/sc/arion/projects/schade01a/scott/anticor/big/1mil_pbmc/"
adata_path = os.path.join(obj_save_path,"adata_obj3.h5ad")
##############################################################################

def gm_mean(x):
    keep_idxs = np.where(x > 0)[0]
    if len(keep_idxs)==0:
        return(0)
    elif len(keep_idxs)==1:
        return(x[keep_idxs])
    else:
        x=x[keep_idxs]
        gm=np.exp(np.sum(np.log(x)) / x.shape[0])
        return(gm)

def new_calcFactorRLE(data):
    #   Scale factors as in Anders et al (2010)
    #   Mark Robinson
    #   Created 16 Aug 2010
    # gm <- exp(rowMeans(log(data)))
    # return(apply(data, 2, function(u) median((u/gm)[gm > 0])))
    sparse_convert = False
    if issparse(data):
        sparse_convert=True
    gm = np.zeros(data.shape[0])
    for i in range(data.shape[0]):
        if sparse_convert:
            gm[i] = gm_mean(data[i,:].toarray().flatten())
        else:
            gm[i] = gm_mean(data[i,:])
    keep_rows = np.where(gm>0)[0]
    gm = gm[keep_rows]
    size_factor = np.zeros(data.shape[1])
    for i in range(data.shape[1]):
        if sparse_convert:
            data_subset = data[:,i].toarray().flatten()
        else:
            data_subset = data[:,i]
        #print("1")
        data_subset = data_subset[keep_rows]
        #print("2")
        active_subset_idxs = np.where(data_subset>0)[0]
        #print("3")
        #print(active_subset_idxs.shape)
        data_subset = data_subset[active_subset_idxs]
        #print("4")
        #print(data_subset)
        size_factor[i] = np.median(data_subset/gm[active_subset_idxs])
    return(size_factor)

def do_depth_normalization(data, size_factors = None):
    if size_factors is None:
        size_factors = new_calcFactorRLE(data)
    for i in range(data.shape[1]):
        if size_factors[i]<=0:
            pass
        else:
            data[:,i]=data[:,i]/size_factors[i]
    return(data)






############################################################################


# import numpy as np
# from scipy.sparse import issparse
# from joblib import Parallel, delayed
# import cython

# @cython.boundscheck(False)
# cdef double gm_mean_cython(double[:] x):
#     cdef int len_keep_idxs = np.count_nonzero(x > 0)
#     if len_keep_idxs == 0:
#         return 0.0
#     elif len_keep_idxs == 1:
#         return np.max(x)
#     else:
#         return np.exp(np.sum(np.log(x[x > 0])) / len_keep_idxs)


# def new_calcFactorRLE(data):
#     sparse_convert = issparse(data)
#     gm = np.zeros(data.shape[0])
#     for i in range(data.shape[0]):
#         if sparse_convert:
#             gm[i] = gm_mean_cython(data[i,:].toarray().flatten())
#         else:
#             gm[i] = gm_mean_cython(data[i,:])
#     keep_rows = gm > 0
#     gm = gm[keep_rows]
#     size_factor = Parallel(n_jobs=-1)(delayed(calc_size_factor)(data, i, gm, keep_rows, sparse_convert) for i in range(data.shape[1]))
#     return np.array(size_factor)


# @cython.boundscheck(False)
# cdef double calc_size_factor(data, int i, double[:] gm, keep_rows, sparse_convert):
#     if sparse_convert:
#         data_subset = data[:,i].toarray().flatten()
#     else:
#         data_subset = data[:,i]
#     data_subset = data_subset[keep_rows]
#     active_subset_idxs = data_subset > 0
#     data_subset = data_subset[active_subset_idxs]
#     if np.sum(active_subset_idxs) == 0:
#         return 0.0
#     else:
#         return np.median(data_subset / gm[active_subset_idxs])


# def do_depth_normalization(data, size_factors = None):
#     if size_factors is None:
#         size_factors = new_calcFactorRLE(data)
#     for i in range(data.shape[1]):
#         if size_factors[i] > 0:
#             data[:,i] /= size_factors[i]
#     return data





################################################################
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


def new_calcFactorRLE(data):
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
        size_factors = new_calcFactorRLE(data)
    norm_start = time.time()
    ## make sure we don't explode anything with zeros!
    size_factors[size_factors==0]=1
    if issparse(data):
        data=csc_matrix(data)
    num_bins = int(np.ceil(len(non_zero_factors) / bin_size))
    for b in range(num_bins):
        print("\t", b*bin_size)
        start_idx = b*bin_size
        end_idx = min(start_idx + bin_size,data.shape[1])
        sub_matrix = data[:, start_idx:end_idx].toarray()
        sub_size_factors = size_factors[start_idx:end_idx]
        # Normalize the submatrix
        sub_matrix /= sub_size_factors
        data[:, start_idx:end_idx] = csc_matrix(sub_matrix)
    if issparse(data):
        data = orig_format(data)
    print("normalization time:",time.time()-norm_start)
    return data


import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, lil_matrix
from count_split.count_split import multi_split
from matplotlib import pyplot as plt

# We'll simulate the 'full transcriptome' of 5000 samples
in_mat = np.random.negative_binomial(.1, .1, size=(20000,5000))

# And normalize it in different formats
in_mat_norm = do_depth_normalization_old(in_mat.astype(float))
in_mat_norm2 = do_depth_normalization(csc_matrix(in_mat.astype(float)))
in_mat_norm3 = do_depth_normalization(csr_matrix(in_mat.astype(float)))




##############################################################################
print("reading")
adata = sc.read(adata_path)

start = time.time()
adata.X=do_depth_normalization(adata.X.T).T
total_time = time.time()-start
print("total_time:",total_time)
print("logging")
sc.pp.log1p(adata)
print("saving")
adata.write(os.path.join(obj_save_path,"adata_log_rle_norm.h5ad"))

##
print("\n\n\nclustering")
sc.tl.leiden(adata)
print(adata.obs.head())
print(adata.obs.columns)

print("\n\n\numapping")
sc.pl.umap(adata)

print(adata.obs.head())
print(adata.obs.columns)


adata.write(os.path.join(obj_save_path,"adata_normed_clustered.h5ad"))

