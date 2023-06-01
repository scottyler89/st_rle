# st_rle

`st_rle` is a Python package for calculating the size factors of a matrix using the Robust Logarithmic Estimation (RLE) method. This method is robust to outliers and is particularly useful for normalizing count data, such as RNA-seq read counts.

## Installation

You can install `st_rle` from PyPI:

```bash
pip install st_rle
```

## Usage

Here is a basic use case example using a negative binomial distributed matrix:

```python
import numpy as np
from st_rle.do_depth_normalization import do_depth_normalization
from scipy.stats import nbinom

# Generate a negative binomial distributed matrix
np.random.seed(0)
n = 100  # number of genes
m = 20  # number of samples
data = nbinom.rvs(10, 0.5, size=(n, m)).astype(np.float64)

# Normalize the data
normalized_data = do_depth_normalization(data)
```

In this example, data is a 100x20 matrix representing the count data for 100 genes across 20 samples. The function new_calcFactorRLE(data) calculates the size factors for the data, and do_depth_normalization(data, size_factors) normalizes the data using these size factors.

## Usage

AGPL-v3


