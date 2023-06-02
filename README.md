# st_rle

`st_rle` is a Python package for calculating the size factors of a matrix using the Relative Log Expression (RLE) method. This method is robust to outliers and is particularly useful for normalizing count data, such as RNA-seq read counts.
This was inspired by the RLE normalization in scater:
https://pubmed.ncbi.nlm.nih.gov/28088763/
http://bioconductor.org/packages/scater
Which was inspired by the normalization in DESeq (at least from what I can tell):
https://pubmed.ncbi.nlm.nih.gov/20979621/

Make sure to cite those if you use this repo!

## Installation
*We're in alpha right now, so this is a placeholder*
You can install `st_rle` from PyPI:

```bash
python3 -m pip install st_rle
```

## Usage

Here is a basic use case example using a negative binomial distributed matrix:

```python
import numpy as np
from st_rle.rle import do_depth_normalization
from scipy.sparse import csc_matrix, csr_matrix, lil_matrix

# We'll create a random NB distributed 'full transcriptome' size dataset of 5000 samples
in_mat = np.random.negative_binomial(.1, .1, size=(20000,5000))

# And normalize it in different formats
in_mat_norm = do_depth_normalization(in_mat.astype(float))
in_mat_norm2 = do_depth_normalization(csc_matrix(in_mat.astype(float)))
in_mat_norm3 = do_depth_normalization(csr_matrix(in_mat.astype(float)))
```


## License

AGPL-v3


