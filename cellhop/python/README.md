cellhop (Python side)

Install locally:

```bash
pip install -e .
```

Usage:

```python
from anndata import AnnData
from cellhop import anndata_to_seurat

# adata = ...
seurat_obj = anndata_to_seurat(adata)
```
