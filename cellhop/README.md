## cellhop

One-liners to hop between Python AnnData and R Seurat without exposing intermediate files.

- In Python: convert AnnData → Seurat using rpy2 + SeuratDisk under the hood
- In R: convert Seurat → AnnData using reticulate + anndata under the hood

**中文版介绍**: [README_zh.md](README_zh.md)

**Author:** Nigmat Rahim (<nigmatrahim@stu.pku.edu.cn>)

### Install

Python (local):

```bash
cd python
pip install -e .
```

R (local):

```r
# from repo root
remotes::install_local("R")
# or using devtools
devtools::install("R")
```

### Requirements

- R: Seurat, SeuratDisk (and their dependencies like hdf5r)
- Python: anndata, rpy2
- Bridges: rpy2 (Python→R), reticulate (R→Python)

On Windows/macOS, ensure R is installed and discoverable (R_HOME for rpy2), and reticulate points to a Python with anndata installed (set `RETICULATE_PYTHON` if needed).

### Usage

Python (AnnData → Seurat):

```python
from anndata import AnnData
from cellhop import anndata_to_seurat

# adata = ...
seurat_obj = anndata_to_seurat(adata)
```

R (Seurat → AnnData):

```r
library(cellhop)
# seurat_obj <- ...
adata <- seurat_to_anndata(seurat_obj)
```

Temporary files are created behind the scenes and cleaned up automatically.

### Testing

Run smoke tests to verify your setup:

**Python side:**
```bash
cd python
python test_smoke.py
```

**R side:**
```bash
cd R
Rscript test_smoke.R
```

Or in R console:
```r
source("R/test_smoke.R")
```

### Demo

See example usage in:
- `demo.py` - Python workflow examples
- `demo.R` - R workflow examples
- `test_roundtrip.ipynb` - Jupyter notebook showing round-trip conversion

### Notes

- Data model differences mean some fields may be mapped approximately. Inspect results for your workflow.
- Large objects will briefly hit disk during conversion; ensure enough space in your temp directory.
- The tests check package structure and basic functionality even without full environment setup.
