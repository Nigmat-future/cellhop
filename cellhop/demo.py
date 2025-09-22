#!/usr/bin/env python3
"""
Demo of cellhop workflow (requires full environment setup)
"""

# This demo shows the intended workflow
# Uncomment and run after setting up dependencies

def demo_workflow():
    print("cellhop Demo")
    print("============")
    print()

    print("This demo shows the intended workflow:")
    print()

    print("1. In Python:")
    print("""
    from anndata import AnnData
    from cellhop import anndata_to_seurat
    import numpy as np

    # Create test data
    adata = AnnData(X=np.random.rand(100, 500))
    adata.obs_names = [f"cell_{i}" for i in range(100)]
    adata.var_names = [f"gene_{i}" for i in range(500)]

    # Convert to Seurat (one-liner!)
    seurat_obj = anndata_to_seurat(adata)
    print(f"Converted to Seurat: {type(seurat_obj)}")
    """)

    print("2. In R:")
    print("""
    library(Seurat)
    library(cellhop)

    # Create test Seurat object
    counts <- matrix(rpois(100*500, lambda=5), nrow=500, ncol=100)
    seurat_obj <- CreateSeuratObject(counts)

    # Convert to AnnData (one-liner!)
    adata <- seurat_to_anndata(seurat_obj)
    print(class(adata))
    print(adata$shape)
    """)

    print("3. Round-trip conversion:")
    print("   AnnData → Seurat → AnnData")
    print("   (verify data integrity after conversion)")
    print()

    print("Setup requirements:")
    print("- Python: pip install anndata rpy2")
    print("- R: install.packages(c('Seurat', 'reticulate'))")
    print("- R: remotes::install_github('mojaveazure/seurat-disk')")
    print("- Ensure R_HOME and RETICULATE_PYTHON are set correctly")

if __name__ == "__main__":
    demo_workflow()
