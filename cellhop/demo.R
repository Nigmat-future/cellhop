# Demo of cellhop workflow (requires full environment setup)

demo_workflow <- function() {
  message("cellhop Demo")
  message("============")
  message()

  message("This demo shows the intended workflow:")
  message()

  message("1. In Python:")
  message("
    from anndata import AnnData
    from cellhop import anndata_to_seurat
    import numpy as np

    # Create test data
    adata = AnnData(X=np.random.rand(100, 500))
    adata.obs_names = [f'cell_{i}' for i in range(100)]
    adata.var_names = [f'gene_{i}' for i in range(500)]

    # Convert to Seurat (one-liner!)
    seurat_obj = anndata_to_seurat(adata)
    print(f'Converted to Seurat: {type(seurat_obj)}')
  ")

  message("2. In R:")
  message("
    library(Seurat)
    library(cellhop)

    # Create test Seurat object
    counts <- matrix(rpois(100*500, lambda=5), nrow=500, ncol=100)
    seurat_obj <- CreateSeuratObject(counts)

    # Convert to AnnData (one-liner!)
    adata <- seurat_to_anndata(seurat_obj)
    print(class(adata))
    print(adata$shape)
  ")

  message("3. Round-trip conversion:")
  message("   AnnData → Seurat → AnnData")
  message("   (verify data integrity after conversion)")
  message()

  message("Setup requirements:")
  message("- Python: pip install anndata rpy2")
  message("- R: install.packages(c('Seurat', 'reticulate'))")
  message("- R: remotes::install_github('mojaveazure/seurat-disk')")
  message("- Ensure R_HOME and RETICULATE_PYTHON are set correctly")
}

# Run demo
demo_workflow()
