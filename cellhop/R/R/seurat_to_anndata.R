#' Convert a Seurat object to a Python AnnData object (one-liner)
#'
#' @param seurat_obj A Seurat object
#' @return A reticulate Python AnnData object
#' @examples
#' # adata <- seurat_to_anndata(seurat_obj)
seurat_to_anndata <- function(seurat_obj) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    stop("Package 'SeuratDisk' is required. Install it in R.")
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required. Install it in R.")
  }

  tmpdir <- tempdir()
  h5s <- tempfile(pattern = "cellhop_", tmpdir = tmpdir, fileext = ".h5seurat")
  h5a <- sub("\\.h5seurat$", ".h5ad", h5s)

  on.exit({
    suppressWarnings(try(unlink(c(h5s, h5a)), silent = TRUE))
  }, add = TRUE)

  # Save to h5Seurat and convert to h5ad
  SeuratDisk::SaveH5Seurat(seurat_obj, filename = h5s, overwrite = TRUE)
  SeuratDisk::Convert(h5s, dest = "h5ad", overwrite = TRUE)

  # Load AnnData from Python
  anndata <- reticulate::import("anndata", delay_load = FALSE)
  ad <- anndata$read_h5ad(h5a)
  return(ad)
}
