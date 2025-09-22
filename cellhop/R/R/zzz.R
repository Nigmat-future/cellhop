.onAttach <- function(libname, pkgname) {
  if (requireNamespace("reticulate", quietly = TRUE)) {
    py <- tryCatch(reticulate::py_config(), error = function(e) NULL)
    if (is.null(py)) {
      packageStartupMessage(
        "reticulate is installed. Set RETICULATE_PYTHON to a Python with 'anndata' installed if needed."
      )
    }
  }
}
