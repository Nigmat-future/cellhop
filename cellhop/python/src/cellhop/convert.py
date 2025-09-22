import os
import pathlib
import tempfile
from typing import Optional

from anndata import AnnData

try:
    import rpy2.robjects as ro
    from rpy2.robjects import vectors as rv
    from rpy2.robjects.packages import importr
    from rpy2.rinterface_lib.embedded import RRuntimeError
except Exception as exc:  # pragma: no cover
    raise ImportError(
        "rpy2 is required to use cellhop.anndata_to_seurat. Install rpy2 and ensure R is available."
    ) from exc


def _ensure_r_packages() -> None:
    try:
        importr("Seurat")
        importr("SeuratDisk")
    except RRuntimeError as err:
        raise RuntimeError(
            "R packages 'Seurat' and 'SeuratDisk' must be installed in your R library.\n"
            "Install in R: install.packages('Seurat'); remotes::install_github('mojaveazure/seurat-disk')\n"
            "Also ensure R is discoverable by rpy2 (R_HOME set)."
        ) from err


def _as_r_path(path: str) -> rv.StrVector:
    # Use forward slashes to avoid Windows backslash escaping issues in R
    posix = pathlib.Path(path).as_posix()
    return rv.StrVector([posix])


def anndata_to_seurat(
    adata: AnnData,
    *,
    tmpdir: Optional[str] = None,
    overwrite: bool = True,
    keep_tmp: bool = False,
) -> ro.vectors.ListVector:
    """
    Convert an AnnData object to an R Seurat object in one call.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    tmpdir : Optional[str]
        Optional temp directory to use for intermediate files.
    overwrite : bool
        Whether to overwrite intermediate files if they exist.
    keep_tmp : bool
        If True, keep temporary files for debugging.

    Returns
    -------
    rpy2.robjects.vectors.ListVector
        An R Seurat object (rpy2 proxy) ready for use with R functions.
    """
    if not isinstance(adata, AnnData):
        raise TypeError("adata must be an AnnData object")

    _ensure_r_packages()
    seuratdisk = importr("SeuratDisk")

    tmp_context = tempfile.TemporaryDirectory(dir=tmpdir)
    tmp_path = tmp_context.name
    h5ad_path = os.path.join(tmp_path, "cellhop_input.h5ad")
    h5seurat_path = os.path.join(tmp_path, "cellhop_input.h5seurat")

    try:
        # Write AnnData to h5ad
        adata.write_h5ad(h5ad_path)

        # Convert to h5seurat using R SeuratDisk
        seuratdisk.Convert(_as_r_path(h5ad_path), dest="h5seurat", overwrite=overwrite)

        # Load Seurat object
        r_obj = seuratdisk.LoadH5Seurat(_as_r_path(h5seurat_path))
    finally:
        if not keep_tmp:
            try:
                tmp_context.cleanup()
            except Exception:
                pass

    return r_obj
