# cellhop

一行代码在 Python AnnData 和 R Seurat 之间无缝转换，无需处理中间文件。

**English version**: [README.md](README.md)

## 解决的问题

在单细胞分析流程中，研究者经常需要在 Python 和 R 之间切换数据格式：
- **Python 生态**: 使用 AnnData 进行高效的数据操作和分析
- **R 生态**: 使用 Seurat 进行专业的单细胞分析和可视化

### 传统方法的痛点

目前的转换流程依赖 SeuratDisk，需要手动处理中间文件：

```python
# Python → R 传统方法
import anndata as ad
from rpy2.robjects import r

# 保存为中间文件
adata.write_h5ad('temp.h5ad')

# 使用 SeuratDisk 转换
r('SeuratDisk::Convert("temp.h5ad", dest="h5seurat")')

# 在 R 中加载
r('seurat_obj <- SeuratDisk::LoadH5Seurat("temp.h5seurat")')

# 手动清理文件
import os
os.remove('temp.h5ad')
os.remove('temp.h5seurat')
```

```r
# R → Python 传统方法
library(SeuratDisk)

# 保存为中间文件
SeuratDisk::SaveH5Seurat(seurat_obj, 'temp.h5seurat')

# 转换格式
SeuratDisk::Convert('temp.h5seurat', dest='h5ad')

# 在 Python 中加载
# 需要手动处理文件路径和清理
```

**问题：**
- 🔴 需要手动管理临时文件
- 🔴 中间文件占用磁盘空间
- 🔴 增加了出错的可能性
- 🔴 代码冗长，工作流程不流畅

### cellhop 的解决方案

**一行代码转换，零中间文件暴露！**

#### Python 中：AnnData → Seurat

```python
from cellhop import anndata_to_seurat

# 你的数据
seurat_obj = anndata_to_seurat(adata)
```

#### R 中：Seurat → AnnData

```r
library(cellhop)

# 你的数据
adata <- seurat_to_anndata(seurat_obj)
```

**优势：**
- ✅ **一行代码**：简洁的API设计
- ✅ **零文件管理**：临时文件自动创建和清理
- ✅ **无缝集成**：直接在内存中转换
- ✅ **环境友好**：不会在你的工作目录中留下垃圾文件

## 作者

**Nigmat Rahim** (<nigmatrahim@stu.pku.edu.cn>)

## 安装

### Python 本地安装

```bash
cd python
pip install -e .
```

### R 本地安装

```r
# 从项目根目录
remotes::install_local("R")
# 或使用 devtools
devtools::install("R")
```

## 系统要求

- **R**: Seurat, SeuratDisk (及其依赖如 hdf5r)
- **Python**: anndata, rpy2
- **桥梁**: rpy2 (Python→R), reticulate (R→Python)

在 Windows/macOS 上，确保 R 可被发现 (rpy2 需要 R_HOME)，并配置 reticulate 指向包含 anndata 的 Python 环境 (如需要可设置 `RETICULATE_PYTHON`)。

## 使用方法

### Python (AnnData → Seurat)

```python
from anndata import AnnData
from cellhop import anndata_to_seurat

# adata = ... 你的 AnnData 对象
seurat_obj = anndata_to_seurat(adata)
```

### R (Seurat → AnnData)

```r
library(cellhop)
# seurat_obj <- ... 你的 Seurat 对象
adata <- seurat_to_anndata(seurat_obj)
```

临时文件在后台自动创建和清理，对最终用户完全透明。

## 测试

运行冒烟测试验证你的环境设置：

**Python 端:**
```bash
cd python
python test_smoke.py
```

**R 端:**
```bash
cd R
Rscript test_smoke.R
```

或在 R 控制台中：
```r
source("R/test_smoke.R")
```

## 示例

查看使用示例：
- `demo.py` - Python 工作流示例
- `demo.R` - R 工作流示例
- `test_roundtrip.ipynb` - 显示往返转换的 Jupyter notebook

## 注意事项

- 数据模型差异可能导致某些字段近似映射。请根据你的工作流程检查结果。
- 大对象会在转换期间短暂占用磁盘空间；确保临时目录有足够空间。
- 测试会在完整环境设置的情况下检查包结构和基本功能。

## 工作原理

cellhop 在后台使用成熟的 SeuratDisk 进行实际转换，但提供了更友好的接口：

1. **Python → R**: 使用 rpy2 调用 SeuratDisk::Convert
2. **R → Python**: 使用 reticulate 调用 anndata 读取
3. **文件管理**: 自动创建和管理临时文件，用户无需关心

这样你就可以专注于数据分析，而不是文件管理！🚀
