import pandas as pd
import numpy as np
import anndata
from scipy.io import mmread
import os
def seurat_to_adata(counts,      # R中导出的counts.mtx文件
                     meta,        # R中导出的metadata.csv文件
                     gene_name):  # R中导出的gene_names.csv文件

    # 加载表达矩阵
    X = mmread(counts)  # 读取稀疏矩阵格式的counts文件
    X = X.transpose().tocsc()  # 转置并转换为CSC格式（适合稀疏矩阵）

    # 创建 AnnData 对象
    adata = anndata.AnnData(X=X)

    # 读取 metadata
    cell_meta = pd.read_csv(meta)
    adata.obs = cell_meta  # 将 metadata 赋给 adata.obs
    adata.obs.index = adata.obs['barcode']  # 使用 'barcode' 作为索引

    # 读取基因名称
    with open(gene_name, 'r') as f:
        gene_names = f.read().splitlines()  # 基因名称列表
    adata.var.index = gene_names  # 将基因名赋给 adata.var.index

    return adata
adata = seurat_to_adata(
    counts='/home/guile2/CD8Tex/重整泛癌数据/Tcell/h5ad/1267/counts.mtx',
    meta='/home/guile2/CD8Tex/重整泛癌数据/Tcell/h5ad/1267/metadata.csv',
    gene_name='/home/guile2/CD8Tex/重整泛癌数据/Tcell/h5ad/1267/gene_names.csv'
)
adata.write_h5ad('/home/guile2/CD8Tex/重整泛癌数据/Tcell/h5ad/1267/adata.h5ad')
