# 细胞亚群响应注释

通过对marker基因进行GO和KEGG注释，可以研究这些细胞在不同实验条件下的行为和响应，发现与某些疾病（如癌症、自身免疫疾病）相关的特定细胞类型及其特征基因表达。



## 加载包

``` r
library(COSG)
library(harmony)
library(ggsci)
library(future)
library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(patchwork)
library(stringr)
library(SingleR)
library(tidyverse)
library(tidydr)
```

## 提取各个亚群的marker基因

``` r
sce.all <- readRDS('sce.all.rds')

Idents(sce.all) <- 'celltype'
marker_cosg <- cosg(
    sce.all,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=100)
```



