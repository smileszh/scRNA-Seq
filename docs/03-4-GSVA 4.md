



# GSVA 分析
> 单细胞数据的GSVA和芯片、bulk转录组的GSVA没有本质区别，就使用AverageExpression获取平均表达量得到新的表达矩阵再计算即可。

## 加载数据和R包

获得每种细胞的平均表达量.

这里的示例数据seu.obj.Rdata是GSE218208降维聚类分群的结果，因为文件太大，没有直接放进文件夹里，如果load报错就自己运行一下隔壁GSE218208的代码得到这个文件再跑。

![](https://upload-images.jianshu.io/upload_images/9475888-2e226458a1b95f77.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


``` r
rm(list = ls())
library(Seurat)
library(GSVA)
library(clusterProfiler)
load("../2.GSE218208/seu.obj.Rdata")
table(Idents(seu.obj))
## 
##  Naive CD4 T   CD14+ Mono            B        CD8 T           NK FCGR3A+ Mono 
##         1675         1206          598          406          337          125 
##     Platelet           DC 
##           48           88
exp  =  AverageExpression(seu.obj)[[1]]
#exp =  AggregateExpression(seu.obj)[[1]]
exp  =  as.matrix(exp)
exp  =  exp[rowSums(exp)>0,] 
exp[1:4,1:4]
##          Naive CD4 T  CD14+ Mono          B      CD8 T
## TSPAN6    0.01890007 0.000000000 0.00000000 0.00446691
## DPM1      0.50764534 0.398461857 0.52602493 0.49951298
## SCYL3     0.10701976 0.049771894 0.10397003 0.12101561
## C1orf112  0.02653607 0.005093801 0.05426134 0.02747031
```

Seurat v5 提示建议用AggregateExpression做伪bulk转录组分析，那个是用来求和的，目前查到的文献和教程都是使用平均值，这里就木有改动.

## 做GSVA

gmt文件下载自GSEA-msigdb官网


``` r
h_df = read.gmt("h.all.v2023.2.Hs.symbols.gmt")[,c(2,1)]
h_list = unstack(h_df)
# ES  =  gsva(exp, h_list) #⭐R 语言版本4.3运行这一句，代替下面的两句
gsvapar <- gsvaParam(exp, h_list, maxDiff=TRUE) 
ES <- gsva(gsvapar)
ES[1:4,1:4]
```


## 热图可视化


``` r
library(pheatmap)
pheatmap(ES, scale = "row",angle_col = "45",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
```

<img src="03-4-GSVA_files/figure-html/unnamed-chunk-3-1.png" width="672" style="display: block; margin: auto;" />





