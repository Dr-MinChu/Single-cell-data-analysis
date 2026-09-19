# 设置工作目录
setwd(dir = "D:/数据分析/单细胞数据分析/FTD_data_202512/wd_R_V3")
# 加载相关库
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(harmony)
library(ggrepel)
library(CellChat)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(NMF)
library(ggalluvial)
library(pheatmap)
library(future)
library(cowplot)
library(scRNAtoolVis)
library(jjPlot)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggupset)
library(biomaRt)
library(Seurat)
library(dplyr)
library(stringr)
library(DOSE)  
library(STRINGdb)
library(igraph)
library(ggraph)



# 查看当前工作目录
getwd()


#####1、NK cell 数据归一化、筛选高变基因与PCA降维（可加载数据）####
{
  #====读取数据，计算=========================================================
  # 重新加载初步注释完成的数据（可直接执行这行，读取既往数据结果，免得再执行前述读取操作）
  RNA_harmony_annotation_first <- readRDS("RNA_harmony_annotation.rds")
  
  sub_NKcell_first <- subset(RNA_harmony_annotation_first, subset = celltype == "NK cell")
  
  # harmony整合是基于PCA降维结果进行的。
  sub_NKcell_first <- NormalizeData(sub_NKcell_first) %>% # 数据归一化处理
    FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% #筛选高变基因
    ScaleData() %>% #数据标准化
    RunPCA(npcs = 30, verbose = T)#npcs：计算和存储的PC数（默认为 50）
  
  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(sub_NKcell_first, file = "sub_NKcell_first_PCA.rds")
  
  #====读取数据，绘图=========================================================
  sub_NKcell_first <- readRDS("sub_NKcell_first_PCA.rds")
  
  a=DimPlot(sub_NKcell_first,reduction = "pca",group.by = "orig.ident")
  #PCA图看到还有一点的批次效应（融合得比较好批次就弱）
  a
  
  ##看下高变基因有哪些,可视化
  # 提取前15个高变基因ID
  top15 <- head(VariableFeatures(sub_NKcell_first), 15) 
  plot1 <- VariableFeaturePlot(sub_NKcell_first) 
  plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
  # 合并图片
  feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
  feat_15
  # 保存图片
  ggsave(file = "NKcell_feat_15.pdf",plot = feat_15,he = 10,wi = 15 )
}


#####2、NK cell 细胞周期评分####
{
  #====读取数据，计算=========================================================
  sub_NKcell_first <- readRDS("sub_NKcell_first_PCA.rds")
  
  # 提取g2m特征向量
  g2m_genes = cc.genes$g2m.genes
  g2m_genes = CaseMatch(search = g2m_genes, match = rownames(sub_NKcell_first))
  # 提取s期特征向量
  s_genes = cc.genes$s.genes
  s_genes = CaseMatch(search = s_genes, match = rownames(sub_NKcell_first))
  # 对细胞周期阶段进行评分
  sub_NKcell_first <- CellCycleScoring(object=sub_NKcell_first,  g2m.features=g2m_genes,  s.features=s_genes)
  sub_NKcell_first=CellCycleScoring(object = sub_NKcell_first, 
                                   s.features = s_genes, 
                                   g2m.features = g2m_genes, 
                                   set.ident = TRUE)#set.ident 是否给每个细胞标注一个细胞周期标签
  sub_NKcell_first <- CellCycleScoring(object=sub_NKcell_first,  g2m.features=g2m_genes,  s.features=s_genes)
  
  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(sub_NKcell_first, file = "sub_NKcell_first_PCA_g2m.rds")
  
  #====读取数据，绘图=========================================================
  sub_NKcell_first <- readRDS("sub_NKcell_first_PCA_g2m.rds")
  
  # p4=VlnPlot(sub_NKcell_first, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
  #            ncol = 2, pt.size = 0.1)
  # p4
  sub_NKcell_first@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
    theme_minimal()
  
}


#####3、NK cell RunHarmony去批次####
{
  #====读取数据，计算=========================================================
  sub_NKcell_first <- readRDS("sub_NKcell_first_PCA_g2m.rds")
  
  # 整合需要指定Seurat对象和metadata中需要整合的变量名。
  NKcell_scRNA_harmony <- RunHarmony(sub_NKcell_first, group.by.vars = "orig.ident")
  NKcell_scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
  
  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(NKcell_scRNA_harmony, file = "NKcell_scRNA_harmony.rds")
  
  #====读取数据，绘图=========================================================
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony.rds")
  
  b=DimPlot(NKcell_scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
  #PCA图看到还有一点的批次效应（融合得比较好批次就弱）
  b
  # 合并图片
  #pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1) #a需要先运行前述的第7章的绘图才会有
  #pca_harmony_integrated
  #去批次前后对比
  
  # 后续都是基于Harmony矫正之后的数据，不是基因表达数据和直接的PCA降维数据。
  # 设置reduction = 'harmony'，后续分析是基于Harmony矫正之后的数据。
}


#####4、NK cell 聚类、umap/tsne降维（可加载数据）####
{
  #====读取数据，计算=========================================================
  # 重新加载（可直接执行这行，读取既往数据结果，免得再执行前述读取操作）
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony.rds")
  
  # 提取全部元数据列名
  meta_cols <- colnames(NKcell_scRNA_harmony@meta.data)
  # 筛选出所有 RNA_snn_res. 开头的分辨率列
  res_cols <- meta_cols[grepl("^RNA_snn_res\\.", meta_cols)]
  # 删除这些列（旧精度聚类列）
  NKcell_scRNA_harmony@meta.data[, res_cols] <- NULL
  
  
  ElbowPlot(NKcell_scRNA_harmony, ndims=50, reduction="harmony") #Harmony降维后的主成分（PC）数量选择图，一般选择拐点数量作为后续聚类dims

  # 先基于最优PC（之前elbow图拐点）构建近邻图，只运行一次
  NKcell_scRNA_harmony <- FindNeighbors(
    object = NKcell_scRNA_harmony,
    reduction = "harmony",
    dims = 1:11
  )
  
  # 批量循环跑分辨率，步长0.1，全部存入对象
  res_range <- seq(from = 0.1, to = 2.0, by = 0.1)
  NKcell_scRNA_harmony <- FindClusters(
    object = NKcell_scRNA_harmony,
    resolution = res_range,
    verbose = TRUE
  )
  
  # 保存带多分辨率聚类信息的对象（避免重复计算）
  saveRDS(NKcell_scRNA_harmony, "NKcell_scRNA_harmony_multiRes.rds")
  
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_multiRes.rds")
  # 绘制不同精度的聚类树
  #install.packages("clustree")
  library(clustree)
  tree_plot <- clustree(NKcell_scRNA_harmony, prefix = "RNA_snn_res.")
  print(tree_plot)
  ggsave("clustree_resolution_tree_NKcell.pdf", plot = tree_plot, width = 22, height = 20, dpi = 300)
  
  # 根据前述聚类树，固定最优分辨率为0.6，在不使用seurat_clusters时的精度可由这个去固定，但重要的是下面那个
  Idents(NKcell_scRNA_harmony) <- "RNA_snn_res.0.5"#后续寻找差异基因要用这个
  
  table(NKcell_scRNA_harmony@meta.data$seurat_clusters)#这个元数据里的seurat_cluster默认保存是最后一个精度值的簇，比如此时是最大的1.2精度对应44类，需要用指定精度列的簇数值覆盖，这样不影响后续使用想选的精度。
  
  # 重要：彻底固定最优分辨率，将指定精度0.6分辨率的分群，直接覆盖到NKcell_scRNA_harmony@meta.data$seurat_clusters列
  NKcell_scRNA_harmony$seurat_clusters <- NKcell_scRNA_harmony$RNA_snn_res.0.5 ##注意这里精度一定要改为自己想用的精度值
  
  # 验证：现在table出来就是设定的精度簇
  table(NKcell_scRNA_harmony$seurat_clusters)
  
  
    ##umap/tsne降维
  NKcell_scRNA_harmony <- RunTSNE(NKcell_scRNA_harmony, reduction = "harmony", dims = 1:11)
  NKcell_scRNA_harmony <- RunUMAP(NKcell_scRNA_harmony, reduction = "harmony", dims = 1:11)

  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(NKcell_scRNA_harmony, file = "NKcell_scRNA_harmony_umap_tsne.rds")
  
  #====读取数据，绘图=========================================================
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_umap_tsne.rds")

  # 按样本绘图
  umap_integrated1 <- DimPlot(NKcell_scRNA_harmony, reduction = "umap", group.by = "orig.ident")
  umap_integrated2 <- DimPlot(NKcell_scRNA_harmony, reduction = "umap", label = TRUE)
  tsne_integrated1 <- DimPlot(NKcell_scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
  tsne_integrated2 <- DimPlot(NKcell_scRNA_harmony, reduction = "tsne", label = TRUE)
  # 合并图片
  umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
  # 将图片输出到画板
  umap_tsne_integrated
  # 保存图片
  ggsave("umap_tsne_integrated_NKcell.pdf",umap_tsne_integrated,wi=25,he=15)
  
  # 按分组去绘图
  umap_integrated1_group <- DimPlot(NKcell_scRNA_harmony, reduction = "umap", group.by = "group")
  tsne_integrated1_group <- DimPlot(NKcell_scRNA_harmony, reduction = "tsne", group.by = "group") 
  # 合并图片
  umap_tsne_integrated_group <- CombinePlots(list(tsne_integrated1_group,tsne_integrated2,umap_integrated1_group,umap_integrated2),ncol=2)
  # 将图片输出到画板
  umap_tsne_integrated_group
  # 保存图片
  ggsave("umap_tsne_integrated_group_NKcell.pdf",umap_tsne_integrated_group,wi=25,he=15)
  
  table(NKcell_scRNA_harmony@meta.data$seurat_clusters)
  
}



#####5、NK cell FindAllMarkers（可加载数据）####
{
  #====读取数据，计算=========================================================
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_umap_tsne.rds")
  table(NKcell_scRNA_harmony@meta.data$seurat_clusters)#这个元数据里的seurat_cluster默认保存是最后一个精度值的簇，比如此时是最大的1.2精度对应44类，需要用指定精度列的簇数值覆盖，这样不影响后续使用想选的精度。
  
  #分别对每个cluster进行与剩下除他之外所有的cluster的差异分析
  markers <- FindAllMarkers(object = NKcell_scRNA_harmony, test.use="wilcox" ,
                            only.pos = TRUE,
                            logfc.threshold = 0.25)
  # 将每个簇与所有其他簇进行比较，以识别潜在的标记基因。每个簇中的细胞被视为重复，本质上是通过一些统计测试进行差异表达分析。
  #FindAllMarkers()函数进行亚细胞注释的原理：分别对每个cluster做差异分析，找到在每个cluster特异性高表达的基因作为每个cluster的marker基因，然后将marker基因与cellmarker等原有数据库进行比对，看marker基因到底属于哪种celltype。
  #对于FindAllMarkers()函数，参数设置如下：
  #1.1. 默认为Wilcoxon Rank Sum检验
  #1.2. lgFC默认值为0.25，lgFC低于0.25的自动删除，这个参数可调节，一般按0.25为阈值
  #1.3. min.diff.pct:这个参数代表该cluster中表达该基因的细胞百分比与其他类群中表达该基因的百分比之间的最小差异
  #1.4. min.pct:最少的比例，即这个基因在这两个cluster中，至少要有多少的检测比例。如设置为0.1,代表在这两个cluster中，每个cluster中至少有10%的细胞表达这个基因，如果小于10%，就会把这个基因排除在外
  #1.5. only.pos:这个参数表示，是否选择只保留表达上调的基因，在找Marker基因时，我们一般只找阳性的marker基因(个人理解为只在这个cluster中特异性高表达的基因)，一般不要阴性的
  
  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(markers, file = "markers_01_NKcell.rds")
  markers <- readRDS("markers_01_NKcell.rds")
  # 对计算好的每cluster的marker基因进行筛选
  all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
  #筛选出P<0.05的marker基因
  top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) #将每个cluster lgFC排在前15的marker基因挑选出来
  # top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
  View(top15)
  write.csv(top15,"cluster_top15_NKcell.csv",row.names = T)
  write.csv(all.markers, "cluster_allmarkers_NKcell.csv", row.names = F)
  
  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(top15, file = "NKcell_scRNA_harmony_top15_markers.rds")
  

  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(NKcell_scRNA_harmony, file = "NKcell_scRNA_harmony_markers.rds")
  
  #====读取数据，绘图=========================================================
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_markers.rds")
  top15 <- readRDS("NKcell_scRNA_harmony_top15_markers.rds")
  
  # 可视化maker
  DoHeatmap(NKcell_scRNA_harmony, features = top15$gene, slot="data") + NoLegend()#slot默认使用scaledata里边只有2k个高变gene表达 这里要使用data数据，不然有些gene会找不到表达量
  VlnPlot(NKcell_scRNA_harmony,features = top15$gene[1:15])#选前15个makergene看看,正好前15行就是cluster0的差异基因
  
  #先看下前15行基因的点图，恰属于cluster0，这里可以任意改top15$gene范围或者具体的基因
  p <- DotPlot(NKcell_scRNA_harmony, features = top15$gene[1:15],
               assay='RNA' ,group.by = 'seurat_clusters' ) + coord_flip()+ggtitle("")
  p
  
  # 或者自己选gene看看
  VlnPlot(NKcell_scRNA_harmony,features = c("CD3D","CD3E"))
  
  View(NKcell_scRNA_harmony@meta.data)#seurat_clusters这一列存放了每个细胞对应的cluster
  table(NKcell_scRNA_harmony@meta.data$seurat_clusters)#可以查看每个cluster有多少个细胞

}




#####6、NK cell 进行细胞注释（使用细化六类NK子群）####
{
  #====读取数据，计算=========================================================
  # 重新加载（可直接执行这行，读取既往数据结果，免得再执行前述读取操作）
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_markers.rds")
  
  library(ggplot2) 
  
  #利用整理的一些markers，例如制作的外周血makers文件：blood_cell_markers_NKcell.rds
  # 使用marker基因库进行细胞类型注释
  blood_NKcell_markers <- readRDS("blood_cell_markers_NKcell_fine_subsets.rds")  #加载提前自行生成的基因库
  
  NKcell_harmony_annotation <- NKcell_scRNA_harmony #重命名要注释的seurat对象，以跟之前未注释状态作区分
  
  ## 自动细胞注释函数（deepseek生成）
  # 基于cluster的细胞类型注释（每个cluster一个类型）
  # @param seurat_obj , Seurat对象
  # @param markers_list , marker基因列表
  # @return 包含cluster_celltypes列的Seurat对象
  # 最简单直接的解决方案
  #优化逻辑：不再单纯取所有 marker 平均，改用 (本簇均值 - 其他簇均值) 差值打分，优先「该套基因在本簇特异性上调」的类型，消除广谱基因干扰
    ## 自动细胞注释函数（消除未知细胞类型版本）
  simple_cluster_annotation <- function(seurat_obj, markers_list,
                                        pct_cut = 0.15    # 基因阳性细胞占比阈值
  ) {
    # 1. 依赖包前置检测
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("依赖缺失：请先加载Seurat包")
    }
    if (!requireNamespace("plyr", quietly = TRUE)) {
      stop("依赖缺失：请先加载plyr包")
    }
    
    # 2. 校验聚类列存在
    if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      stop("Seurat对象中无seurat_clusters列，请先运行聚类分析！")
    }
    
    # 3. 统一转为数值型
    seurat_obj@meta.data$seurat_clusters <- as.numeric(
      as.character(seurat_obj@meta.data$seurat_clusters)
    )
    
    # 4. 提取有效cluster
    clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters[!is.na(seurat_obj@meta.data$seurat_clusters)]))
    cluster_celltypes <- list()
    
    # =====================【修正权重：NK1C压低，其余亚群提高竞争力】=====================
    weight_rule <- list(
      "NK1A" = list(
        core = c("CXCR4", "JUN", "JUNB"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      ),
      "NK1B" = list(
        core = c("CD160", "NEAT1", "IFITM1"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      ),
      "NK1C" = list(
        core = c("PRF1", "GZMA", "GZMB", "PTGDS", "ACTB", "CFL1"),
        aux = c(),
        w_core = 1.5,  # 压低毒性基因权重，避免全簇归为NK1C
        w_aux = 1
      ),
      "NKint" = list(
        core = c("CXCR4", "ZFP36", "IER2"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      ),
      "NK2" = list(
        core = c("IL7R", "SELL", "LTB", "FLT3LG"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      ),
      "NK3" = list(
        core = c("KLRC2", "GZMH", "IL32", "ASCL2"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      )
    )
    
    
    # 打印权重配置用于核对
    cat("\n================ 当前NK亚群权重配置一览 ================\n")
    for (cell_name in names(weight_rule)) {
      rule <- weight_rule[[cell_name]]
      cat(sprintf("【%s】\n", cell_name))
      cat(sprintf("  core基因：%s\n", paste0(rule$core, collapse = ",")))
      cat(sprintf("  aux 基因：%s\n", ifelse(length(rule$aux)==0, "无", paste0(rule$aux, collapse = ","))))
      cat(sprintf("  核心权重w_core = %.2f , 辅助权重w_aux = %.2f\n\n", rule$w_core, rule$w_aux))
    }
    cat("==================================================\n\n")
    
    # 5. 强制校验marker列表与权重列表名称完全匹配
    if (!identical(sort(names(markers_list)), sort(names(weight_rule)))) {
      stop("错误：markers_list名称必须为 NK1A/NK1B/NK1C/NKint/NK2_sub/NK3_sub，请核对！")
    }
    
    cat("========= NK亚群加权注释启动 =========\n")
    
    # 6. 逐簇循环打分
    for (cluster in clusters) {
      cell_idx <- which(seurat_obj@meta.data$seurat_clusters == cluster)
      cell_bar <- colnames(seurat_obj)[cell_idx]
      n_cell <- length(cell_bar)
      
      # 空簇跳过
      if (n_cell == 0) {
        cluster_celltypes[[as.character(cluster)]] <- "空Cluster"
        cat(sprintf("Cluster %d：空Cluster，细胞数0\n", cluster))
        next
      }
      
      expr_all <- Seurat::GetAssayData(seurat_obj, layer = "data", assay = "RNA")
      all_gene <- rownames(expr_all)
      expr_cl <- expr_all[, cell_bar, drop = F]
      rm(expr_all)
      gc()
      
      score_list <- c()
      for (ctype in names(markers_list)) {
        rule <- weight_rule[[ctype]]
        all_marker <- markers_list[[ctype]]
        marker_use <- intersect(all_marker, all_gene)
        
        if (length(marker_use) == 0) {
          score_list[ctype] <- 0
          next
        }
        
        core_g <- intersect(marker_use, rule$core)
        aux_g <- intersect(marker_use, rule$aux)
        
        valid_core <- c()
        valid_aux <- c()
        # 筛选阳性率达标核心基因
        if (length(core_g) > 0) {
          for (g in core_g) {
            p <- sum(expr_cl[g,] > 0) / n_cell
            if (p >= pct_cut) valid_core <- c(valid_core, g)
          }
        }
        # 筛选阳性率达标辅助基因
        if (length(aux_g) > 0) {
          for (g in aux_g) {
            p <- sum(expr_cl[g,] > 0) / n_cell
            if (p >= pct_cut) valid_aux <- c(valid_aux, g)
          }
        }
        
        core_score <- 0
        if (length(valid_core) > 0) {
          core_avg <- mean(rowMeans(expr_cl[valid_core, , drop=F], na.rm=T), na.rm=T)
          core_score <- core_avg * rule$w_core
        }
        aux_score <- 0
        if (length(valid_aux) > 0) {
          aux_avg <- mean(rowMeans(expr_cl[valid_aux, , drop=F]), na.rm=T)
          aux_score <- aux_avg * rule$w_aux
        }
        total_score <- core_score + aux_score
        score_list[ctype] <- total_score
      }
      
      # =========【核心改动：取消阈值判断，永远取最高分亚型】=========
      max_score <- max(score_list)
      final_ct <- names(which.max(score_list))
      # =============================================================
      
      cluster_celltypes[[as.character(cluster)]] <- final_ct
      cat(sprintf("Cluster %d | %s | 最高分：%.3f | 细胞数：%d\n",
                  cluster, final_ct, max_score, n_cell))
      
      rm(expr_cl)
      gc()
    }
    
    # 映射赋值celltype
    map_df <- data.frame(
      cluster = as.numeric(names(cluster_celltypes)),
      celltype = unlist(cluster_celltypes),
      stringsAsFactors = F
    )
    
    seurat_obj$celltype <- plyr::mapvalues(
      x = seurat_obj$seurat_clusters,
      from = map_df$cluster,
      to = map_df$celltype,
      warn_missing = F
    )
    # 理论不会触发NA兜底，保留安全防线
    seurat_obj$celltype[is.na(seurat_obj$celltype)] <- "NK1_main"
    
    # 汇总输出
    cat("\n========= NK亚群注释汇总结果 =========\n")
    for (i in 1:nrow(map_df)) {
      cl <- map_df$cluster[i]
      ct <- map_df$celltype[i]
      cnt <- sum(seurat_obj$seurat_clusters == cl, na.rm=T)
      cat(sprintf("Cluster %s：%s，共%d细胞\n", cl, ct, cnt))
    }
    
    return(seurat_obj)
  }
  
  
  
  # 使用示例,运行这个简单版本进行自动注释
  NKcell_harmony_annotation <- simple_cluster_annotation(
    seurat_obj = NKcell_harmony_annotation,
    markers_list = blood_NKcell_markers
  )
  
  table(NKcell_harmony_annotation@meta.data$celltype,NKcell_harmony_annotation@meta.data$seurat_clusters)
  
  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(NKcell_harmony_annotation, file = "NKcell_harmony_annotation_fine_subsets.rds")
  
  #====读取数据，绘图=========================================================
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  ##注释结果绘图TSNE，并保存文件
  p.dim.cell=DimPlot(NKcell_harmony_annotation, reduction = "tsne", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
  p.dim.cell
  ggsave(plot=p.dim.cell,filename="DimPlot_tsne_NK_celltype_fine_subsets.pdf",width=9, height=7) #PDF无法加载中文名称，故细胞类型名可改为英文，或者保存图片格式而非PDF
  
  ##注释结果绘图UMAP，并保存文件
  p.dim.cell=DimPlot(NKcell_harmony_annotation, reduction = "umap", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
  p.dim.cell
  ggsave(plot=p.dim.cell,filename="DimPlot_umap_NK_celltype_fine_subsets.pdf",width=9, height=7) #PDF无法加载中文名称，故细胞类型名可改为英文，或者保存图片格式而非PDF
  
  
  # 1. 统计“细胞类型 × 样本”的绝对数量（核心交叉表）
  # 格式：table(细胞类型列, 样本列)
  cell_count_cross <- table(
    celltype = NKcell_harmony_annotation$celltype,  # 替换为你的细胞类型列名（如sce$custom_celltype）
    sample = NKcell_harmony_annotation$orig.ident   # 替换为你的样本列名（如sce$sample_id）
  )
  
  # 2. 查看交叉表（直观看到所有样本-细胞类型的数量分布）
  print("=== 所有样本的细胞数量交叉表 ===")
  print(cell_count_cross)

  #生成表格比较细胞比例
  # 构建交叉计数矩阵
  ct <- table(NKcell_harmony_annotation$orig.ident, NKcell_harmony_annotation$celltype)
  
  # 转成数据框，行名是样本
  df_wide <- as.data.frame.matrix(ct)
  # 新增样本ID列放到最前面
  df_wide$sample_id <- rownames(df_wide)
  df_wide <- df_wide[, c("sample_id", colnames(df_wide)[1:(ncol(df_wide)-1)])]
  
  # 输出
  write.csv(df_wide, file.path("D:/数据分析/单细胞数据分析/FTD_data_202512/wd_R_V3/NK", "Sub_NKcell_count_wide_fine_subsets.csv"), row.names = F, fileEncoding = "UTF-8")
  
  # 打印预览核对
  print(head(df_wide, 3))
  
  
  
  # 1. 提取患者组细胞ID
  patient_cells <- rownames(NKcell_harmony_annotation@meta.data)[NKcell_harmony_annotation$group == "patient"]
  
  # 2. 单独绘制患者组TSNE（保留你原有参数：label=T、repel=TRUE、pt.size=1）
  p.dim.cell.patient.tsne <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = patient_cells  # 关键参数：仅绘制患者组细胞
  ) 
  p.dim.cell.patient.tsne
  
  # 3. 保存患者组TSNE图（文件名区分患者组）
  ggsave(plot=p.dim.cell.patient.tsne, filename="DimPlot_tsne_sub_celltype_fine_subsets_patient.pdf", width=9, height=7)
  
  # 4. 单独绘制患者组UMAP（同逻辑）
  p.dim.cell.patient.umap <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = patient_cells
  ) 
  p.dim.cell.patient.umap
  ggsave(plot=p.dim.cell.patient.umap, filename="DimPlot_umap_sub_celltype_fine_subsets_patient.pdf", width=9, height=7)
  
  
  # 1. 提取对照组细胞ID
  control_cells <- rownames(NKcell_harmony_annotation@meta.data)[NKcell_harmony_annotation$group == "control"]
  
  # 2. 单独绘制对照组TSNE
  p.dim.cell.control.tsne <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = control_cells  # 关键参数：仅绘制对照组细胞
  ) 
  p.dim.cell.control.tsne
  ggsave(plot=p.dim.cell.control.tsne, filename="DimPlot_tsne_sub_celltype_fine_subsets_control.pdf", width=9, height=7)
  
  # 3. 单独绘制对照组UMAP
  p.dim.cell.control.umap <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = control_cells
  ) 
  p.dim.cell.control.umap
  ggsave(plot=p.dim.cell.control.umap, filename="DimPlot_umap_sub_celltype_fine_subsets_control.pdf", width=9, height=7)
  
  
  # 1. 提取carrier组细胞ID
  carrier_cells <- rownames(NKcell_harmony_annotation@meta.data)[NKcell_harmony_annotation$group == "carrier"]
  
  # 2. 单独绘制carrier组TSNE
  p.dim.cell.carrier.tsne <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = carrier_cells  # 仅绘制carrier组细胞
  ) 
  p.dim.cell.carrier.tsne
  ggsave(plot=p.dim.cell.carrier.tsne, filename="DimPlot_tsne_sub_celltype_fine_subsets_carrier.pdf", width=9, height=7)
  
  # 3. 单独绘制carrier组UMAP
  p.dim.cell.carrier.umap <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = carrier_cells
  ) 
  p.dim.cell.carrier.umap
  ggsave(plot=p.dim.cell.carrier.umap, filename="DimPlot_umap_sub_celltype_fine_subsets_carrier.pdf", width=9, height=7)
  
  
}



#####6.1 气泡图：NK1B NK1C NK3亚群特征Marker基因 DotPlot（使用副本，原始对象不改动）####

{
  # 读取NK注释对象
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  # 建立绘图中间副本，原始对象不改动
  NK_plot <- NKcell_harmony_annotation
  
  # marker基因
  marker_list <- list(
    NK1B = c("CD160", "NEAT1", "IFITM1"),
    NK1C = c("PRF1", "GZMA", "GZMB", "PTGDS", "ACTB", "CFL1"),
    NK3  = c("KLRC2", "GZMH", "IL32", "ASCL2")
  )
  marker_genes <- unlist(marker_list, use.names = FALSE)
  
  # 设置ident为细胞类型列celltype
  Idents(NK_plot) <- "celltype"
  
  # 副本取子集，只保留目标3个亚群
  NK_sub <- subset(NK_plot, idents = c("NK1B","NK1C","NK3"))
  # X轴显示顺序：NK1B，NK1C，NK3
  Idents(NK_sub) <- factor(Idents(NK_sub), levels = c("NK1B","NK1C","NK3"))
  
  p_dot <- DotPlot(
    NK_sub,
    features = marker_genes,
    dot.scale = 8,
    scale.by = "size"
  ) +
    coord_flip() +  # 坐标轴翻转：X=细胞亚群，Y=marker基因
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) +
    labs(x = "NK cell subsets",
         y = "Marker genes",
         title = "Marker gene expression of NK1B, NK1C, NK3 subsets")
  
  ggsave("dotplot_NK_subsets_marker_flip.pdf", plot = p_dot, width = 7, height = 10, dpi = 300)
  cat("✅坐标轴翻转气泡图输出完成：dotplot_NK_subsets_marker_flip.pdf\n")
  
  # 输出统计表格
  dot_data <- DotPlot(NK_sub, features = marker_genes)$data
  write.csv(dot_data, "dotplot_NK_subsets_stat.csv", row.names = FALSE)
  cat("✅气泡图统计数据表输出完成：dotplot_NK_subsets_stat.csv\n")
  
}



#####6.2 桑基冲积图：左侧全部NK总群，右侧NK1B/NK1C/NK3 + Other_NK(其余NK亚群合并展示)####

{
  #读取对象
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  #绘图副本，原始对象不修改
  NK_plot <- NKcell_harmony_annotation
  
  meta_df <- NK_plot@meta.data
  
  #筛选全部NK开头的细胞，完整NK总群
  nk_mask <- grepl("^NK", meta_df$celltype, ignore.case = FALSE)
  meta_NK_all <- meta_df[nk_mask, , drop = FALSE]
  
  #分组：保留NK1B/NK1C/NK3，其余NK亚群归为Other_NK
  meta_NK_all$group_plot <- ifelse(
    meta_NK_all$celltype %in% c("NK1B","NK1C","NK3"),
    meta_NK_all$celltype,
    "Other_NK"
  )
  
  #统计各组细胞数
  sankey_count <- as.data.frame(table(target = meta_NK_all$group_plot))
  cat("====NK分组细胞统计====\n")
  print(sankey_count)
  cat("✅NK总细胞数：", sum(sankey_count$Freq),"\n")
  
  #桑基输入：虚拟源头=全部NK总群
  sankey_input <- data.frame(
    source = "NK_total",
    target = sankey_count$target,
    count  = sankey_count$Freq,
    stringsAsFactors = FALSE
  )
  
  write.csv(sankey_input, "sankey_NKtotal_3sub_plusOtherNK.csv", row.names = FALSE)
  cat("✅桑基计数表输出：sankey_NKtotal_3sub_plusOtherNK.csv\n")
  
  p_sankey <- ggplot(sankey_input,
                     aes(y = count,
                         axis1 = source,
                         axis2 = target)) +
    geom_alluvium(aes(fill = target), width = 0.3) +
    geom_stratum(width = 0.3, fill = "gray85", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(
      limits = c("Before sub‑clustering\n(All NK population)", "NK fine subsets"),
      expand = c(0.05,0.05)
    ) +
    scale_fill_manual(values = c("NK1B"="#E64B35","NK1C"="#4DBBD5","NK3"="#3C5488","Other_NK"="#999999")) +
    labs(y = "Cell number",
         title = "Cell distribution: total NK population to fine NK subsets") +
    theme_bw(base_size =12)+
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  ggsave("sankey_NKtotal_with_OtherNK.pdf", plot = p_sankey, width =10, height =6.5, dpi=300)
  cat("✅桑基图输出完成：sankey_NKtotal_with_OtherNK.pdf\n")
  
}



#####6.3 FeaturePlot：GZMB、CD160、KLRC2（对象本身全部为NK细胞，无需再筛选）####

{
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  # 绘图副本，原始对象不受修改
  NK_plot <- NKcell_harmony_annotation
  
  # 对象本身全部是NK细胞，不再做细胞筛选
  feat_genes <- c("GZMB","CD160","KLRC2")
  
  # 绘制多张UMAP FeaturePlot，一张图多面板
  p_feature <- FeaturePlot(
    NK_plot,
    features = feat_genes,
    ncol = 3,          # 一行放3张子图
    pt.size = 0.2,     # 点大小，细胞多可以调小，发黑就改成0.15
    order = TRUE       # 高表达点画在上层，不会被掩盖
  ) &
    theme_bw(base_size =10) &
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  
  ggsave("featurePlot_NK_markers_GZMB_CD160_KLRC2.pdf",
         plot = p_feature, width = 14, height = 5, dpi = 300)
  cat("✅FeaturePlot组合图输出完成：featurePlot_NK_markers_GZMB_CD160_KLRC2.pdf\n")
  
  # =========可选：分别输出每张基因单独的PDF文件=========
  for(g in feat_genes){
    p_single <- FeaturePlot(NK_plot, features = g, pt.size=0.2, order=TRUE) +
      theme_bw(base_size =11)+
      theme(plot.title = element_text(hjust =0.5), panel.grid=element_blank())
    ggsave(paste0("FeaturePlot_",g,"_NK.pdf"), plot=p_single, width=6, height=5, dpi=300)
  }
  cat("✅各基因独立FeaturePlot输出完毕\n")
  
}



#####6.4 NK1B/NK1C/NK3 箱线‑散点叠加 + Wilcoxon秩和检验 + 导出完整P值表####

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
{
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  # 绘图副本，原始对象不修改
  NK_plot <- NKcell_harmony_annotation
  
  gene_list <- c("FOS","JUN","GZMB","PTGDS","CX3CR1")
  Idents(NK_plot) <- "celltype"
  
  subset_list <- c("NK1B","NK1C","NK3")
  # 两两比较组合
  comp_list <- list(
    c("control","patient"),
    c("control","carrier"),
    c("patient","carrier")
  )
  comp_name_vec <- c("control_vs_patient","control_vs_carrier","patient_vs_carrier")
  
  # 用于存储全部wilcoxon结果
  all_pvalue_res <- data.frame()
  
  for(sub in subset_list){
    
    obj_sub <- subset(NK_plot, idents = sub)
    cat(sprintf("\n========== %s 亚群，细胞数：%d ==========\n", sub, ncol(obj_sub)))
    
    expr_mat <- FetchData(obj_sub, vars = gene_list)
    meta_sub <- obj_sub@meta.data[, c("group"), drop = FALSE]
    plot_df <- cbind(meta_sub, expr_mat)
    
    plot_long <- plot_df %>%
      pivot_longer(cols = all_of(gene_list),
                   names_to = "gene",
                   values_to = "expression")
    
    plot_long$group <- factor(plot_long$group, levels = c("control","patient","carrier"))
    
    #=====================批量计算wilcoxon p值=====================
    for(g in gene_list){
      for(i in seq_along(comp_list)){
        g1 <- comp_list[[i]][1]
        g2 <- comp_list[[i]][2]
        cname <- comp_name_vec[i]
        
        d1 <- filter(plot_long, gene == g, group == g1)$expression
        d2 <- filter(plot_long, gene == g, group == g2)$expression
        
        # 两组都要有至少2个细胞才做检验，否则p记为NA
        if(length(d1)>=2 && length(d2)>=2){
          wt <- wilcox.test(d1, d2)
          pval <- wt$p.value
        }else{
          pval <- NA
        }
        
        row_tmp <- data.frame(
          subset = sub,
          gene = g,
          comparison = cname,
          group1 = g1,
          group2 = g2,
          p_value = pval,
          stringsAsFactors = FALSE
        )
        all_pvalue_res <- rbind(all_pvalue_res, row_tmp)
      }
    }
    
    #=====================绘图=====================
    p <- ggplot(plot_long, aes(x = group, y = expression, fill = group)) +
      geom_boxplot(outlier = FALSE, width = 0.5, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 0.3, alpha = 0.4, color = "black") +
      scale_fill_manual(values = c("control"="#66c2a5","patient"="#fc8d62","carrier"="#8da0cb")) +
      geom_signif(comparisons = comp_list,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  textsize = 3,
                  step_increase = 0.12) +
      facet_wrap(~gene, scales = "free_y", ncol = 5) +
      labs(x = "Group", y = "Normalized expression", fill = "Group",
           title = sprintf("Gene expression in %s subset across groups", sub)) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(sprintf("box_jitter_%s_5genes_wilcoxon.pdf",sub),
           plot = p, width = 17, height = 5.2, dpi = 300)
    
    write.csv(plot_long, sprintf("boxjitter_%s_gene_expression_table.csv",sub), row.names = FALSE)
    cat(sprintf("✅%s 图与数据表输出完成\n", sub))
  }
  
  # 输出全部wilcoxon检验p值总表
  write.csv(all_pvalue_res, "wilcoxon_all_pvalue.csv", row.names = FALSE)
  cat("\n✅全部亚群‑基因‑两两比较P值总表输出：wilcoxon_all_pvalue.csv\n")
  print(head(all_pvalue_res,10))
  
  
}













