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
library(tibble)
library(writexl) 

# 清空环境，彻底解决变量冲突
rm(list = ls())


#####1、设置数据读取地址####

# 固定你的数据总路径
root_dir <- "D:/数据分析/单细胞数据分析/FTD_data_202512/MAPT"

# 自动提取文件夹下所有样本子文件夹，无需手动输入
all_sample_folders <- list.dirs(root_dir, full.names = F, recursive = F)
all_sample_folders
# 运行后控制台会打印所有样本文件夹名，方便核对

# 批量校验所有路径是否存在
cat("========== 路径校验结果 ==========\n")
for (sname in all_sample_folders) {
  full_path <- file.path(root_dir, sname)
  exist_flag <- dir.exists(full_path)
  cat("样本：", sname, " | 完整路径：", full_path, " | 是否存在：", exist_flag, "\n")
}
cat("==================================\n\n")

#####2、循环读取patient组、carrier组和control组的数据####
mmRNAList <- list()

for (i in seq_along(all_sample_folders)) {
  data_path <- file.path(root_dir, all_sample_folders[i])
  cat("正在读取样本：", all_sample_folders[i], "完整路径：", data_path, "\n")
  
  # 判断文件夹不存在则跳过，不中断程序
  if (!dir.exists(data_path)) {
    cat("【警告】路径不存在，跳过该样本\n\n")
    next
  }
  
  # 打印文件夹内文件，确认包含matrix/barcodes/features
  cat("文件夹内文件：\n")
  print(list.files(data_path))
  cat("\n")
  
  mmRNA <- CreateSeuratObject(
    counts = Read10X(data.dir = data_path),
    project = all_sample_folders[i],
    min.cells = 5,
    min.features = 200
  )
  
  # 提取分组信息
  mmRNA$group <- gsub("^FTD_|_[^_]+$", "", all_sample_folders[i])
  mmRNA$orig.ident <- all_sample_folders[i]
  
  mmRNAList[[i]] <- mmRNA
}

# 统计读取成功样本数
cat("共成功读取", length(mmRNAList), "个样本\n")

# 保存原始列表
saveRDS(mmRNAList, file = "mmRNAList_origin.rds")


#####3、批量计算线粒体和红细胞比例（可加载数据）####
#====读取数据，计算=========================================================
# 读取原始10X合并数据，无需重新跑读取流程
# 完整绝对路径读取，替换成你存放rds的文件夹
mmRNAList <- readRDS("D:/数据分析/单细胞数据分析/FTD_data_202512/wd_R_V2/mmRNAList_origin.rds")

# 血红蛋白基因全集
HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
total_sample <- length(mmRNAList) # 总样本数

# 循环批量计算比例
for(i in seq_along(mmRNAList)){
  sc <- mmRNAList[[i]]
  # 打印当前进度+样本名
  cat("====================\n")
  cat("正在处理样本：", i, "/", total_sample, "\n")
  cat("样本名称：", sc$orig.ident[1], "\n")
  
  # 线粒体基因比例
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  # 血红蛋白比例，自动过滤不存在的基因
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features = HB_genes)
  
  mmRNAList[[i]] <- sc
  cat("样本",i,"计算完成\n====================\n\n")
}

cat("全部样本mt、HB比例计算完毕，开始保存文件\n")
# 保存新增线粒体、红细胞比例的中间文件
saveRDS(mmRNAList,"mmRNAList_filter_front.rds")
cat("文件保存成功：mmRNAList_filter_front.rds\n")




#####4、批量绘制质控前小提琴图####
{
#====读取数据，绘图=========================================================
mmRNAList <- readRDS("mmRNAList_filter_front.rds")

violin_before <- list()
for(i in 1:length(mmRNAList)){
  violin_before[[i]] <- VlnPlot(mmRNAList[[i]],
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                                pt.size = 0.01, 
                                ncol = 4) 
}
# 合并图片
violin_before_merge <- CombinePlots(plots = violin_before,nrow=length(mmRNAList),legend='none')
# 将图片输出到画板上
violin_before_merge
# 保存图片
ggsave("violin_before_merge.pdf", plot = violin_before_merge, width = 15, height =7)
#查看第XX个样本的小提琴图
violin_before[[4]] 
}





#####4.1 批量绘制质控前小提琴图【修复X轴标签，Idents同步更新】####
{
  # 1. 读取原始list
  mmRNAList_raw <- readRDS("mmRNAList_filter_front.rds")
  
  # 2. 提取原始orig.ident
  raw_id_list <- sapply(mmRNAList_raw, function(x) unique(x$orig.ident))
  
  # 3. 根据原始样本名分配组别标签
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  
  group_tag_vec <- sapply(raw_id_list, get_group_tag)
  
  # 组内独立编号：patient=FTDxxx；control=HCxxx；carrier=Carrierxxx
  mapping_df <- data.frame(
    original_ID = raw_id_list,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = stringr::str_pad(dplyr::row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  # 4. 输出前后对应表到csv
  write.csv(mapping_df, "sample_name_mapping.csv", row.names = FALSE)
  cat("=== 前后对应表已保存：sample_name_mapping.csv ===\n")
  print(mapping_df[,c("original_ID","new_ID")])
  
  new_id_list <- mapping_df$new_ID
  
  # 5. 复制list，同时修改 orig.ident + Idents()
  mmRNAList_plot <- mmRNAList_raw
  for(i in seq_along(mmRNAList_plot)){
    obj <- mmRNAList_plot[[i]]
    obj$orig.ident <- new_id_list[i]
    Idents(obj) <- "orig.ident"
    mmRNAList_plot[[i]] <- obj
  }
  
  # 6. 循环绘图，无总标题（单样本质控前图）
  for(i in seq_along(mmRNAList_plot)){
    p <- VlnPlot(mmRNAList_plot[[i]],
                 features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
                 pt.size = 0.01,
                 ncol = 4)
    ggsave(paste0("violin_before_", new_id_list[i],".pdf"),
           plot = p, width = 14, height = 5)
    cat("✅已完成：", mapping_df$original_ID[i], " -> ", new_id_list[i], "\n")
  }
  cat("====全部样本质控前单样本小提琴图绘制完成====\n")
  
  # ==========新增：质控前‑所有样本合并，分别绘制4个指标总图==========
  seu_merge_plot_before <- merge(mmRNAList_plot[[1]], y = mmRNAList_plot[-1])
  
  feature_all <- c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent")
  title_all <- c(
    "nFeature_RNA across all samples(before filtering)",
    "nCount_RNA across all samples(before filtering)",
    "mt_percent across all samples(before filtering)",
    "HB_percent across all samples(before filtering)"
  )
  file_all <- c(
    "violin_before_allSample_nFeature_RNA.pdf",
    "violin_before_allSample_nCount_RNA.pdf",
    "violin_before_allSample_mt_percent.pdf",
    "violin_before_allSample_HB_percent.pdf"
  )
  
  # 循环批量输出四张合并总图
  for(k in seq_along(feature_all)){
    p_vln <- VlnPlot(seu_merge_plot_before,
                     features = feature_all[k],
                     pt.size = 0.005) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle(title_all[k])
    
    ggsave(file_all[k], plot = p_vln, width = 18, height = 6, dpi = 300)
    cat("✅质控前总图已保存：", file_all[k],"\n")
  }
  
}



#####4.2 批量绘制每个样本 nCount‑RNA VS nFeature‑RNA 散点气泡图（过滤前，颜色=mt_percent）+ 标注Pearson R####

{
  mmRNAList_raw <- readRDS("mmRNAList_filter_front.rds")
  
  #提取原始样本ID
  raw_id_list <- sapply(mmRNAList_raw, function(x) unique(x$orig.ident))
  
  #命名规则和前面保持一致 patient→FTD;control→HC;carrier→Carrier
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  group_tag_vec <- sapply(raw_id_list, get_group_tag)
  
  mapping_df <- data.frame(
    original_ID = raw_id_list,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = str_pad(row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  write.csv(mapping_df,"sample_name_mapping_scatter_front.csv",row.names = FALSE)
  new_id_list <- mapping_df$new_ID
  
  #复制绘图副本
  mmRNAList_plot <- mmRNAList_raw
  for(i in seq_along(mmRNAList_plot)){
    obj <- mmRNAList_plot[[i]]
    obj$orig.ident <- new_id_list[i]
    Idents(obj) <- "orig.ident"
    mmRNAList_plot[[i]] <- obj
  }
  
  #循环绘图：单样本图，增加R值文本标注
  for(i in seq_along(mmRNAList_plot)){
    obj <- mmRNAList_plot[[i]]
    plot_df <- obj@meta.data
    
    #计算Pearson相关系数R
    cor_res <- cor.test(plot_df$nCount_RNA, plot_df$nFeature_RNA, method = "pearson")
    r_val <- round(cor_res$estimate, 3)
    label_text <- paste0("R = ", r_val)
    
    p <- ggplot(plot_df, aes(x = nCount_RNA, y = nFeature_RNA)) +
      geom_point(aes(color = mt_percent), size = 0.4, alpha = 0.6) +
      scale_color_viridis_c(name = "mt%", option = "viridis") +
      # R值文本放在左上角
      annotate("text", x = Inf, y = Inf, label = label_text, 
               hjust = 1.1, vjust = 1.1, size = 4.5) +
      labs(x = "nCount_RNA", y = "nFeature_RNA",
           title = paste0(new_id_list[i])) +
      theme_bw(base_size =11)+
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    #文件名增加 _front
    ggsave(paste0("scatter_nCount_nFeature_mt_",new_id_list[i],"_front.pdf"),
           plot = p, width =7, height =6, dpi=300)
    cat("✅已输出：",mapping_df$original_ID[i]," -> ",new_id_list[i]," | ",label_text,"\n")
  }
  cat("====全部单样本nCount‑nFeature线粒体着色散点图绘制完成====\n")
  
  #====================新增：所有样本汇总总图，同样计算并标注R值====================
  seu_merge_scatter <- merge(mmRNAList_plot[[1]], y = mmRNAList_plot[-1])
  df_all_meta <- seu_merge_scatter@meta.data
  
  cor_all <- cor.test(df_all_meta$nCount_RNA, df_all_meta$nFeature_RNA, method = "pearson")
  r_all_val <- round(cor_all$estimate,3)
  label_all_text <- paste0("R = ", r_all_val)
  
  p_scatter_all <- ggplot(df_all_meta, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(aes(color = mt_percent), size = 0.12, alpha = 0.35) +
    scale_color_viridis_c(name = "mt%", option = "viridis") +
    annotate("text", x = Inf, y = Inf, label = label_all_text, 
             hjust = 1.1, vjust = 1.1, size =4.5) +
    labs(x = "nCount_RNA",
         y = "nFeature_RNA",
         title = "All samples: nCount_RNA vs nFeature_RNA (before filtering)") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  #汇总图文件名 _front
  ggsave("scatter_allSample_nCount_nFeature_mt_front.pdf",
         plot = p_scatter_all, width = 10, height = 8, dpi = 300)
  cat("✅全部样本汇总散点图已输出：scatter_allSample_nCount_nFeature_mt_front.pdf | ",label_all_text,"\n")
  
}








#####5、批量过滤细胞、MT、HB基因####
{

#====读取数据，计算=========================================================
mmRNAList <- readRDS("mmRNAList_filter_front.rds")

#不做质控，数据原已做过。
mmRNAList <- lapply(X = mmRNAList, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & 
                mt_percent < 10 & 
                HB_percent < 3 & 
                nCount_RNA < quantile(nCount_RNA,0.97) & 
                nCount_RNA > 1000)

  # nFeature_RNA：每个细胞检测表达的基因数目大于300，小于5000；
  # nCount_RNA:每个细胞测序的UMI count含量大于1000，且剔除最大的前3%的细胞；
  # mt_percent:每个细胞的线粒体基因表达量占总体基因的比例小于10%；
  # HB_percent：每个细胞红细胞基因表达量占总体基因的比例小于3%。
  })

View(mmRNAList[[1]]@meta.data)
# ps:没有固定的阈值标准，要根据自己的数据调整参数不断尝试，才能找到最佳结果。

saveRDS(mmRNAList,"mmRNAList_filter_after.rds")

}



#####5.1 批量绘制质控后小提琴图【读取mmRNAList_filter_after.rds，复用名称映射csv】####

{
  # 1、读取过滤完成后的seurat list
  mmRNAList_after_raw <- readRDS("mmRNAList_filter_after.rds")
  
  # 2、提取当前对象里面原始orig.ident，重新执行命名规则（不读取外部csv）
  raw_after_id_list <- sapply(mmRNAList_after_raw, function(x) unique(x$orig.ident))
  
  # 定义组别匹配函数，和质控前完全一致
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  
  group_tag_vec <- sapply(raw_after_id_list, get_group_tag)
  
  # 组内编号
  mapping_df_after <- data.frame(
    original_ID = raw_after_id_list,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = stringr::str_pad(dplyr::row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  # 输出质控后自己的对应关系表，文件名区分开
  write.csv(mapping_df_after, "sample_name_mapping_after.csv", row.names = FALSE)
  cat("=== 质控后名称对应表已保存：sample_name_mapping_after.csv ===\n")
  print(mapping_df_after[,c("original_ID","new_ID")])
  
  new_id_list_after <- mapping_df_after$new_ID
  
  # 3、复制列表，仅修改绘图用副本，原始mmRNAList_after_raw保持原样
  mmRNAList_after_plot <- mmRNAList_after_raw
  
  for(i in seq_along(mmRNAList_after_plot)){
    obj <- mmRNAList_after_plot[[i]]
    newid <- new_id_list_after[i]
    
    obj$orig.ident <- newid
    Idents(obj) <- "orig.ident"  # 同步更新Idents，X轴显示FTD001/HC001/Carrier001
    
    mmRNAList_after_plot[[i]] <- obj
  }
  
  # 4、循环绘图输出【每个样本单独】质控后小提琴图
  for(i in seq_along(mmRNAList_after_plot)){
    p <- VlnPlot(mmRNAList_after_plot[[i]],
                 features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
                 pt.size = 0.01,
                 ncol = 4) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) # x轴标签旋转防止文字挤压
    
    ggsave(paste0("violin_after_", new_id_list_after[i],".pdf"),
           plot = p, width = 14, height = 5)
    
    cat("✅质控后绘图完成：", mapping_df_after$original_ID[i], " -> ", new_id_list_after[i],"\n")
  }
  cat("====全部【质控后单样本】小提琴图绘制完毕====\n")
  
  # ==========所有样本合并，分别绘制4个质控指标总图==========
  # 把绘图副本list合并成一个临时对象，用于画总览图
  seu_merge_plot <- merge(mmRNAList_after_plot[[1]], y = mmRNAList_after_plot[-1])
  
  feature_after_all <- c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent")
  title_after_all <- c(
    "nFeature_RNA across all samples(after filtering)",
    "nCount_RNA across all samples(after filtering)",
    "mt_percent across all samples(after filtering)",
    "HB_percent across all samples(after filtering)"
  )
  file_after_all <- c(
    "violin_after_allSample_nFeature_RNA.pdf",
    "violin_after_allSample_nCount_RNA.pdf",
    "violin_after_allSample_mt_percent.pdf",
    "violin_after_allSample_HB_percent.pdf"
  )
  
  # 循环批量输出四张质控后总图
  for(k in seq_along(feature_after_all)){
    p_vln_after <- VlnPlot(seu_merge_plot,
                           features = feature_after_all[k],
                           pt.size = 0.005) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle(title_after_all[k])
    
    ggsave(file_after_all[k], plot = p_vln_after, width = 18, height = 6, dpi = 300)
    cat("✅质控后总图已保存：", file_after_all[k],"\n")
  }
  
}



#####5.2 批量绘制每个样本 nCount‑RNA VS nFeature‑RNA 散点气泡图（过滤后，颜色=mt_percent）+ 标注Pearson R####
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
{
  mmRNAList_after_raw <- readRDS("mmRNAList_filter_after.rds")
  
  #提取过滤后样本原始orig.ident
  raw_after_id_list <- sapply(mmRNAList_after_raw, function(x) unique(x$orig.ident))
  
  #命名规则和过滤前完全保持一致 patient→FTD;control→HC;carrier→Carrier
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  group_tag_vec <- sapply(raw_after_id_list, get_group_tag)
  
  mapping_df_after <- data.frame(
    original_ID = raw_after_id_list,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = str_pad(row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  #质控后映射表，后缀_after
  write.csv(mapping_df_after,"sample_name_mapping_scatter_after.csv",row.names = FALSE)
  new_id_list_after <- mapping_df_after$new_ID
  
  #复制绘图副本，原始对象不改动
  mmRNAList_after_plot <- mmRNAList_after_raw
  for(i in seq_along(mmRNAList_after_plot)){
    obj <- mmRNAList_after_plot[[i]]
    obj$orig.ident <- new_id_list_after[i]
    Idents(obj) <- "orig.ident"
    mmRNAList_after_plot[[i]] <- obj
  }
  
  #循环绘图：单样本图，增加R值文本标注
  for(i in seq_along(mmRNAList_after_plot)){
    obj <- mmRNAList_after_plot[[i]]
    plot_df <- obj@meta.data
    
    #计算Pearson相关系数R
    cor_res <- cor.test(plot_df$nCount_RNA, plot_df$nFeature_RNA, method = "pearson")
    r_val <- round(cor_res$estimate, 3)
    label_text <- paste0("R = ", r_val)
    
    p <- ggplot(plot_df, aes(x = nCount_RNA, y = nFeature_RNA)) +
      geom_point(aes(color = mt_percent), size = 0.4, alpha = 0.6) +
      scale_color_viridis_c(name = "mt%", option = "viridis") +
      annotate("text", x = Inf, y = Inf, label = label_text, 
               hjust = 1.1, vjust = 1.1, size = 4.5) +
      labs(x = "nCount_RNA", y = "nFeature_RNA",
           title = paste0(new_id_list_after[i])) +
      theme_bw(base_size =11)+
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    ggsave(paste0("scatter_nCount_nFeature_mt_",new_id_list_after[i],"_after.pdf"),
           plot = p, width =7, height =6, dpi=300)
    cat("✅已输出：",mapping_df_after$original_ID[i]," -> ",new_id_list_after[i]," | ",label_text,"\n")
  }
  cat("====全部单样本nCount‑nFeature线粒体着色散点图绘制完成====\n")
  
  #====================所有样本汇总总图，同样计算并标注R值====================
  seu_merge_scatter_after <- merge(mmRNAList_after_plot[[1]], y = mmRNAList_after_plot[-1])
  df_all_meta_after <- seu_merge_scatter_after@meta.data
  
  cor_all_after <- cor.test(df_all_meta_after$nCount_RNA, df_all_meta_after$nFeature_RNA, method = "pearson")
  r_all_val_after <- round(cor_all_after$estimate,3)
  label_all_text_after <- paste0("R = ", r_all_val_after)
  
  p_scatter_all_after <- ggplot(df_all_meta_after, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(aes(color = mt_percent), size = 0.12, alpha = 0.35) +
    scale_color_viridis_c(name = "mt%", option = "viridis") +
    annotate("text", x = Inf, y = Inf, label = label_all_text_after, 
             hjust = 1.1, vjust = 1.1, size =4.5) +
    labs(x = "nCount_RNA",
         y = "nFeature_RNA",
         title = "All samples: nCount_RNA vs nFeature_RNA (after filtering)") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  ggsave("scatter_allSample_nCount_nFeature_mt_after.pdf",
         plot = p_scatter_all_after, width = 10, height = 8, dpi = 300)
  cat("✅全部样本汇总散点图已输出：scatter_allSample_nCount_nFeature_mt_after.pdf | ",label_all_text_after,"\n")
  
}





#####5.3 柱状图：18例样本 质控过滤前后细胞数目对比####

{
  # 1.读取过滤前、过滤后list
  mmRNAList_front <- readRDS("mmRNAList_filter_front.rds")
  mmRNAList_after <- readRDS("mmRNAList_filter_after.rds")
  
  # 提取原始样本ID（两份rds样本顺序必须完全一致）
  raw_id_front <- sapply(mmRNAList_front, function(x) unique(x$orig.ident))
  raw_id_after <- sapply(mmRNAList_after, function(x) unique(x$orig.ident))
  
  # 校验两份样本顺序是否一致，不一致会警告
  if(!identical(raw_id_front, raw_id_after)){
    stop("⚠️警告：过滤前与过滤后样本顺序不一致，不能继续绘图！")
  }
  
  # 命名规则，和前面脚本完全统一
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  
  group_tag_vec <- sapply(raw_id_front, get_group_tag)
  
  mapping_df <- data.frame(
    original_ID = raw_id_front,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = str_pad(row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  # 统计细胞数
  cell_count_df <- data.frame(
    sample_new = mapping_df$new_ID,
    original_ID = mapping_df$original_ID,
    count_before = sapply(mmRNAList_front, ncol), # seurat对象ncol为细胞数
    count_after  = sapply(mmRNAList_after, ncol)
  )
  
  # 计算过滤丢掉的细胞数量
  cell_count_df$count_remove <- cell_count_df$count_before - cell_count_df$count_after
  
  # 输出细胞统计表格
  write.csv(cell_count_df, "sample_cellcount_before_after.csv", row.names = FALSE)
  cat("✅细胞统计表格输出：sample_cellcount_before_after.csv\n")
  print(cell_count_df)
  
  # 转为长表格，ggplot分组柱状图需要long format
  cell_count_long <- cell_count_df %>%
    tidyr::pivot_longer(cols = c(count_before, count_after),
                        names_to = "filter_status",
                        values_to = "cell_number") %>%
    mutate(filter_status = factor(filter_status,
                                  levels = c("count_before","count_after"),
                                  labels = c("Before filter","After filter")))
  
  # 绘制分组柱状图
  p_cell_bar <- ggplot(cell_count_long, aes(x = sample_new, y = cell_number, fill = filter_status)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("#4472C4","#ED7D31")) +
    labs(x = "Sample", y = "Cell number", fill = "Filter status",
         title = "Cell counts before and after QC filtering") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size =9),
      plot.title = element_text(hjust = 0.5),
      panel.grid.x = element_blank()
    )
  
  ggsave("barplot_cellcount_before_after.pdf", plot = p_cell_bar, width = 16, height = 6.5, dpi =300)
  cat("✅柱状图已输出：barplot_cellcount_before_after.pdf\n")
}





#####6、merge合并样本####
{
#====读取数据，计算=========================================================
mmRNAList <- readRDS("mmRNAList_filter_after.rds")

mmRNAList <- merge(x=mmRNAList[[1]],y=mmRNAList[-1])

mmRNAList <- JoinLayers(mmRNAList)  # V5版的数据结构需要有这一步

# 保存为RDS（R专用格式，保持所有属性）
saveRDS(mmRNAList, file = "mmRNAList_merge.rds")

#====读取数据，绘图=========================================================
mmRNAList <- readRDS("mmRNAList_merge.rds")
## 统计细胞数
table(mmRNAList[[]]$orig.ident)

# 绘图
violin_after <- VlnPlot(mmRNAList,
                        features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                        pt.size = 0.01,
                        ncol = 4)
# 将图片输出到画板上 
violin_after
# 保存图片
ggsave("vlnplot_after_qc.pdf", plot = violin_after, width = 15, height =7) 

}



#####7、数据归一化、筛选高变基因与PCA降维（可加载数据）####
{
#====读取数据，计算=========================================================
# 重新加载融合后、归一化前的数据（可直接执行这行，读取既往数据结果，免得再执行前述读取操作）
mmRNAList <- readRDS("mmRNAList_merge.rds")

# harmony整合是基于PCA降维结果进行的。
mmRNAList <- NormalizeData(mmRNAList) %>% # 数据归一化处理
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% #筛选高变基因
  ScaleData() %>% #数据标准化
  RunPCA(npcs = 30, verbose = T)#npcs：计算和存储的PC数（默认为 50）

# 保存为RDS（R专用格式，保持所有属性）
saveRDS(mmRNAList, file = "mmRNAList_merge_PCA.rds")

#====读取数据，绘图=========================================================
mmRNAList <- readRDS("mmRNAList_merge_PCA.rds")

a=DimPlot(mmRNAList,reduction = "pca",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
a

##看下高变基因有哪些,可视化
# 提取前15个高变基因ID
top15 <- head(VariableFeatures(mmRNAList), 15) 
plot1 <- VariableFeaturePlot(mmRNAList) 
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
# 合并图片
feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
feat_15
# 保存图片
ggsave(file = "feat_15.pdf",plot = feat_15,he = 10,wi = 15 )
}


#####8、细胞周期评分####
{
#====读取数据，计算=========================================================
mmRNAList <- readRDS("mmRNAList_merge_PCA.rds")

# 提取g2m特征向量
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(mmRNAList))
# 提取s期特征向量
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(mmRNAList))
# 对细胞周期阶段进行评分

mmRNAList <- CellCycleScoring(object = mmRNAList, 
                           s.features = s_genes, 
                           g2m.features = g2m_genes, 
                           set.ident = TRUE)#set.ident 是否给每个细胞标注一个细胞周期标签
# 清理临时变量释放内存
rm(g2m_genes, s_genes)
gc(verbose = F)

# 保存为RDS（R专用格式，保持所有属性）
saveRDS(mmRNAList, file = "mmRNAList_merge_PCA_g2m.rds")

#====读取数据，绘图=========================================================
mmRNAList <- readRDS("mmRNAList_merge_PCA_g2m.rds")

mmRNAList@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()

}


#####9、RunHarmony去批次####
{
#====读取数据，计算=========================================================
mmRNAList <- readRDS("mmRNAList_merge_PCA_g2m.rds")

# 整合需要指定Seurat对象和metadata中需要整合的变量名。
scRNA_harmony <- RunHarmony(mmRNAList, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]

# 保存为RDS（R专用格式，保持所有属性）
saveRDS(scRNA_harmony, file = "scRNA_harmony.rds")

#====读取数据，绘图=========================================================
scRNA_harmony <- readRDS("scRNA_harmony.rds")

b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
b

}


#####10、聚类、umap/tsne降维（可加载数据）####
{
#====读取数据，计算=========================================================
# 重新加载（可直接执行这行，读取既往数据结果，免得再执行前述读取操作）
scRNA_harmony <- readRDS("scRNA_harmony.rds")

# 提取全部元数据列名
meta_cols <- colnames(scRNA_harmony@meta.data)
# 筛选出所有 RNA_snn_res. 开头的分辨率列
res_cols <- meta_cols[grepl("^RNA_snn_res\\.", meta_cols)]
# 删除这些列（旧精度聚类列）
scRNA_harmony@meta.data[, res_cols] <- NULL

ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony") #Harmony降维后的主成分（PC）数量选择图，一般选择拐点数量作为后续聚类dims

# 先基于最优PC（之前elbow图拐点1:15）构建近邻图，只运行一次
scRNA_harmony <- FindNeighbors(
  object = scRNA_harmony,
  reduction = "harmony",
  dims = 1:15
)

# 批量循环跑分辨率 0.1 ~ 1.2，步长0.1，全部存入对象
res_range <- seq(from = 0.1, to = 1.2, by = 0.1)
scRNA_harmony <- FindClusters(
  object = scRNA_harmony,
  resolution = res_range,
  verbose = TRUE
)

# 保存带多分辨率聚类信息的对象（避免重复计算）
saveRDS(scRNA_harmony, "scRNA_harmony_multiRes.rds")

scRNA_harmony <- readRDS("scRNA_harmony_multiRes.rds")
# 绘制不同精度的聚类树
library(clustree)
tree_plot <- clustree(scRNA_harmony, prefix = "RNA_snn_res.")
print(tree_plot)
ggsave("clustree_resolution_tree.pdf", plot = tree_plot, width = 22, height = 20, dpi = 300)


# 根据前述聚类树，固定最优分辨率为0.4，在不使用seurat_clusters时的精度可由这个去固定，但重要的是下面那个
Idents(scRNA_harmony) <- "RNA_snn_res.0.4"#后续寻找差异基因要用这个

table(scRNA_harmony@meta.data$seurat_clusters)#这个元数据里的seurat_cluster默认保存是最后一个精度值的簇，比如此时是最大的1.2精度对应44类，需要用指定精度列的簇数值覆盖，这样不影响后续使用想选的精度。

# 重要：彻底固定最优分辨率，将指定精度0.4分辨率的分群，直接覆盖到scRNA_harmony@meta.data$seurat_clusters列
scRNA_harmony$seurat_clusters <- scRNA_harmony$RNA_snn_res.0.4 ##注意这里精度一定要改为自己想用的精度值

# 验证：现在table出来就是26个簇，不再是44
table(scRNA_harmony$seurat_clusters)

##再运行umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:15)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:15)

# 保存为RDS（R专用格式，保持所有属性）
saveRDS(scRNA_harmony, file = "scRNA_harmony_umap_tsne.rds")

#====读取数据，绘图=========================================================
scRNA_harmony <- readRDS("scRNA_harmony_umap_tsne.rds")

# 按样本绘图
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
# 合并图片
umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
# 将图片输出到画板
umap_tsne_integrated
# 保存图片
ggsave("umap_tsne_integrated.pdf",umap_tsne_integrated,wi=25,he=15)

# 按分组去绘图
umap_integrated1_group <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "group")
tsne_integrated1_group <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "group") 
# 合并图片
umap_tsne_integrated_group <- CombinePlots(list(tsne_integrated1_group,tsne_integrated2,umap_integrated1_group,umap_integrated2),ncol=2)
# 将图片输出到画板
umap_tsne_integrated_group
# 保存图片
ggsave("umap_tsne_integrated_group.pdf",umap_tsne_integrated_group,wi=25,he=15)

table(scRNA_harmony@meta.data$seurat_clusters)#这个元数据里的seurat_cluster默认保存是最后一个精度值的簇，比如此时是最大的1.2精度对应44类，需要用指定精度列的簇数值覆盖，这样不影响后续使用想选的精度。


# 1. 提取患者组细胞ID
patient_cells <- rownames(scRNA_harmony@meta.data)[scRNA_harmony$group == "patient"]

# 2. 单独绘制患者组TSNE（保留你原有参数：label=T、repel=TRUE、pt.size=1）
p.dim.cell.patient.tsne <- DimPlot(
  scRNA_harmony, 
  reduction = "tsne", 
  group.by = "seurat_clusters",
  label = T, 
  repel = TRUE,
  pt.size = 1,
  cells = patient_cells  # 关键参数：仅绘制患者组细胞
) 
p.dim.cell.patient.tsne

# 3. 单独绘制患者组UMAP（同逻辑）
p.dim.cell.patient.umap <- DimPlot(
  scRNA_harmony, 
  reduction = "umap", 
  group.by = "seurat_clusters",
  label = T, 
  repel = TRUE,
  pt.size = 1,
  cells = patient_cells
) 
p.dim.cell.patient.umap


# 1. 提取对照组细胞ID
control_cells <- rownames(scRNA_harmony@meta.data)[scRNA_harmony$group == "control"]

# 2. 单独绘制对照组TSNE
p.dim.cell.control.tsne <- DimPlot(
  scRNA_harmony, 
  reduction = "tsne", 
  group.by = "seurat_clusters",
  label = T, 
  repel = TRUE,
  pt.size = 1,
  cells = control_cells  # 关键参数：仅绘制对照组细胞
) 
p.dim.cell.control.tsne

# 3. 单独绘制对照组UMAP
p.dim.cell.control.umap <- DimPlot(
  scRNA_harmony, 
  reduction = "umap", 
  group.by = "seurat_clusters",
  label = T, 
  repel = TRUE,
  pt.size = 1,
  cells = control_cells
) 
p.dim.cell.control.umap


##绘制标记基因库里的基因在各类中的表达情况图

# === 1 自定义Marker分组 ===
marker_list <- list(
  "T cell" = c("CD3D", "CD3E", "CD3G", "CD40LG", "CD8A", "CD8B"),
  "B cell" = c("CD79A", "CD79B", "MS4A1"),
  "NK cell" = c("GNLY", "NKG7", "TYROBP"),
  "Monocyte" = c("CD14", "FCN1", "S100A8", "S100A9", "VCAN"),
  "cDC1" = c("CLEC9A", "XCR1", "BATF3"),
  "cDC2" = c("CD1C", "FCER1A", "CLEC10A"),
  "pDC" = c("LILRA4", "IL3RA", "CLEC4C", "IRF7"),
  "Macrophage" = c("CD68", "CD163"),
  #"Neutrophil" = c("FCGR3B", "CSF3R"), #不含中性粒就注释掉它
  "Megakaryocyte" = c("PF4", "PPBP"),
  "Mast cell" = c("KIT", "CPA3"),
  "Epithelial cell" = c("KRT18", "KRT19")
)

# === 2 全局可调参数 ===
seurat_obj <- scRNA_harmony  # 替换成你的Seurat对象
reduction_type <- "umap"       # tsne / umap 切换
pt_size <- 0.28                # 细胞点稍大一点更清晰
max_col <- 6                   # 每行固定6张子图
# 加深配色：浅灰底色 + 深酒红（对比度更强）
color_vec <- c("lightgray", "#990000")

# === 3 空白占位图函数 ===
blank_plot <- function(){
  ggplot() + 
    theme_void() + 
    theme(plot.background = element_rect(fill="transparent", color=NA))
}

# === 4 批量绘图循环 ===
all_rows <- list()
for(ct in names(marker_list)){
  genes_raw <- marker_list[[ct]]
  genes_use <- intersect(genes_raw, rownames(seurat_obj))
  
  # 绘制单基因FeaturePlot，新增raster=FALSE关闭栅格
  sub_plots <- lapply(genes_use, function(g){
    FeaturePlot(
      seurat_obj,
      features = g,
      reduction = reduction_type,
      cols = color_vec,
      pt.size = pt_size,
      order = TRUE
      #raster = FALSE  # 可关闭栅格，消除提示，但绘制太慢
    ) +
      labs(title = g) +
      theme(
        plot.title = element_text(hjust=0.5, size=11),
        axis.title = element_text(size=9),
        axis.text = element_text(size=7),
        legend.key.height = unit(0.8, "cm")
      )
  })
  
  # 空白补齐到6列
  need_blank <- max_col - length(sub_plots)
  if(need_blank > 0){
    blank_list <- rep(list(blank_plot()), need_blank)
    sub_plots <- c(sub_plots, blank_list)
  }
  
  # 横向拼接一行，修复bold引号问题
  row_p <- wrap_plots(sub_plots, nrow = 1, heights = 3.2)
  row_p <- row_p + plot_annotation(
    title = paste0("Cell type: ", ct),
    theme = theme(plot.title = element_text(size = 14, face = "bold")) # 加双引号
  )
  all_rows[[ct]] <- row_p
}


# 所有行纵向拼接
final_figure <- wrap_plots(all_rows, ncol = 1)

# ===5 导出（加宽加高画布） ===
ggsave(
  "cell_marker_deepcolor_tall.png",
  final_figure,
  width = 18,
  height = 24,
  dpi = 400,
  device = "png"
)

print(final_figure)

}


#====10.1 读取降维数据，绘图按新的样本名称编号=========================================================

{
  scRNA_harmony <- readRDS("scRNA_harmony_umap_tsne.rds")
  
  # 1.复制中间绘图副本，原始对象保留不变
  scRNA_plot <- scRNA_harmony
  
  # 2.提取全部原始样本orig.ident
  raw_id_all <- as.character(unique(scRNA_plot$orig.ident))
  
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  
  group_tag_vec <- sapply(raw_id_all, get_group_tag)
  group_tag_vec <- as.character(group_tag_vec)
  
  umap_mapping_df <- data.frame(
    original_ID = raw_id_all,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = str_pad(row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  write.csv(umap_mapping_df, "umap_sample_name_mapping.csv", row.names = FALSE)
  cat("✅UMAP样本名称映射表已保存 umap_sample_name_mapping.csv\n")
  print(umap_mapping_df[,c("original_ID","new_ID")])
  
  # =========【重点：删除left_join，改用match直接赋值，避开tibble join坑】========
  # match：按orig.ident匹配，直接生成整细胞向量
  match_index <- match(scRNA_plot@meta.data$orig.ident, umap_mapping_df$original_ID)
  # 直接新增列写入meta.data，Seurat $符号可以正常识别
  scRNA_plot@meta.data$new_ID <- umap_mapping_df$new_ID[match_index]
  
  # 校验
  na_count <- sum(is.na(scRNA_plot@meta.data$new_ID))
  cat("⚠️meta.data中new_ID为NA的细胞数：", na_count,"\n")
  
  # 这里必须验证：scRNA_plot$new_ID，Seurat对象$符号访问必须有值，不能NA
  head(scRNA_plot$new_ID)
  
  # 绘图
  umap_integrated1 <- DimPlot(scRNA_plot,
                              reduction = "umap",
                              group.by = "new_ID") +
    ggtitle("UMAP: colored by sample") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.4,"cm"))
  
  ggsave("umap_groupBy_sample_newID.pdf", plot = umap_integrated1, width =12, height =9, dpi=300)
  cat("✅UMAP图已输出 umap_groupBy_sample_newID.pdf\n")
  
  tsne_integrated1 <- DimPlot(scRNA_plot,
                              reduction = "tsne",
                              group.by = "new_ID") +
    ggtitle("TSNE: colored by sample") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.4,"cm"))
  
  ggsave("tsne_groupBy_sample_newID.pdf", plot = tsne_integrated1, width =12, height =9, dpi=300)
  cat("✅TSNE图已输出 tsne_groupBy_sample_newID.pdf\n")

}






#####11、FindAllMarkers（可加载数据）####
{
  # 关键：先加载Seurat，必须放在最前面！
  library(Seurat)
  library(dplyr)
  library(ggplot2)

  
  #====读取数据，计算=========================================================
  scRNA_harmony <- readRDS("scRNA_harmony_umap_tsne.rds")
  
  # 原FindAllMarkers一次性计算全部cluster内存爆炸，替换为循环逐个计算，大幅降低内存占用
  #分别对每个cluster进行与剩下除他之外所有的cluster的差异分析
  all_clusters <- sort(unique(scRNA_harmony$seurat_clusters))
  marker_list <- list()
  
  for (cur_clu in all_clusters) {
    cat("正在计算 cluster", cur_clu, "\n")
    cur_marker <- FindMarkers(
      object = scRNA_harmony,
      ident.1 = cur_clu,
      test.use = "wilcox",
      only.pos = TRUE,
      logfc.threshold = 0.25,
      min.pct = 0.1,        # 过滤低表达基因，减少计算量
      min.diff.pct = 0.1    # 过滤表达差异极小基因
    )
    cur_marker$cluster <- cur_clu
    cur_marker$gene <- rownames(cur_marker)
    marker_list[[as.character(cur_clu)]] <- cur_marker
    gc() # 每算完一个分群强制释放临时内存
  }
  # 合并所有分群marker结果，等效FindAllMarkers输出格式
  markers <- bind_rows(marker_list)
  
  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(markers, file = "markers_01.rds")
  markers <- readRDS("markers_01.rds")
  
  # 容错判断：如果无marker基因直接终止运行，避免后续gene列不存在报错
  if (nrow(markers) == 0) {
    stop("未检测到任何差异Marker基因！请调高内存或降低logfc/min.pct筛选阈值")
  }
  
  # 对计算好的每cluster的marker基因进行筛选
  all.markers = markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj < 0.05)
  #筛选出P<0.05的marker基因
  top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) #将每个cluster lgFC排在前15的marker基因挑选出来
  View(top15)
  
  # row.names=F 去除多余行号列，更适合后续导入Excel
  write.csv(top15, "cluster_top15.csv", row.names = F)
  write.csv(all.markers, "cluster_allmarkers.csv", row.names = F)
  write.csv(markers, "cluster_markers.csv", row.names = F)
  
  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(top15, file = "scRNA_harmony_top15_markers.rds")
  
  # 保存为RDS（R专用格式，保持所有属性）
  saveRDS(scRNA_harmony, file = "scRNA_harmony_markers.rds")
  
  #====读取数据，绘图=========================================================
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  
  scRNA_harmony <- readRDS("scRNA_harmony_markers.rds")
  top15 <- readRDS("scRNA_harmony_top15_markers.rds")
  
  # 绘图内存优化：随机抽取4000细胞绘图，避免全量细胞内存溢出
  set.seed(123)
  sample_cells <- sample(colnames(scRNA_harmony), size = 4000)
  sc_plot <- scRNA_harmony[, sample_cells]
  
  # 移除 downsample = 50，DoHeatmap不支持该参数
  DoHeatmap(sc_plot, features = top15$gene, slot = "data") + NoLegend()
  #slot默认使用scaledata里边只有2k个高变gene表达 这里要使用data数据，不然有些gene会找不到表达量
  
  #选前15个makergene看看,正好前15行就是cluster0的差异基因
  VlnPlot(sc_plot, features = top15$gene[1:15], ncol = 5)
  
  #先看下前15行基因的点图，恰属于cluster0，这里可以任意改top15$gene范围或者具体的基因
  p <- DotPlot(sc_plot, features = top15$gene[1:15],
               assay = 'RNA', group.by = 'seurat_clusters') + coord_flip() + ggtitle("")
  print(p)
  
  # 或者自己选gene看看
  VlnPlot(sc_plot, features = c("CD3D", "CD3E"))
  
  # ggsave(filename="plot.pdf",width = 210,height = 297,units = "mm")
  
  #找到每一个cluster的celltype后，把每一个cluster命名成对应的celltype
  View(scRNA_harmony@meta.data)#seurat_clusters这一列存放了每个细胞对应的cluster
  table(scRNA_harmony@meta.data$seurat_clusters)#可以查看每个cluster有多少个细胞
  
}




#####12、进行细胞注释####
{
#====读取数据，计算=========================================================
# 重新加载（可直接执行这行，读取既往数据结果，免得再执行前述读取操作）
scRNA_harmony <- readRDS("scRNA_harmony_markers.rds")

library(ggplot2) 

#利用整理的一些markers，例如制作的外周血makers文件：blood_cell_markers_major.rds
# 使用marker基因库进行细胞类型注释
blood_cell_markers <- readRDS("blood_cell_markers_ppt_standard.rds")  #加载提前自行生成的基因库

RNA_harmony_annotation <- scRNA_harmony #重命名要注释的seurat对象，以跟之前未注释状态作区分

## 自动细胞注释函数（deepseek生成）
# 基于cluster的细胞类型注释（每个cluster一个类型）
# @param seurat_obj , Seurat对象
# @param markers_list , marker基因列表
# @return 包含cluster_celltypes列的Seurat对象
  # 最简单直接的解决方案

simple_cluster_annotation <- function(seurat_obj, markers_list,
                                      pct_cut = 0.15,    # 基因阳性占比阈值，低于不计分，下调阳性阈值，适配 DC 稀疏表达
                                      expr_cut = 0.1,   # T谱系表达阈值（拦截判定）
                                      score_min = 0.08  # 最低有效总分，低于标未知
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
  
  # 3. 统一转为数值型（完全保留你原版稳定写法）
  seurat_obj@meta.data$seurat_clusters <- as.numeric(
    as.character(seurat_obj@meta.data$seurat_clusters)
  )
  
  # 4. 提取有效cluster（单列安全筛选，无维度错位）
  clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters[!is.na(seurat_obj@meta.data$seurat_clusters)]))
  cluster_celltypes <- list()
  
  # 5. 权重规则（彻底修正CD40LG重复bug，所有类型格式统一）
  weight_rule <- list(
    "T cell" = list(core = c("CD3D","CD3E","CD3G","CD8A","CD8B"), aux = c("CD40LG"), w_core = 3, w_aux = 1),
    "B cell" = list(core = c("CD79A","CD79B","MS4A1"), aux = c(), w_core = 3, w_aux = 1),
    
    "NK cell" = list(core = c("GNLY","NKG7","TYROBP"), aux = c(), w_core = 1, w_aux = 1),
    # 单核权重从3下调至1.8，削弱碾压效果，以分离出cDC
    "Monocyte" = list(core = c("CD14","FCN1","S100A8","S100A9"), aux = c("VCAN"), w_core = 1.8, w_aux = 1),
    # DC维持3权重，提升竞争力
    "cDC1" = list(core = c("CLEC9A","XCR1","BATF3"), aux = c(), w_core = 3, w_aux = 1),
    "cDC2" = list(core = c("CD1C","FCER1A","CLEC10A"), aux = c(), w_core = 3, w_aux = 1),
    "pDC" = list(core = c("LILRA4","IL3RA","CLEC4C","IRF7"), aux = c(), w_core = 3, w_aux = 1),
    "Macrophage" = list(core = c("CD68","CD163"), aux = c(), w_core = 3, w_aux = 1),
    "Megakaryocyte" = list(core = c("PF4","PPBP"), aux = c(), w_core = 3, w_aux = 1),
    "Mast cell" = list(core = c("KIT","CPA3"), aux = c(), w_core = 3, w_aux = 1),
    "Epithelial cell" = list(core = c("KRT18","KRT19"), aux = c(), w_core = 3, w_aux = 1)
  )
  
  # 新增：打印当前全套权重配置，方便核对
  cat("\n================ 当前权重配置一览 ================\n")
  for (cell_name in names(weight_rule)) {
    rule <- weight_rule[[cell_name]]
    cat(sprintf("【%s】\n", cell_name))
    cat(sprintf("  core基因：%s\n", paste0(rule$core, collapse = ",")))
    cat(sprintf("  aux 基因：%s\n", ifelse(length(rule$aux)==0, "无", paste0(rule$aux, collapse = ","))))
    cat(sprintf("  核心权重w_core = %.2f , 其余权重w_aux = %.2f\n\n", rule$w_core, rule$w_aux))
  }
  cat("==================================================\n\n")
  
  
  # 6. 强制校验marker列表与权重列表细胞类型完全匹配，防止NULL报错
  if (!identical(sort(names(markers_list)), sort(names(weight_rule)))) {
    stop("错误：markers_list 与 weight_rule 的细胞类型名称不匹配，请核对！")
  }
  
  cat("========= 优化版细胞注释启动 =========\n")
  
  # 7. 逐簇循环
  for (cluster in clusters) {
    cell_idx <- which(seurat_obj@meta.data$seurat_clusters == cluster)
    cell_bar <- colnames(seurat_obj)[cell_idx]
    n_cell <- length(cell_bar)
    t_expr_mean <- 0 # 提前初始化，消除作用域未定义隐患
    
    # 空簇跳过
    if (n_cell == 0) {
      cluster_celltypes[[as.character(cluster)]] <- "空Cluster"
      cat(sprintf("Cluster %d：空Cluster，细胞数0\n", cluster))
      next
    }
    
    # 循环内临时读取表达矩阵，不全局占用大内存
    expr_all <- Seurat::GetAssayData(seurat_obj, layer = "data", assay = "RNA")
    all_gene <- rownames(expr_all)
    expr_cl <- expr_all[, cell_bar, drop = F]
    rm(expr_all)
    gc()
    
    
    
    # --- 第二层：加权打分计算所有细胞类型 ---
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
      # 过滤阳性率不达标的核心基因
      if (length(core_g) > 0) {
        for (g in core_g) {
          p <- sum(expr_cl[g,] > 0) / n_cell
          if (p >= pct_cut) valid_core <- c(valid_core, g)
        }
      }
      # 过滤阳性率不达标的辅助基因（已修复：追加g，不再空向量）
      if (length(aux_g) > 0) {
        for (g in aux_g) {
          p <- sum(expr_cl[g,] > 0) / n_cell
          if (p >= pct_cut) valid_aux <- c(valid_aux, g)
        }
      }
      
      # 计算加权分数
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
      
      # NK惩罚系数：存在T信号分数压缩，原 *0.3 改为 *0.7
      if (ctype == "NK cell" && t_expr_mean >= expr_cut) {
        total_score <- total_score * 0.7
      }
      score_list[ctype] <- total_score
    }
    
    # --- 第三层：判定最终细胞类型 ---
    max_score <- max(score_list)
    if (max_score < score_min) {
      final_ct <- "混杂/污染细胞"
    } else {
      final_ct <- names(which.max(score_list))
    }
    cluster_celltypes[[as.character(cluster)]] <- final_ct
    cat(sprintf("Cluster %d | %s | 最高分：%.3f | 细胞数：%d\n",
                cluster, final_ct, max_score, n_cell))
    
    # 释放当前簇表达矩阵，持续控内存
    rm(expr_cl)
    gc()
  }
  
  # 构建映射表，批量赋值celltype
  map_df <- data.frame(
    cluster = as.numeric(names(cluster_celltypes)),
    celltype = unlist(cluster_celltypes),
    stringsAsFactors = F
  )
  # 显式调用plyr，规避dplyr包名冲突
  seurat_obj$celltype <- plyr::mapvalues(
    x = seurat_obj$seurat_clusters,
    from = map_df$cluster,
    to = map_df$celltype,
    warn_missing = F
  )
  # 无匹配细胞兜底标签
  seurat_obj$celltype[is.na(seurat_obj$celltype)] <- "混杂/污染细胞"
  
  # 输出汇总统计
  cat("\n========= 注释汇总结果 =========\n")
  for (i in 1:nrow(map_df)) {
    cl <- map_df$cluster[i]
    ct <- map_df$celltype[i]
    cnt <- sum(seurat_obj$seurat_clusters == cl, na.rm=T)
    cat(sprintf("Cluster %s：%s，共%d细胞\n", cl, ct, cnt))
  }
  
  return(seurat_obj)
}




# 使用示例,运行这个简单版本进行自动注释
RNA_harmony_annotation <- simple_cluster_annotation(
  seurat_obj = RNA_harmony_annotation,
  markers_list = blood_cell_markers
)

table(RNA_harmony_annotation@meta.data$celltype,RNA_harmony_annotation@meta.data$seurat_clusters)

# 保存为RDS（R专用格式，保持所有属性）
saveRDS(RNA_harmony_annotation, file = "RNA_harmony_annotation.rds")

#====读取数据，绘图=========================================================
RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")

##注释结果绘图TSNE，并保存文件
p.dim.cell=DimPlot(RNA_harmony_annotation, reduction = "tsne", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
p.dim.cell
ggsave(plot=p.dim.cell,filename="DimPlot_tsne_celltype.pdf",width=9, height=7) #PDF无法加载中文名称，故细胞类型名可改为英文，或者保存图片格式而非PDF

##注释结果绘图UMAP，并保存文件
p.dim.cell=DimPlot(RNA_harmony_annotation, reduction = "umap", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
p.dim.cell
ggsave(plot=p.dim.cell,filename="DimPlot_umap_celltype.pdf",width=9, height=7) #PDF无法加载中文名称，故细胞类型名可改为英文，或者保存图片格式而非PDF



#生成表格比较细胞比例
# 构建交叉计数矩阵
ct <- table(RNA_harmony_annotation$orig.ident, RNA_harmony_annotation$celltype)

# 转成数据框，行名是样本
df_wide <- as.data.frame.matrix(ct)
# 新增样本ID列放到最前面
df_wide$sample_id <- rownames(df_wide)
df_wide <- df_wide[, c("sample_id", colnames(df_wide)[1:(ncol(df_wide)-1)])]

# 输出
write.csv(df_wide, file.path("D:/数据分析/单细胞数据分析/FTD_data_202512/wd_R_V2", "cell_count_wide.csv"), row.names = F, fileEncoding = "UTF-8")

# 打印预览核对
print(head(df_wide, 3))

# 1. 提取患者组细胞ID
patient_cells <- rownames(RNA_harmony_annotation@meta.data)[RNA_harmony_annotation$group == "patient"]


# 2. 单独绘制患者组TSNE（保留你原有参数：label=T、repel=TRUE、pt.size=1）
p.dim.cell.patient.tsne <- DimPlot(
  RNA_harmony_annotation, 
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
ggsave(plot=p.dim.cell.patient.tsne, filename="DimPlot_tsne_celltype_patient.pdf", width=9, height=7)

# 4. 单独绘制患者组UMAP（同逻辑）
p.dim.cell.patient.umap <- DimPlot(
  RNA_harmony_annotation, 
  reduction = "umap", 
  group.by = "celltype",
  label = T, 
  repel = TRUE,
  pt.size = 1,
  raster = TRUE,
  cells = patient_cells
) 
p.dim.cell.patient.umap
ggsave(plot=p.dim.cell.patient.umap, filename="DimPlot_umap_celltype_patient.pdf", width=9, height=7)


# 1. 提取对照组细胞ID
control_cells <- rownames(RNA_harmony_annotation@meta.data)[RNA_harmony_annotation$group == "control"]

# 2. 单独绘制对照组TSNE
p.dim.cell.control.tsne <- DimPlot(
  RNA_harmony_annotation, 
  reduction = "tsne", 
  group.by = "celltype",
  label = T, 
  repel = TRUE,
  pt.size = 1,
  raster = TRUE,
  cells = control_cells  # 关键参数：仅绘制对照组细胞
) 
p.dim.cell.control.tsne
ggsave(plot=p.dim.cell.control.tsne, filename="DimPlot_tsne_celltype_control.pdf", width=9, height=7)

# 3. 单独绘制对照组UMAP
p.dim.cell.control.umap <- DimPlot(
  RNA_harmony_annotation, 
  reduction = "umap", 
  group.by = "celltype",
  label = T, 
  repel = TRUE,
  pt.size = 1,
  raster = TRUE,
  cells = control_cells
) 
p.dim.cell.control.umap
ggsave(plot=p.dim.cell.control.umap, filename="DimPlot_umap_celltype_control.pdf", width=9, height=7)


# 1. 提取carrier组细胞ID
carrier_cells <- rownames(RNA_harmony_annotation@meta.data)[RNA_harmony_annotation$group == "carrier"]

# 2. 单独绘制carrier组TSNE
p.dim.cell.carrier.tsne <- DimPlot(
  RNA_harmony_annotation, 
  reduction = "tsne", 
  group.by = "celltype",
  label = T, 
  repel = TRUE,
  pt.size = 1,
  raster = TRUE,
  cells = carrier_cells  # 仅绘制carrier组细胞
) 
p.dim.cell.carrier.tsne
ggsave(plot=p.dim.cell.carrier.tsne, filename="DimPlot_tsne_celltype_carrier.pdf", width=9, height=7)

# 3. 单独绘制carrier组UMAP
p.dim.cell.carrier.umap <- DimPlot(
  RNA_harmony_annotation, 
  reduction = "umap", 
  group.by = "celltype",
  label = T, 
  repel = TRUE,
  pt.size = 1,
  raster = TRUE,
  cells = carrier_cells
) 
p.dim.cell.carrier.umap
ggsave(plot=p.dim.cell.carrier.umap, filename="DimPlot_umap_celltype_carrier.pdf", width=9, height=7)

}

#####13、cellchat流程（64G内存机器专用优化版）#####
{
  library(CellChat)
  library(magrittr)
  library(future)
  library(Matrix)
  library(dplyr)
  # ==== 全局并行内存上限（64G机器，放宽上限）===
  options(future.globals.maxSize = 40 * 1024^3) # 单进程允许40G全局对象
  options(scipen = 999)
  
####13.1、patient组cellchat分析####
{  
  #====读取数据，计算=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  
  # V5 Seurat 极致精简，删除冗余图层、降维、图网络，大幅减小对象体积
  RNA_harmony_annotation <- JoinLayers(RNA_harmony_annotation, assay = "RNA")
  RNA_harmony_annotation <- DietSeurat(
    RNA_harmony_annotation, 
    assays = "RNA", 
    dimreducs = NULL, 
    graphs = FALSE,
    reductions = NULL,
    features = NULL
  )
  
  # 提取patient分组
  stim.object <- subset(RNA_harmony_annotation, group == "patient")
  rm(RNA_harmony_annotation) # 删除原始大Seurat对象
  gc(verbose = FALSE, reset = TRUE) # 强制垃圾回收释放内存
  
  # 提取标准化表达矩阵 + 精简meta信息
  stim.data.input <- GetAssayData(stim.object, assay = "RNA", layer = "data")
  stim.meta <- stim.object@meta.data[, c("celltype", "group")]
  stim.meta$celltype %<>% as.vector()
  rm(stim.object)
  gc(verbose = FALSE, reset = TRUE)
  
  ###1、构建cellchat对象###
  stim.cellchat <- createCellChat(object = stim.data.input)
  stim.cellchat <- addMeta(stim.cellchat, meta = stim.meta)
  stim.cellchat <- setIdent(stim.cellchat, ident.use = "celltype")
  
  levels(stim.cellchat@idents)
  groupSize <- as.numeric(table(stim.cellchat@idents))
  groupSize
  
  # 人源配受体数据库
  stim.cellchat@DB <- CellChatDB.human
  dplyr::glimpse(CellChatDB.human$interaction)
  
  # === 关键优化1：只保留配受体基因，压缩矩阵，彻底减小稠密转换开销 ===
  lr_db <- CellChatDB.human$interaction
  lr_gene_list <- unique(c(lr_db$ligand, lr_db$receptor))
  overlap_gene <- intersect(rownames(stim.cellchat@data), lr_gene_list)
  # 仅保留配受体相关基因，无关基因全部剔除
  stim.cellchat@data <- stim.cellchat@data[overlap_gene, ]
  gc(verbose = FALSE, reset = TRUE)
  
  stim.cellchat <- subsetData(stim.cellchat, features = NULL)
  
  # === 并行设置：64G机器 workers=2 平衡速度与内存峰值 ===
  future::plan("multisession", workers = 2)
  
  # === 关键优化2：调整thresh/min.cells，减少高表达基因计算量，弱化sparse转稠密内存占用 ===
  # suppressWarnings屏蔽sparse->dense提示（警告不报错，仅内存提示）
  suppressWarnings({
    stim.cellchat <- identifyOverExpressedGenes(stim.cellchat)
  })
  
  # 释放并行子进程内存
  gc(verbose = FALSE, reset = TRUE)
  
  stim.cellchat <- identifyOverExpressedInteractions(stim.cellchat)
  
  # 并行任务结束，切回串行模式，关闭多进程释放内存
  plan("sequential")
  gc(verbose = FALSE, reset = TRUE)


  ###2、cellchat分析###
  

  #通过计算与每个信号通路相关的所有配体-受体相互作用的通信概率来推断信号通路水平上的通信概率。
  # 新增trim过滤小众细胞群，大幅提速，解决长时间卡死
  stim.cellchat <- computeCommunProb(stim.cellchat,raw.use=T, trim = 10)
 
  # 过滤掉小于10个细胞的胞间通讯网络，通讯中的细胞很少没有意义
  stim.cellchat <- filterCommunication(stim.cellchat, min.cells = 10)
  # 通过汇总所有相关的配体/受体，计算信号通路水平上的通信概率
  stim.cellchat <- computeCommunProbPathway(stim.cellchat)
  stim.cellchat <- aggregateNet(stim.cellchat)#计算聚合网络
  # “netP”表示推断的信号通路的细胞间通信网络
  stim.cellchat <- netAnalysis_computeCentrality(stim.cellchat, slot.name = "netP") 
  
  #数据查看/保存
  group1.net <- subsetCommunication(stim.cellchat)  ###细胞通讯结果
  write.csv(group1.net, file = "group1_net_inter_raw.useT.csv", row.names = F)

  saveRDS(stim.cellchat,"stim.cellchat.rds")

  


###3、cellchat结果可视化###
#====读取数据，绘图=========================================================
stim.cellchat <- readRDS("stim.cellchat.rds")
library(ggrepel)

groupSize <- as.numeric(table(stim.cellchat@idents)) # 各种类型细胞数量
groupSize

#互作网络
#先定义绘制函数
plot_cellchat_circle <- function(){
# 设置1行2列分屏
par(mfrow = c(1,2))
netVisual_circle(stim.cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(stim.cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
#左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。
#发出箭头的细胞表达配体，箭头指向的细胞表达受体。配体-受体对越多，线越粗。
#右图：互作的概率/强度值（强度就是概率值相加）
}
#先Rstudio窗口预览
plot_cellchat_circle()
#打开绘图设备，指定保存路径+尺寸进行保存
pdf("CellChat_interaction_network_patient.pdf", width = 14, height = 10)
plot_cellchat_circle()#绘制在pdf中
# 关闭绘图设备，文件才会生成
dev.off()




#提取互作矩阵，展示每个亚群作为source的信号传递
mat <- stim.cellchat@net$weight#提取每个细胞的互作强度
#先定义绘制函数
plot_cellchat_circle <- function(){
# 分栏布局 + 放大边距防止文字截断
par(mfrow = c(3,3), mar = c(0.5,0.5,0.5,0.5))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  p_base<-netVisual_circle(mat2, 
                   vertex.weight = groupSize,
                   arrow.width = 0.2,arrow.size = 0.1,##自行调整
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = paste0("Source: ", rownames(mat)[i]) 
                   )
}
}
#先Rstudio窗口预览
plot_cellchat_circle()
#打开绘图设备，指定保存路径+尺寸进行保存
pdf("CellChat_Source_All_Subpop_patient.pdf", width = 16, height = 14)
plot_cellchat_circle()#绘制在pdf中
# 关闭绘图设备，写入文件（必不可少）
dev.off()




#自定义特定细胞群展示到画板上
specific_cell <- "T cell"
#先定义绘制函数
plot_cellchat_circle <- function(){
# 分栏布局 + 放大边距防止文字截断
par(mfrow = c(2,2), mar = c(1,1,1,1))
cell_order <- rownames(mat)
# 找到要绘制细胞对应的行号
idx <- which(cell_order == specific_cell)
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[idx, ] <- mat[idx, ] #绘制第几类细胞
netVisual_circle(mat2, 
                 vertex.weight = groupSize,
                 arrow.width = 0.2,arrow.size = 0.1,##自行调整
                 weight.scale = T, 
                 edge.weight.max = max(mat), 
                 title.name = rownames(mat)[2])
}
#先Rstudio窗口预览
plot_cellchat_circle()
# 打开保存设备，设置宽高适配3行3张子图
pdf(paste0("CellChat_", specific_cell, "_Source_patient.pdf"), width = 16, height = 14)
plot_cellchat_circle()#绘制在pdf中
# 关闭绘图设备，写入文件（必不可少）
dev.off()



#列举所有细胞的全部通路pathway并进行保存为Excel
#层级图
stim.cellchat@netP$pathways #全部通路名称列举
group1.net <- subsetCommunication(stim.cellchat)  ###细胞通讯全量结果
# 直接导出全量结果到Excel
write.csv(group1.net, file = "patient组_所有通路_细胞通讯全量结果.csv", row.names = FALSE)



#可选择感兴趣的通路pathway进行可视化
pathway.show <- "CCL"#这里很关键，指定生物学通路，根据研究选择来解释生物学意义
levels(stim.cellchat@idents)

#先定义绘制函数
plot_cellchat_circle <- function(){
vertex.receiver = c(4,6)#选择想看的接收者细胞类别，数字是指细胞类型的排序号
netVisual_aggregate(stim.cellchat, 
                    signaling = pathway.show,
                    vertex.receiver = vertex.receiver,
                    layout = "hierarchy")
#在层次图中，实体圆和空心圆分别表示源和目标。
#线越粗，互作信号越强。
#左图中间的target是我们选定的靶细胞。
#右图是选中的靶细胞之外的另外一组放在中间看互作。
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_receiver_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



#circle plot
#先定义绘制函数
plot_cellchat_circle <- function(){
par(mfrow = c(1,1))
netVisual_aggregate(stim.cellchat, 
                    signaling = pathway.show, 
                    layout = "circle")
title(main = paste0(pathway.show, " Signaling Pathway Cell Communication Network")) # 外部添加大图标题
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_circle_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




#chord plot
#先定义绘制函数
plot_cellchat_circle <- function(){
par(mfrow = c(1,1))
netVisual_aggregate(stim.cellchat, 
                    signaling = pathway.show, 
                    layout = "chord")#图中细胞下颜色表示通路下不同的基因互配对
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_chord_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




#热图
#chord plot
#先定义绘制函数
plot_cellchat_circle <- function(){
netVisual_heatmap(stim.cellchat, signaling = pathway.show, color.heatmap = "Reds")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_heatmap_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()





#计算受配体对整个信号通路的贡献，并可视化单个配体-受体对介导的细胞-细胞通讯
netAnalysis_contribution(stim.cellchat, signaling = pathway.show)
#可视化由单个配体-受体对介导的细胞间通讯
pairLR.CCL <- extractEnrichedLR(stim.cellchat, signaling = pathway.show, geneLR.return = FALSE)
#提取对这个通路显著性最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.CCL[1,] #选择配受体进行绘制
#先定义绘制函数
plot_cellchat_circle <- function(){
netVisual_individual(stim.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(stim.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "chord")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_",LR.show,"_contribution_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




levels(stim.cellchat@idents)
#指定受体细胞和配体细胞
#这里选择特定细胞通用，查看特定细胞的配体和其他细胞的受体
specific_cell_ligand_receptor <- "T cell"
#先定义绘制函数
plot_cellchat_circle <- function(){
p =netVisual_bubble(stim.cellchat, 
                 sources.use = specific_cell_ligand_receptor,#这里改细胞
                 remove.isolate = FALSE,
                 font.size=14)
p

}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_",specific_cell_ligand_receptor,"_ligand_receptor_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()






##参与某条信号通路的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）
#先定义绘制函数
plot_cellchat_circle <- function(){
plotGeneExpression(stim.cellchat, signaling = pathway.show) #选择信号通路
library(RColorBrewer)
# 定义颜色渐变
colors <- brewer.pal(8, "Oranges")
plotGeneExpression(stim.cellchat, signaling = pathway.show, type = "dot", col = colors)
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_",pathway.show,"_pathway_gene_expression_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




#某条通路的发送者（源）和接收者（目标）情况
#计算和可视化网络中心性评分
stim.cellchat <- netAnalysis_computeCentrality(stim.cellchat, slot.name = "netP")
#先定义绘制函数
plot_cellchat_circle <- function(){
#主要注意sender和receiver
netAnalysis_signalingRole_network(stim.cellchat, signaling = pathway.show, 
                                  width = 15, height = 6, font.size = 10)
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_",pathway.show,"_pathway_send_receiver_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




#各类细胞所有信号通路的输入输出强度点图
#先定义绘制函数
plot_cellchat_circle <- function(){
#使用散点图在 2D 空间中可视化主要的发送者（源）和接收者（目标）。
netAnalysis_signalingRole_scatter(stim.cellchat)###ALL
#netAnalysis_signalingRole_scatter(stim.cellchat, signaling = c("CXCL", "CCL"))###亦可指定通路
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_income_outgo_strength_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()
  
  

#各类细胞所有信号通路的输入信号模式图
#先定义绘制函数
plot_cellchat_circle <- function(){
  #识别对某些细胞类群的输出和输入信号贡献最大的信号
  netAnalysis_signalingRole_heatmap(stim.cellchat, pattern = "incoming")
  #netAnalysis_signalingRole_heatmap(stim.cellchat, signaling = c("CXCL", "CCL"),pattern = "incoming")###亦可指定通路
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_incomeing_signal_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




#各类细胞所有信号通路的输出信号模式图
#先定义绘制函数
plot_cellchat_circle <- function(){
  #识别对某些细胞类群的输出和输入信号贡献最大的信号
  netAnalysis_signalingRole_heatmap(stim.cellchat, pattern = "outgoing")
  #netAnalysis_signalingRole_heatmap(stim.cellchat, signaling = c("CXCL", "CCL"),pattern = "outgoing")###亦可指定通路
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_outgoing_signal_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




# 细胞的信号流出通讯模式的鉴定和可视化（细胞通讯的聚类）
selectK(stim.cellchat, pattern = "outgoing")#时长久
nPatterns = 6 # 挑选曲线中第一个出现下降的点(这里都可以试试结果)
stim.cellchat <- identifyCommunicationPatterns(stim.cellchat, pattern = "outgoing", k = nPatterns, 
                                               width = 5, height = 9, font.size = 6)
dev.off()

##查看dotplot每个细胞对通路贡献具体数值

# 1. 读取数值矩阵、标签
out_mat = stim.cellchat@netP$pattern$outgoing$data
cell_vec = levels(stim.cellchat@idents)
path_vec = stim.cellchat@netP$pathways

# 2. 构建完整数据表：行=7种细胞，列=CellType + 33条通路
dot_full_df = as.data.frame(out_mat)
colnames(dot_full_df) = path_vec # 通路赋值为列名
dot_full_df$CellType = cell_vec  # 新增细胞名称列
# 调整列顺序，CellType放第一列
dot_full_df = dot_full_df[, c("CellType", path_vec)]

# 3. 查看全部数据（7行34列，第一列细胞，后面33列通路贡献值）
print(dot_full_df, row.names = FALSE)

# 4. 导出全部细胞所有通路分数CSV
write.csv(dot_full_df, "Dotplot_Outgoing_AllCell_AllPath_Value_patient.csv", row.names = F)



#先定义绘制函数
plot_cellchat_circle <- function(){
##riverplot
netAnalysis_river(stim.cellchat, pattern = "outgoing")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_riverplot_outgoing_communication_patterns_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



#先定义绘制函数
plot_cellchat_circle <- function(){
  ##dotplot
  netAnalysis_dot(stim.cellchat, pattern = "outgoing")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_dotplot_outgoing_communication_patterns_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()






selectK(stim.cellchat, pattern = "incoming")#时长久
nPatterns = 5 # 挑选曲线中第一个出现下降的点(这里都可以试试结果)
stim.cellchat <- identifyCommunicationPatterns(stim.cellchat, pattern = "incoming", k = nPatterns, 
                                               width = 5, height = 9, font.size = 6)
dev.off()

##查看dotplot每个细胞对通路贡献具体数值

# 1. 读取数值矩阵、标签
out_mat = stim.cellchat@netP$pattern$incoming$data
cell_vec = levels(stim.cellchat@idents)
path_vec = stim.cellchat@netP$pathways

# 2. 构建完整数据表：行=7种细胞，列=CellType + 33条通路
dot_full_df = as.data.frame(out_mat)
colnames(dot_full_df) = path_vec # 通路赋值为列名
dot_full_df$CellType = cell_vec  # 新增细胞名称列
# 调整列顺序，CellType放第一列
dot_full_df = dot_full_df[, c("CellType", path_vec)]

# 3. 查看全部数据（7行34列，第一列细胞，后面33列通路贡献值）
print(dot_full_df, row.names = FALSE)

# 4. 导出全部细胞所有通路分数CSV
write.csv(dot_full_df, "Dotplot_Incoming_AllCell_AllPath_Value_patient.csv", row.names = F)



#先定义绘制函数
plot_cellchat_circle <- function(){
  ##riverplot
  netAnalysis_river(stim.cellchat, pattern = "incoming")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_riverplot_incoming_communication_patterns_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



#先定义绘制函数
plot_cellchat_circle <- function(){
  ##dotplot
  netAnalysis_dot(stim.cellchat, pattern = "incoming")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_dotplot_incoming_communication_patterns_patient.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()


}

####13.2、control组cellchat分析####
{
#====读取数据，计算=========================================================
RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")

# V5分层矩阵合并+精简，减少subset内存占用
RNA_harmony_annotation <- JoinLayers(RNA_harmony_annotation, assay = "RNA")
RNA_harmony_annotation <- DietSeurat(RNA_harmony_annotation, assays = "RNA", dimreducs = NULL, graphs = FALSE)

ctrl.object <- subset(RNA_harmony_annotation,group=="control") #提取对照组control
rm(RNA_harmony_annotation) # 释放原始完整Seurat对象，节省内存
gc(verbose=F)

ctrl.data.input <- GetAssayData(ctrl.object, assay = "RNA", layer = "data")
ctrl.meta = ctrl.object@meta.data[,c("celltype", "group")] 
ctrl.meta$CellType %<>% as.vector(.)
rm(ctrl.object)
gc(verbose=F)

ctrl.cellchat <- createCellChat(object = ctrl.data.input)
ctrl.cellchat <- addMeta(ctrl.cellchat, meta = ctrl.meta)
ctrl.cellchat <- setIdent(ctrl.cellchat, ident.use = "celltype")

levels(ctrl.cellchat@idents)#查看下细胞类型
groupSize1 <- as.numeric(table(ctrl.cellchat@idents)) # 各种类型细胞数量
groupSize1


# 人类样本使用人配受体数据库
ctrl.cellchat@DB <- CellChatDB.human 
ctrl.cellchat <- subsetData(ctrl.cellchat) 

ctrl.cellchat <- identifyOverExpressedGenes(ctrl.cellchat)

ctrl.cellchat <- identifyOverExpressedInteractions(ctrl.cellchat)

# 并行计算完成切换串行，释放并行进程内存
plan("sequential")


# 增加trim=10过滤小众细胞群，大幅提速避免长时间卡死；raw.use=T和patient保持统一
ctrl.cellchat <- computeCommunProb(ctrl.cellchat, raw.use = T, trim = 10)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
ctrl.cellchat <- filterCommunication(ctrl.cellchat, min.cells = 10)

ctrl.cellchat <- computeCommunProbPathway(ctrl.cellchat)
ctrl.cellchat <- aggregateNet(ctrl.cellchat)
ctrl.cellchat <- netAnalysis_computeCentrality(ctrl.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#数据查看/保存
group2.net <- subsetCommunication(ctrl.cellchat)  ###细胞通讯结果
write.csv(group2.net, file = "group2_net_inter_raw.useT.csv", row.names = F)

saveRDS(ctrl.cellchat,"ctrl.cellchat.rds")

#====读取数据，绘图=========================================================
ctrl.cellchat <- readRDS("ctrl.cellchat.rds")
library(ggrepel)

groupSize <- as.numeric(table(ctrl.cellchat@idents)) # 各种类型细胞数量
groupSize


#互作网络
#先定义绘制函数
plot_cellchat_circle <- function(){
  # 设置1行2列分屏
  par(mfrow = c(1,2))
  netVisual_circle(ctrl.cellchat@net$count, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   label.edge= F, 
                   title.name = "Number of interactions")
  
  netVisual_circle(ctrl.cellchat@net$weight, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   label.edge= F, 
                   title.name = "Interaction weights/strength")
  #左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。
  #发出箭头的细胞表达配体，箭头指向的细胞表达受体。配体-受体对越多，线越粗。
  #右图：互作的概率/强度值（强度就是概率值相加）
}
#先Rstudio窗口预览
plot_cellchat_circle()
#打开绘图设备，指定保存路径+尺寸进行保存
pdf("CellChat_interaction_network_control.pdf", width = 14, height = 10)
plot_cellchat_circle()#绘制在pdf中
# 关闭绘图设备，文件才会生成
dev.off()




#提取互作矩阵，展示每个亚群作为source的信号传递
mat <- ctrl.cellchat@net$weight#提取每个细胞的互作强度
#先定义绘制函数
plot_cellchat_circle <- function(){
  # 分栏布局 + 放大边距防止文字截断
  par(mfrow = c(3,3), mar = c(0.5,0.5,0.5,0.5))
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    p_base<-netVisual_circle(mat2, 
                             vertex.weight = groupSize,
                             arrow.width = 0.2,arrow.size = 0.1,##自行调整
                             weight.scale = T, 
                             edge.weight.max = max(mat), 
                             title.name = paste0("Source: ", rownames(mat)[i]) 
    )
  }
}
#先Rstudio窗口预览
plot_cellchat_circle()
#打开绘图设备，指定保存路径+尺寸进行保存
pdf("CellChat_Source_All_Subpop_control.pdf", width = 16, height = 14)
plot_cellchat_circle()#绘制在pdf中
# 关闭绘图设备，写入文件（必不可少）
dev.off()




#自定义特定细胞群展示到画板上
specific_cell <- "T cell"
#先定义绘制函数
plot_cellchat_circle <- function(){
  # 分栏布局 + 放大边距防止文字截断
  par(mfrow = c(2,2), mar = c(1,1,1,1))
  cell_order <- rownames(mat)
  # 找到要绘制细胞对应的行号
  idx <- which(cell_order == specific_cell)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[idx, ] <- mat[idx, ] #绘制第几类细胞
  netVisual_circle(mat2, 
                   vertex.weight = groupSize,
                   arrow.width = 0.2,arrow.size = 0.1,##自行调整
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[2])
}
#先Rstudio窗口预览
plot_cellchat_circle()
# 打开保存设备，设置宽高适配3行3张子图
pdf(paste0("CellChat_", specific_cell, "_Source_control.pdf"), width = 16, height = 14)
plot_cellchat_circle()#绘制在pdf中
# 关闭绘图设备，写入文件（必不可少）
dev.off()



#列举所有细胞的全部通路pathway并进行保存为Excel
#层级图
ctrl.cellchat@netP$pathways #全部通路名称列举
group1.net <- subsetCommunication(ctrl.cellchat)  ###细胞通讯全量结果
# 直接导出全量结果到Excel
write.csv(group1.net, file = "control组_所有通路_细胞通讯全量结果.csv", row.names = FALSE)



#可选择感兴趣的通路pathway进行可视化
pathway.show <- "CCL"#这里很关键，指定生物学通路，根据研究选择来解释生物学意义
levels(ctrl.cellchat@idents)

#先定义绘制函数
plot_cellchat_circle <- function(){
  vertex.receiver = c(4,6)#选择想看的接收者细胞类别，数字是指细胞类型的排序号
  netVisual_aggregate(ctrl.cellchat, 
                      signaling = pathway.show,
                      vertex.receiver = vertex.receiver,
                      layout = "hierarchy")
  #在层次图中，实体圆和空心圆分别表示源和目标。
  #线越粗，互作信号越强。
  #左图中间的target是我们选定的靶细胞。
  #右图是选中的靶细胞之外的另外一组放在中间看互作。
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_receiver_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



#circle plot
#先定义绘制函数
plot_cellchat_circle <- function(){
  par(mfrow = c(1,1))
  netVisual_aggregate(ctrl.cellchat, 
                      signaling = pathway.show, 
                      layout = "circle")
  title(main = paste0(pathway.show, " Signaling Pathway Cell Communication Network")) # 外部添加大图标题
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_circle_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




#chord plot
#先定义绘制函数
plot_cellchat_circle <- function(){
  par(mfrow = c(1,1))
  netVisual_aggregate(ctrl.cellchat, 
                      signaling = pathway.show, 
                      layout = "chord")#图中细胞下颜色表示通路下不同的基因互配对
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_chord_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




#热图
#chord plot
#先定义绘制函数
plot_cellchat_circle <- function(){
  netVisual_heatmap(ctrl.cellchat, signaling = pathway.show, color.heatmap = "Reds")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_heatmap_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()





#计算受配体对整个信号通路的贡献，并可视化单个配体-受体对介导的细胞-细胞通讯
netAnalysis_contribution(ctrl.cellchat, signaling = pathway.show)
#可视化由单个配体-受体对介导的细胞间通讯
pairLR.CCL <- extractEnrichedLR(ctrl.cellchat, signaling = pathway.show, geneLR.return = FALSE)
#提取对这个通路显著性最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.CCL[1,] #选择配受体进行绘制
#先定义绘制函数
plot_cellchat_circle <- function(){
  netVisual_individual(ctrl.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "circle")
  netVisual_individual(ctrl.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "chord")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_", pathway.show, "_pathway_",LR.show,"_contribution_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




levels(ctrl.cellchat@idents)
#指定受体细胞和配体细胞
#这里选择特定细胞通用，查看特定细胞的配体和其他细胞的受体
specific_cell_ligand_receptor <- "T cell"
#先定义绘制函数
plot_cellchat_circle <- function(){
  p =netVisual_bubble(ctrl.cellchat, 
                      sources.use = specific_cell_ligand_receptor,#这里改细胞
                      remove.isolate = FALSE,
                      font.size=14)
  p

}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_",specific_cell_ligand_receptor,"_ligand_receptor_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



##参与某条信号通路的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）
#先定义绘制函数
plot_cellchat_circle <- function(){
  plotGeneExpression(ctrl.cellchat, signaling = pathway.show) #选择信号通路
  library(RColorBrewer)
  # 定义颜色渐变
  colors <- brewer.pal(8, "Oranges")
  plotGeneExpression(ctrl.cellchat, signaling = pathway.show, type = "dot", col = colors)
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_",pathway.show,"_pathway_gene_expression_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




#某条通路的发送者（源）和接收者（目标）情况
#计算和可视化网络中心性评分
ctrl.cellchat <- netAnalysis_computeCentrality(ctrl.cellchat, slot.name = "netP")
#先定义绘制函数
plot_cellchat_circle <- function(){
  #主要注意sender和receiver
  netAnalysis_signalingRole_network(ctrl.cellchat, signaling = pathway.show, 
                                    width = 15, height = 6, font.size = 10)
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_",pathway.show,"_pathway_send_receiver_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



#各类细胞所有信号通路的输入输出强度点图
#先定义绘制函数
plot_cellchat_circle <- function(){
  #使用散点图在 2D 空间中可视化主要的发送者（源）和接收者（目标）。
  netAnalysis_signalingRole_scatter(ctrl.cellchat)###ALL
  #netAnalysis_signalingRole_scatter(ctrl.cellchat, signaling = c("CXCL", "CCL"))###亦可指定通路
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_income_outgo_strength_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



#各类细胞所有信号通路的输入信号模式图（（这个地方总部不行提示里面有0））
#先定义绘制函数
plot_cellchat_circle <- function(){
  #识别对某些细胞类群的输出和输入信号贡献最大的信号
  netAnalysis_signalingRole_heatmap(ctrl.cellchat, pattern = "incoming")
  #netAnalysis_signalingRole_heatmap(ctrl.cellchat, signaling = c("CXCL", "CCL"),pattern = "incoming")###亦可指定通路
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_incomeing_signal_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()


#各类细胞所有信号通路的输出信号模式图
#先定义绘制函数
plot_cellchat_circle <- function(){
  #识别对某些细胞类群的输出和输入信号贡献最大的信号
  netAnalysis_signalingRole_heatmap(ctrl.cellchat, pattern = "outgoing")
  #netAnalysis_signalingRole_heatmap(ctrl.cellchat, signaling = c("CXCL", "CCL"),pattern = "outgoing")###亦可指定通路
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_outgoing_signal_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()




# 细胞的信号流出通讯模式的鉴定和可视化（细胞通讯的聚类）
selectK(ctrl.cellchat, pattern = "outgoing")#时长久
nPatterns = 6 # 挑选曲线中第一个出现下降的点(这里都可以试试结果)
ctrl.cellchat <- identifyCommunicationPatterns(ctrl.cellchat, pattern = "outgoing", k = nPatterns, 
                                               width = 5, height = 9, font.size = 6)
dev.off()

##查看dotplot每个细胞对通路贡献具体数值

# 1. 读取数值矩阵、标签
out_mat = ctrl.cellchat@netP$pattern$outgoing$data
cell_vec = levels(ctrl.cellchat@idents)
path_vec = ctrl.cellchat@netP$pathways

# 2. 构建完整数据表：行=7种细胞，列=CellType + 33条通路
dot_full_df = as.data.frame(out_mat)
colnames(dot_full_df) = path_vec # 通路赋值为列名
dot_full_df$CellType = cell_vec  # 新增细胞名称列
# 调整列顺序，CellType放第一列
dot_full_df = dot_full_df[, c("CellType", path_vec)]

# 3. 查看全部数据（7行34列，第一列细胞，后面33列通路贡献值）
print(dot_full_df, row.names = FALSE)

# 4. 导出全部细胞所有通路分数CSV
write.csv(dot_full_df, "Dotplot_Outgoing_AllCell_AllPath_Value_control.csv", row.names = F)



#先定义绘制函数
plot_cellchat_circle <- function(){
  ##riverplot
  netAnalysis_river(ctrl.cellchat, pattern = "outgoing")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_riverplot_outgoing_communication_patterns_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



#先定义绘制函数
plot_cellchat_circle <- function(){
  ##dotplot
  netAnalysis_dot(ctrl.cellchat, pattern = "outgoing")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_dotplot_outgoing_communication_patterns_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



selectK(ctrl.cellchat, pattern = "incoming")#时长久
nPatterns = 5 # 挑选曲线中第一个出现下降的点(这里都可以试试结果)
ctrl.cellchat <- identifyCommunicationPatterns(ctrl.cellchat, pattern = "incoming", k = nPatterns, 
                                               width = 5, height = 9, font.size = 6)
dev.off()

##查看dotplot每个细胞对通路贡献具体数值

# 1. 读取数值矩阵、标签
out_mat = ctrl.cellchat@netP$pattern$incoming$data
cell_vec = levels(ctrl.cellchat@idents)
path_vec = ctrl.cellchat@netP$pathways

# 2. 构建完整数据表：行=7种细胞，列=CellType + 33条通路
dot_full_df = as.data.frame(out_mat)
colnames(dot_full_df) = path_vec # 通路赋值为列名
dot_full_df$CellType = cell_vec  # 新增细胞名称列
# 调整列顺序，CellType放第一列
dot_full_df = dot_full_df[, c("CellType", path_vec)]

# 3. 查看全部数据（7行34列，第一列细胞，后面33列通路贡献值）
print(dot_full_df, row.names = FALSE)

# 4. 导出全部细胞所有通路分数CSV
write.csv(dot_full_df, "Dotplot_Incoming_AllCell_AllPath_Value_control.csv", row.names = F)



#先定义绘制函数
plot_cellchat_circle <- function(){
  ##riverplot
  netAnalysis_river(ctrl.cellchat, pattern = "incoming")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_riverplot_incoming_communication_patterns_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()



#先定义绘制函数
plot_cellchat_circle <- function(){
  ##dotplot
  netAnalysis_dot(ctrl.cellchat, pattern = "incoming")
}
#先Rstudio窗口预览
plot_cellchat_circle()
#先开启保存设备（放在所有绘图代码最前面）
pdf(paste0("CellChat_allcell_dotplot_incoming_communication_patterns_control.pdf"), width = 10, height = 7)
plot_cellchat_circle()#绘制在pdf中
#必须关闭设备，文件才会生成
dev.off()


}

####13.3、carrier组cellchat分析####
{  
    #====读取数据，计算=========================================================
    RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
    
    # V5 Seurat 极致精简，删除冗余图层、降维、图网络，大幅减小对象体积
    RNA_harmony_annotation <- JoinLayers(RNA_harmony_annotation, assay = "RNA")
    RNA_harmony_annotation <- DietSeurat(
      RNA_harmony_annotation, 
      assays = "RNA", 
      dimreducs = NULL, 
      graphs = FALSE,
      reductions = NULL,
      features = NULL
    )
    
    # 提取carrier分组
    carrier.object <- subset(RNA_harmony_annotation, group == "carrier")
    rm(RNA_harmony_annotation) # 删除原始大Seurat对象
    gc(verbose = FALSE, reset = TRUE) # 强制垃圾回收释放内存
    
    # 提取标准化表达矩阵 + 精简meta信息
    carrier.data.input <- GetAssayData(carrier.object, assay = "RNA", layer = "data")
    carrier.meta <- carrier.object@meta.data[, c("celltype", "group")]
    carrier.meta$celltype %<>% as.vector()
    rm(carrier.object)
    gc(verbose = FALSE, reset = TRUE)
    
    ###1、构建cellchat对象###
    carrier.cellchat <- createCellChat(object = carrier.data.input)
    carrier.cellchat <- addMeta(carrier.cellchat, meta = carrier.meta)
    carrier.cellchat <- setIdent(carrier.cellchat, ident.use = "celltype")
    
    levels(carrier.cellchat@idents)
    groupSize <- as.numeric(table(carrier.cellchat@idents))
    groupSize
    
    # 人源配受体数据库
    carrier.cellchat@DB <- CellChatDB.human
    dplyr::glimpse(CellChatDB.human$interaction)
    
    # === 关键优化1：只保留配受体基因，压缩矩阵，彻底减小稠密转换开销 ===
    lr_db <- CellChatDB.human$interaction
    lr_gene_list <- unique(c(lr_db$ligand, lr_db$receptor))
    overlap_gene <- intersect(rownames(carrier.cellchat@data), lr_gene_list)
    # 仅保留配受体相关基因，无关基因全部剔除
    carrier.cellchat@data <- carrier.cellchat@data[overlap_gene, ]
    gc(verbose = FALSE, reset = TRUE)
    
    carrier.cellchat <- subsetData(carrier.cellchat, features = NULL)
    
    # === 并行设置：64G机器 workers=2 平衡速度与内存峰值 ===
    future::plan("multisession", workers = 2)
    
    # === 关键优化2：调整thresh/min.cells，减少高表达基因计算量，弱化sparse转稠密内存占用 ===

    suppressWarnings({
      carrier.cellchat <- identifyOverExpressedGenes(carrier.cellchat)
    })
    
    # 释放并行子进程内存
    gc(verbose = FALSE, reset = TRUE)
    
    carrier.cellchat <- identifyOverExpressedInteractions(carrier.cellchat)
    
    # 并行任务结束，切回串行模式，关闭多进程释放内存
    plan("sequential")
    gc(verbose = FALSE, reset = TRUE)
    
    
    ###2、cellchat分析###
    
   
    #通过计算与每个信号通路相关的所有配体-受体相互作用的通信概率来推断信号通路水平上的通信概率。
    # 新增trim过滤小众细胞群，大幅提速，解决长时间卡死
    carrier.cellchat <- computeCommunProb(carrier.cellchat,raw.use=T, trim = 10)
   
    # 过滤掉小于10个细胞的胞间通讯网络，通讯中的细胞很少没有意义
    carrier.cellchat <- filterCommunication(carrier.cellchat, min.cells = 10)
    # 通过汇总所有相关的配体/受体，计算信号通路水平上的通信概率
    carrier.cellchat <- computeCommunProbPathway(carrier.cellchat)
    carrier.cellchat <- aggregateNet(carrier.cellchat)#计算聚合网络
    # “netP”表示推断的信号通路的细胞间通信网络
    carrier.cellchat <- netAnalysis_computeCentrality(carrier.cellchat, slot.name = "netP") 
    
    #数据查看/保存
    group1.net <- subsetCommunication(carrier.cellchat)  ###细胞通讯结果
    write.csv(group1.net, file = "group1_net_inter_raw.useT_carrier.csv", row.names = F)
    saveRDS(carrier.cellchat,"carrier.cellchat.rds")

    
    
    ###3、cellchat结果可视化###
    #====读取数据，绘图=========================================================
    carrier.cellchat <- readRDS("carrier.cellchat.rds")
    library(ggrepel)
    
    groupSize <- as.numeric(table(carrier.cellchat@idents)) # 各种类型细胞数量
    groupSize
    
    #互作网络
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      # 设置1行2列分屏
      par(mfrow = c(1,2))
      netVisual_circle(carrier.cellchat@net$count, 
                       vertex.weight = groupSize, 
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = "Number of interactions")
      
      netVisual_circle(carrier.cellchat@net$weight, 
                       vertex.weight = groupSize, 
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = "Interaction weights/strength")
      #左图：外周各种颜色圆圈的大小表示细胞的数量，圈越大，细胞数越多。
      #发出箭头的细胞表达配体，箭头指向的细胞表达受体。配体-受体对越多，线越粗。
      #右图：互作的概率/强度值（强度就是概率值相加）
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_interaction_network_carrier.pdf", width = 14, height = 10)
    plot_cellchat_circle()#绘制在pdf中
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    
    
    #提取互作矩阵，展示每个亚群作为source的信号传递
    mat <- carrier.cellchat@net$weight#提取每个细胞的互作强度
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      # 分栏布局 + 放大边距防止文字截断
      par(mfrow = c(3,3), mar = c(0.5,0.5,0.5,0.5))
      for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        p_base<-netVisual_circle(mat2, 
                                 vertex.weight = groupSize,
                                 arrow.width = 0.2,arrow.size = 0.1,##自行调整
                                 weight.scale = T, 
                                 edge.weight.max = max(mat), 
                                 title.name = paste0("Source: ", rownames(mat)[i]) 
        )
      }
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_Source_All_Subpop_carrier.pdf", width = 16, height = 14)
    plot_cellchat_circle()#绘制在pdf中
    # 关闭绘图设备，写入文件（必不可少）
    dev.off()
    
    
    
    
    #自定义特定细胞群展示到画板上
    specific_cell <- "NK cell"
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      # 分栏布局 + 放大边距防止文字截断
      par(mfrow = c(2,2), mar = c(1,1,1,1))
      cell_order <- rownames(mat)
      # 找到要绘制细胞对应的行号
      idx <- which(cell_order == specific_cell)
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[idx, ] <- mat[idx, ] #绘制第几类细胞
      netVisual_circle(mat2, 
                       vertex.weight = groupSize,
                       arrow.width = 0.2,arrow.size = 0.1,##自行调整
                       weight.scale = T, 
                       edge.weight.max = max(mat), 
                       title.name = specific_cell)
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    # 打开保存设备，设置宽高适配3行3张子图
    pdf(paste0("CellChat_", specific_cell, "_Source_carrier.pdf"), width = 16, height = 14)
    plot_cellchat_circle()#绘制在pdf中
    # 关闭绘图设备，写入文件（必不可少）
    dev.off()
    
    
    
    #列举所有细胞的全部通路pathway并进行保存为Excel
    #层级图
    carrier.cellchat@netP$pathways #全部通路名称列举
    group1.net <- subsetCommunication(carrier.cellchat)  ###细胞通讯全量结果
    # 直接导出全量结果到Excel
    write.csv(group1.net, file = "carrier组_所有通路_细胞通讯全量结果.csv", row.names = FALSE)
    
    
    
    #可选择感兴趣的通路pathway进行可视化
    pathway.show <- "CCL"#这里很关键，指定生物学通路，根据研究选择来解释生物学意义
    levels(carrier.cellchat@idents)
    
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      vertex.receiver = c(4,6)#【重要！运行后根据levels(carrier.cellchat@idents)修改数字】
      netVisual_aggregate(carrier.cellchat, 
                          signaling = pathway.show,
                          vertex.receiver = vertex.receiver,
                          layout = "hierarchy")
      #在层次图中，实体圆和空心圆分别表示源和目标。
      #线越粗，互作信号越强。
      #左图中间的target是我们选定的靶细胞。
      #右图是选中的靶细胞之外的另外一组放在中间看互作。
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf(paste0("CellChat_", pathway.show, "_pathway_receiver_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    #circle plot
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      par(mfrow = c(1,1))
      netVisual_aggregate(carrier.cellchat, 
                          signaling = pathway.show, 
                          layout = "circle")
      title(main = paste0(pathway.show, " Signaling Pathway Cell Communication Network")) # 外部添加大图标题
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf(paste0("CellChat_", pathway.show, "_pathway_circle_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    
    #chord plot
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      par(mfrow = c(1,1))
      netVisual_aggregate(carrier.cellchat, 
                          signaling = pathway.show, 
                          layout = "chord")#图中细胞下颜色表示通路下不同的基因互配对
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf(paste0("CellChat_", pathway.show, "_pathway_chord_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    
    #热图
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      netVisual_heatmap(carrier.cellchat, signaling = pathway.show, color.heatmap = "Reds")
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf(paste0("CellChat_", pathway.show, "_pathway_heatmap_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    #计算受配体对整个信号通路的贡献，并可视化单个配体-受体对介导的细胞-细胞通讯
    netAnalysis_contribution(carrier.cellchat, signaling = pathway.show)
    #可视化由单个配体-受体对介导的细胞间通讯
    pairLR.CCL <- extractEnrichedLR(carrier.cellchat, signaling = pathway.show, geneLR.return = FALSE)
    #提取对这个通路显著性最大的配体受体对来展示（也可以选择其他的配体受体对）
    LR.show <- pairLR.CCL[1,] #选择配受体进行绘制
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      netVisual_individual(carrier.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "circle")
      netVisual_individual(carrier.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "chord")
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf(paste0("CellChat_", pathway.show, "_pathway_",LR.show,"_contribution_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    
    levels(carrier.cellchat@idents)
    #指定受体细胞和配体细胞
    #这里选择特定细胞通用，查看特定细胞的配体和其他细胞的受体
    specific_cell_ligand_receptor <- "NK cell"
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      p =netVisual_bubble(carrier.cellchat, 
                          sources.use = specific_cell_ligand_receptor,#这里改细胞
                          remove.isolate = FALSE,
                          font.size=14)
      p
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf(paste0("CellChat_",specific_cell_ligand_receptor,"_ligand_receptor_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    
    
    
    ##参与某条信号通路的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      plotGeneExpression(carrier.cellchat, signaling = pathway.show) #选择信号通路
      library(RColorBrewer)
      # 定义颜色渐变
      colors <- brewer.pal(8, "Oranges")
      plotGeneExpression(carrier.cellchat, signaling = pathway.show, type = "dot", col = colors)
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf(paste0("CellChat_",pathway.show,"_pathway_gene_expression_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    
    #某条通路的发送者（源）和接收者（目标）情况
    #计算和可视化网络中心性评分
    carrier.cellchat <- netAnalysis_computeCentrality(carrier.cellchat, slot.name = "netP")
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      #主要注意sender和receiver
      netAnalysis_signalingRole_network(carrier.cellchat, signaling = pathway.show, 
                                        width = 15, height = 6, font.size = 10)
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf(paste0("CellChat_",pathway.show,"_pathway_send_receiver_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    
    #各类细胞所有信号通路的输入输出强度点图
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      #使用散点图在 2D 空间中可视化主要的发送者（源）和接收者（目标）。
      netAnalysis_signalingRole_scatter(carrier.cellchat)###ALL
      #netAnalysis_signalingRole_scatter(carrier.cellchat, signaling = c("CXCL", "CCL"))###亦可指定通路
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf("CellChat_allcell_income_outgo_strength_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    #各类细胞所有信号通路的输入信号模式图
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      #识别对某些细胞类群的输出和输入信号贡献最大的信号
      netAnalysis_signalingRole_heatmap(carrier.cellchat, pattern = "incoming")
      #netAnalysis_signalingRole_heatmap(carrier.cellchat, signaling = c("CXCL", "CCL"),pattern = "incoming")###亦可指定通路
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf("CellChat_allcell_incomeing_signal_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    
    #各类细胞所有信号通路的输出信号模式图
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      #识别对某些细胞类群的输出和输入信号贡献最大的信号
      netAnalysis_signalingRole_heatmap(carrier.cellchat, pattern = "outgoing")
      #netAnalysis_signalingRole_heatmap(carrier.cellchat, signaling = c("CXCL", "CCL"),pattern = "outgoing")###亦可指定通路
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf("CellChat_allcell_outgoing_signal_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    
    # 细胞的信号流出通讯模式的鉴定和可视化（细胞通讯的聚类）
    selectK(carrier.cellchat, pattern = "outgoing")#时长久
    nPatterns = 6 # 挑选曲线中第一个出现下降的点(这里都可以试试结果)
    carrier.cellchat <- identifyCommunicationPatterns(carrier.cellchat, pattern = "outgoing", k = nPatterns, 
                                                      width = 5, height = 9, font.size = 6)
    dev.off()
    
    ##查看dotplot每个细胞对通路贡献具体数值
    
    # 1. 读取数值矩阵、标签
    out_mat = carrier.cellchat@netP$pattern$outgoing$data
    cell_vec = levels(carrier.cellchat@idents)
    path_vec = carrier.cellchat@netP$pathways
    
    # 2. 构建完整数据表：行=7种细胞，列=CellType + 33条通路
    dot_full_df = as.data.frame(out_mat)
    colnames(dot_full_df) = path_vec # 通路赋值为列名
    dot_full_df$CellType = cell_vec  # 新增细胞名称列
    # 调整列顺序，CellType放第一列
    dot_full_df = dot_full_df[, c("CellType", path_vec)]
    
    # 3. 查看全部数据
    print(dot_full_df, row.names = FALSE)
    
    # 4. 导出全部细胞所有通路分数CSV
    write.csv(dot_full_df, "Dotplot_Outgoing_AllCell_AllPath_Value_carrier.csv", row.names = F)
    
    
    
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      ##riverplot
      netAnalysis_river(carrier.cellchat, pattern = "outgoing")
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf("CellChat_allcell_riverplot_outgoing_communication_patterns_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      ##dotplot
      netAnalysis_dot(carrier.cellchat, pattern = "outgoing")
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf("CellChat_allcell_dotplot_outgoing_communication_patterns_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    selectK(carrier.cellchat, pattern = "incoming")#时长久
    nPatterns = 4 # 挑选曲线中第一个出现下降的点(这里都可以试试结果)
    carrier.cellchat <- identifyCommunicationPatterns(carrier.cellchat, pattern = "incoming", k = nPatterns, 
                                                      width = 5, height = 9, font.size = 6)
    dev.off()
    
    ##查看dotplot每个细胞对通路贡献具体数值
    
    # 1. 读取数值矩阵、标签
    out_mat = carrier.cellchat@netP$pattern$incoming$data
    cell_vec = levels(carrier.cellchat@idents)
    path_vec = carrier.cellchat@netP$pathways
    
    # 2. 构建完整数据表
    dot_full_df = as.data.frame(out_mat)
    colnames(dot_full_df) = path_vec # 通路赋值为列名
    dot_full_df$CellType = cell_vec  # 新增细胞名称列
    # 调整列顺序，CellType放第一列
    dot_full_df = dot_full_df[, c("CellType", path_vec)]
    
    # 3. 查看全部数据
    print(dot_full_df, row.names = FALSE)
    
    # 4. 导出全部细胞所有通路分数CSV
    write.csv(dot_full_df, "Dotplot_Incoming_AllCell_AllPath_Value_carrier.csv", row.names = F)
    
    
    
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      ##riverplot
      netAnalysis_river(carrier.cellchat, pattern = "incoming")
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf("CellChat_allcell_riverplot_incoming_communication_patterns_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
    
    #先定义绘制函数
    plot_cellchat_circle <- function(){
      ##dotplot
      netAnalysis_dot(carrier.cellchat, pattern = "incoming")
    }
    #先Rstudio窗口预览
    plot_cellchat_circle()
    #先开启保存设备（放在所有绘图代码最前面）
    pdf("CellChat_allcell_dotplot_incoming_communication_patterns_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()#绘制在pdf中
    #必须关闭设备，文件才会生成
    dev.off()
    
    
  }
  

####13.4、patient和control合并后比较cellchat结果可视化####
  {
    #====读取数据，计算=========================================================
    stim.cellchat <- readRDS("stim.cellchat.rds")
    ctrl.cellchat <- readRDS("ctrl.cellchat.rds")
    
    ###合并两组结果###
    object.list <- list(control = ctrl.cellchat, patient = stim.cellchat)
    
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    
    saveRDS(cellchat,"cellchat_patient_control.rds")
    
    #====读取数据，绘图=========================================================
    cellchat <- readRDS("cellchat_patient_control.rds")
    ## 细胞间互作次数bar图
    #比较两组互作数目
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    gg_all <-gg1 + gg2
    gg_all
    # 直接ggsave保存，无需pdf/dev.off
    ggsave("Compare_Interactions_bar_patient&control.pdf", gg_all, width = 12, height = 5, dpi = 300)
    dev.off()
    
    
    library(writexl)  # 用于导出xlsx文件（与你的write_xlsx对应）
    # 明确分组名称（固定为control和patient）
    group1 <- "control"  # 对照组
    group2 <- "patient"  # 患者组
    # 提取互作数目（count）和权重（weight）矩阵
    control_count_mat <- cellchat@net[[group1]]$count
    patient_count_mat <- cellchat@net[[group2]]$count
    control_weight_mat <- cellchat@net[[group1]]$weight
    patient_weight_mat <- cellchat@net[[group2]]$weight
    # 验证矩阵有效性
    cat("control组互作数目矩阵维度：", dim(control_count_mat), "\n")
    cat("patient组互作数目矩阵维度：", dim(patient_count_mat), "\n")
    cat("control组互作权重矩阵维度：", dim(control_weight_mat), "\n")
    cat("patient组互作权重矩阵维度：", dim(patient_weight_mat), "\n")
    
    # 通用数据整理函数
    tidy_interaction_mat <- function(mat1, mat2, group1_name, group2_name, measure_type) {
      # 处理组1矩阵
      mat1_tidy <- as.data.frame(mat1) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group1_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      # 处理组2矩阵
      mat2_tidy <- as.data.frame(mat2) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group2_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      # 合并两组数据
      tidy_df <- merge(mat1_tidy, mat2_tidy, by = c("source", "target"), all = TRUE)
      tidy_df[is.na(tidy_df)] <- 0
      return(tidy_df)
    }
    
    # ---- 处理互作数目（count）数据，计算统计指标 ---
    # 1. 整理count数据为长格式
    count_tidy_df <- tidy_interaction_mat(
      mat1 = control_count_mat,
      mat2 = patient_count_mat,
      group1_name = "control",
      group2_name = "patient",
      measure_type = "count"
    )
    
    # 2. 逐个源-靶对计算统计指标（修正select函数，处理wilcox警告）
    count_stats_df <- count_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        control_count = count_control,
        patient_count = count_patient,
        # 计算logFC（加1e-6避免除0）
        logFC = log2((patient_count + 1e-6) / (control_count + 1e-6)),
        # 计算原始p值：添加exact=FALSE减少连结警告，不影响结果
        pvalue = wilcox.test(
          x = c(control_count, rep(0, 9)),
          y = c(patient_count, rep(0, 9)),
          exact = FALSE  # 关键：关闭精确检验，减少警告
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      # FDR校正
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      # 明确使用dplyr::select，避免包冲突
      dplyr::select(source, target, pvalue, p.adjust, logFC, control_count, patient_count)
    
    # 3. 重命名为中文列名
    p_count_table <- count_stats_df %>%
      dplyr::rename(
        源细胞类型 = source,
        靶细胞类型 = target,
        原始p值_count = pvalue,
        FDR校正p值_count = p.adjust,
        差异倍数_count = logFC,
        control组互作数目 = control_count,
        patient组互作数目 = patient_count
      )
    
    # ---- 处理互作权重（weight）数据，计算统计指标 ---
    # 1. 整理weight数据为长格式
    weight_tidy_df <- tidy_interaction_mat(
      mat1 = control_weight_mat,
      mat2 = patient_weight_mat,
      group1_name = "control",
      group2_name = "patient",
      measure_type = "weight"
    )
    
    # 2. 逐个源-靶对计算统计指标
    weight_stats_df <- weight_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        control_weight = weight_control,
        patient_weight = weight_patient,
        logFC = log2((patient_weight + 1e-6) / (control_weight + 1e-6)),
        # 关闭精确检验，减少警告
        pvalue = wilcox.test(
          x = c(control_weight, rep(0, 9)),
          y = c(patient_weight, rep(0, 9)),
          exact = FALSE
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      # 明确使用dplyr::select
      dplyr::select(source, target, pvalue, p.adjust, logFC, control_weight, patient_weight)
    
    # 3. 重命名为中文列名
    p_weight_table <- weight_stats_df %>%
      dplyr::rename(
        源细胞类型 = source,
        靶细胞类型 = target,
        原始p值_weight = pvalue,
        FDR校正p值_weight = p.adjust,
        差异倍数_weight = logFC,
        control组互作权重 = control_weight,
        patient组互作权重 = patient_weight
      )
    
    # ---- 合并表格+导出Excel ---
    combine_p_table <- merge(
      p_count_table,
      p_weight_table,
      by = c("源细胞类型", "靶细胞类型"),
      all = TRUE
    )
    
    # 导出Excel
    write_xlsx(combine_p_table, "gg1_gg2_互作网络p值综合表（patient_vs_control）.xlsx") #有问题20260906，这里p值做的有问题，应该用每个样本人分别跑，不能混在一起或者混在同一个组
    cat("表格导出成功！\n")
    
    
    
    ## 细胞间互作次数网络图
    # 红色为STIM组相比CTRL组互作次数和互作强度增加，线越粗表示差异越大；蓝色则表示减少。
    # 两个数据集之间的细胞-细胞通信网络中的交互或交互强度的差异数量可以使用圆形图来可视化，其中红色边表示第二个数据集相比于第一个数据集增加的信号。
    
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_compare_interaction_network_patient&control.pdf", width = 14, height = 10)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_compare_interaction_heatmap_patient&control.pdf", width = 14, height = 10)
    #数量与强度差异热图(主要看互作强度，看质不看量)
    par(mfrow = c(1,1))
    h1 <- netVisual_heatmap(cellchat)
    h2 <- netVisual_heatmap(cellchat, measure = "weight")
    h1+h2
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    
    #保守和特异性信号通路的识别与可视化
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
    gg_all <- gg1 + gg2
    gg_all
    # 直接保存
    ggsave("CellChat_rankNet_comparison_patient&control.pdf", gg_all, width = 14, height = 6, dpi = 300)
    dev.off()
    
    
    
    diff.count <- cellchat@net$patient$count - cellchat@net$control$count
    write.csv(cellchat@net$patient$count, "output_STIMcount_patient&control.csv", quote = F)
    write.csv(cellchat@net$control$count, "output_CTRLcount_patient&control.csv", quote = F)
    
    library(pheatmap)
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_diff_count_heatmap_patient&control.pdf", width = 14, height = 10)
    pheatmap(diff.count,
             treeheight_row = "0",treeheight_col = "0",#不画树
             cluster_rows=T, 
             cluster_cols=T)
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    # cellchat里有丰富的“受体-配体”分析结果。结合自己的假设，设计合适的图形展示数据。
    View(cellchat)
    View(cellchat@LR[["patient"]][["LRsig"]])
    
  }
  
  
####13.5、carrier和control合并后cellchat结果可视化####
  {
    #===加载必需包===
    library(CellChat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(pheatmap)
    
    #====读取数据，计算=========================================================
    carrier.cellchat <- readRDS("carrier.cellchat.rds")
    ctrl.cellchat <- readRDS("ctrl.cellchat.rds")
    
    ###合并两组结果###
    object.list <- list(control = ctrl.cellchat, carrier = carrier.cellchat)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    saveRDS(cellchat,"cellchat_carrier_control.rds")
    
    
    
    #====读取数据，绘图=========================================================
    cellchat <- readRDS("cellchat_carrier_control.rds")
    
    ## 1. 细胞间互作总数对比柱状图
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    gg_all <- gg1 + gg2
    gg_all
    # 直接ggsave保存，无需pdf/dev.off
    ggsave("Compare_Interactions_bar_carrier&control.pdf", gg_all, width = 12, height = 5, dpi = 300)
    dev.off()
    
    # 安装/加载导出Excel包
    if (!require(writexl)) {
      install.packages("writexl", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
      library(writexl)
    }
    
    # 定义分组
    group1 <- "control"
    group2 <- "carrier"
    
    # 提取count与weight矩阵
    control_count_mat <- cellchat@net[[group1]]$count
    carrier_count_mat <- cellchat@net[[group2]]$count
    control_weight_mat <- cellchat@net[[group1]]$weight
    carrier_weight_mat <- cellchat@net[[group2]]$weight
    
    # 打印矩阵维度检查
    cat("control组互作数目矩阵维度：", dim(control_count_mat), "\n")
    cat("carrier组互作数目矩阵维度：", dim(carrier_count_mat), "\n")
    cat("control组互作权重矩阵维度：", dim(control_weight_mat), "\n")
    cat("carrier组互作权重矩阵维度：", dim(carrier_weight_mat), "\n")
    
    # 矩阵转长数据通用函数
    tidy_interaction_mat <- function(mat1, mat2, group1_name, group2_name, measure_type) {
      mat1_tidy <- as.data.frame(mat1) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group1_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      mat2_tidy <- as.data.frame(mat2) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group2_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      tidy_df <- merge(mat1_tidy, mat2_tidy, by = c("source", "target"), all = TRUE)
      tidy_df[is.na(tidy_df)] <- 0
      return(tidy_df)
    }
    
    # ---- count互作数量统计 ---
    count_tidy_df <- tidy_interaction_mat(
      mat1 = control_count_mat,
      mat2 = carrier_count_mat,
      group1_name = "control",
      group2_name = "carrier",
      measure_type = "count"
    )
    
    count_stats_df <- count_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        control_count = count_control,
        carrier_count = count_carrier,
        logFC = log2((carrier_count + 1e-6) / (control_count + 1e-6)),
        pvalue = wilcox.test(
          x = c(control_count, rep(0, 9)),
          y = c(carrier_count, rep(0, 9)),
          exact = FALSE
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      dplyr::select(source, target, pvalue, p.adjust, logFC, control_count, carrier_count)
    
    p_count_table <- count_stats_df %>%
      dplyr::rename(
        源细胞类型 = source,
        靶细胞类型 = target,
        原始p值_count = pvalue,
        FDR校正p值_count = p.adjust,
        差异倍数_count = logFC,
        control组互作数目 = control_count,
        carrier组互作数目 = carrier_count
      )
    
    # ---- weight互作强度统计 ---
    weight_tidy_df <- tidy_interaction_mat(
      mat1 = control_weight_mat,
      mat2 = carrier_weight_mat,
      group1_name = "control",
      group2_name = "carrier",
      measure_type = "weight"
    )
    
    weight_stats_df <- weight_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        control_weight = weight_control,
        carrier_weight = weight_carrier,
        logFC = log2((carrier_weight + 1e-6) / (control_weight + 1e-6)),
        pvalue = wilcox.test(
          x = c(control_weight, rep(0, 9)),
          y = c(carrier_weight, rep(0, 9)),
          exact = FALSE
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      dplyr::select(source, target, pvalue, p.adjust, logFC, control_weight, carrier_weight)
    
    p_weight_table <- weight_stats_df %>%
      dplyr::rename(
        源细胞类型 = source,
        靶细胞类型 = target,
        原始p值_weight = pvalue,
        FDR校正p值_weight = p.adjust,
        差异倍数_weight = logFC,
        control组互作权重 = control_weight,
        carrier组互作权重 = carrier_weight
      )
    
    # 合并count+weight结果并导出Excel
    combine_p_table <- merge(
      p_count_table,
      p_weight_table,
      by = c("源细胞类型", "靶细胞类型"),
      all = TRUE
    )
    write_xlsx(combine_p_table, "互作网络p值综合表（carrier_vs_control）.xlsx")  #有问题20260906，这里p值做的有问题，应该用每个样本人分别跑，不能混在一起或者混在同一个组
    cat("Excel表格导出成功！\n")
    
    
    ## 2. 差异互作圆形网络图
    # 红色=carrier组上调，蓝色=carrier组下调
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_compare_interaction_network_carrier&control.pdf", width = 14, height = 10)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    
    
    ## 3. 互作差异热图
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_compare_interaction_heatmap_carrier&control.pdf", width = 14, height = 10)
    par(mfrow = c(1,1))
    h1 <- netVisual_heatmap(cellchat)
    h2 <- netVisual_heatmap(cellchat, measure = "weight")
    h1 + h2
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    
    
    ## 4. 通路富集对比柱状图
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
    gg_all <- gg1 + gg2
    gg_all
    # 直接保存
    ggsave("CellChat_rankNet_comparison_carrier&control.pdf", gg_all, width = 14, height = 6, dpi = 300)
    dev.off()
    
    
    ## 5. 两组互作差值热图 & 原始矩阵导出
    diff.count <- cellchat@net$carrier$count - cellchat@net$control$count
    write.csv(cellchat@net$carrier$count, "output_carrier_countc_carrier&control.csv", quote = F)
    write.csv(cellchat@net$control$count, "output_CTRLcount_carrier&control.csv", quote = F)
    
    
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_diff_count_heatmap_carrier&control.pdf", width = 14, height = 10)
    pheatmap(diff.count,
             treeheight_row = "0",
             treeheight_col = "0",
             cluster_rows = T,
             cluster_cols = T)
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    # 查看配体受体信号
    View(cellchat)
    View(cellchat@LR[["carrier"]][["LRsig"]])
  }
  
####13.6、patient和carrier合并后比较cellchat结果可视化####
  {
    #====读取数据，计算=========================================================
    stim.cellchat <- readRDS("stim.cellchat.rds")
    carrier.cellchat <- readRDS("carrier.cellchat.rds")
    
    ###合并两组结果###
    object.list <- list(carrier = carrier.cellchat, patient = stim.cellchat)
    
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    
    saveRDS(cellchat,"cellchat_patient_carrier.rds")
    
    #====读取数据，绘图=========================================================
    cellchat <- readRDS("cellchat_patient_carrier.rds")
    ## 细胞间互作次数bar图
    #比较两组互作数目
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    gg_all <-gg1 + gg2
    gg_all
    # 直接ggsave保存，无需pdf/dev.off
    ggsave("Compare_Interactions_bar_patient&carrier.pdf", gg_all, width = 12, height = 5, dpi = 300)
    dev.off()
    
    
    library(writexl)  # 用于导出xlsx文件（与你的write_xlsx对应）
    # 明确分组名称（固定为carrier和patient）
    group1 <- "carrier"  # 携带者组
    group2 <- "patient"  # 患者组
    # 提取互作数目（count）和权重（weight）矩阵
    carrier_count_mat <- cellchat@net[[group1]]$count
    patient_count_mat <- cellchat@net[[group2]]$count
    carrier_weight_mat <- cellchat@net[[group1]]$weight
    patient_weight_mat <- cellchat@net[[group2]]$weight
    # 验证矩阵有效性
    cat("carrier组互作数目矩阵维度：", dim(carrier_count_mat), "\n")
    cat("patient组互作数目矩阵维度：", dim(patient_count_mat), "\n")
    cat("carrier组互作权重矩阵维度：", dim(carrier_weight_mat), "\n")
    cat("patient组互作权重矩阵维度：", dim(patient_weight_mat), "\n")
    
    # 通用数据整理函数
    tidy_interaction_mat <- function(mat1, mat2, group1_name, group2_name, measure_type) {
      # 处理组1矩阵
      mat1_tidy <- as.data.frame(mat1) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group1_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      # 处理组2矩阵
      mat2_tidy <- as.data.frame(mat2) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group2_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      # 合并两组数据
      tidy_df <- merge(mat1_tidy, mat2_tidy, by = c("source", "target"), all = TRUE)
      tidy_df[is.na(tidy_df)] <- 0
      return(tidy_df)
    }
    
    # ---- 处理互作数目（count）数据，计算统计指标 ---
    # 1. 整理count数据为长格式
    count_tidy_df <- tidy_interaction_mat(
      mat1 = carrier_count_mat,
      mat2 = patient_count_mat,
      group1_name = "carrier",
      group2_name = "patient",
      measure_type = "count"
    )
    
    # 2. 逐个源-靶对计算统计指标（修正select函数，处理wilcox警告）
    count_stats_df <- count_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        carrier_count = count_carrier,
        patient_count = count_patient,
        # 计算logFC（加1e-6避免除0）
        logFC = log2((patient_count + 1e-6) / (carrier_count + 1e-6)),
        # 计算原始p值：添加exact=FALSE减少连结警告，不影响结果
        pvalue = wilcox.test(
          x = c(carrier_count, rep(0, 9)),
          y = c(patient_count, rep(0, 9)),
          exact = FALSE  # 关键：关闭精确检验，减少警告
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      # FDR校正
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      # 明确使用dplyr::select，避免包冲突
      dplyr::select(source, target, pvalue, p.adjust, logFC, carrier_count, patient_count)
    
    # 3. 重命名为中文列名
    p_count_table <- count_stats_df %>%
      dplyr::rename(
        源细胞类型 = source,
        靶细胞类型 = target,
        原始p值_count = pvalue,
        FDR校正p值_count = p.adjust,
        差异倍数_count = logFC,
        carrier组互作数目 = carrier_count,
        patient组互作数目 = patient_count
      )
    
    # ---- 处理互作权重（weight）数据，计算统计指标 ---
    # 1. 整理weight数据为长格式
    weight_tidy_df <- tidy_interaction_mat(
      mat1 = carrier_weight_mat,
      mat2 = patient_weight_mat,
      group1_name = "carrier",
      group2_name = "patient",
      measure_type = "weight"
    )
    
    # 2. 逐个源-靶对计算统计指标
    weight_stats_df <- weight_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        carrier_weight = weight_carrier,
        patient_weight = weight_patient,
        logFC = log2((patient_weight + 1e-6) / (carrier_weight + 1e-6)),
        # 关闭精确检验，减少警告
        pvalue = wilcox.test(
          x = c(carrier_weight, rep(0, 9)),
          y = c(patient_weight, rep(0, 9)),
          exact = FALSE
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      # 明确使用dplyr::select
      dplyr::select(source, target, pvalue, p.adjust, logFC, carrier_weight, patient_weight)
    
    # 3. 重命名为中文列名
    p_weight_table <- weight_stats_df %>%
      dplyr::rename(
        源细胞类型 = source,
        靶细胞类型 = target,
        原始p值_weight = pvalue,
        FDR校正p值_weight = p.adjust,
        差异倍数_weight = logFC,
        carrier组互作权重 = carrier_weight,
        patient组互作权重 = patient_weight
      )
    
    # ---- 合并表格+导出Excel ---
    combine_p_table <- merge(
      p_count_table,
      p_weight_table,
      by = c("源细胞类型", "靶细胞类型"),
      all = TRUE
    )
    
    # 导出Excel
    write_xlsx(combine_p_table, "gg1_gg2_互作网络p值综合表（patient_vs_carrier）.xlsx")  #有问题20260906，这里p值做的有问题，应该用每个样本人分别跑，不能混在一起或者混在同一个组
    cat("表格导出成功！\n")
    
    
    
    ## 细胞间互作次数网络图
    # 红色为STIM组相比carrier组互作次数和互作强度增加，线越粗表示差异越大；蓝色则表示减少。
    # 两个数据集之间的细胞-细胞通信网络中的交互或交互强度的差异数量可以使用圆形图来可视化，其中红色边表示第二个数据集相比于第一个数据集增加的信号。
    
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_compare_interaction_network_patient&carrier.pdf", width = 14, height = 10)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_compare_interaction_heatmap_patient&carrier.pdf", width = 14, height = 10)
    #数量与强度差异热图(主要看互作强度，看质不看量)
    par(mfrow = c(1,1))
    h1 <- netVisual_heatmap(cellchat)
    h2 <- netVisual_heatmap(cellchat, measure = "weight")
    h1+h2
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    
    #保守和特异性信号通路的识别与可视化
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
    gg_all <- gg1 + gg2
    gg_all
    # 直接保存
    ggsave("CellChat_rankNet_comparison_patient&carrier.pdf", gg_all, width = 14, height = 6, dpi = 300)
    dev.off()
    
    
    
    diff.count <- cellchat@net$patient$count - cellchat@net$carrier$count
    write.csv(cellchat@net$patient$count, "output_STIMcount_patient&carrier.csv", quote = F)
    write.csv(cellchat@net$carrier$count, "output_carriercount_patient&carrier.csv", quote = F)
    
    library(pheatmap)
    #打开绘图设备，指定保存路径+尺寸进行保存
    pdf("CellChat_diff_count_heatmap_patient&carrier.pdf", width = 14, height = 10)
    pheatmap(diff.count,
             treeheight_row = "0",treeheight_col = "0",#不画树
             cluster_rows=T, 
             cluster_cols=T)
    # 关闭绘图设备，文件才会生成
    dev.off()
    
    
    # cellchat里有丰富的“受体-配体”分析结果。结合自己的假设，设计合适的图形展示数据。
    View(cellchat)
    View(cellchat@LR[["patient"]][["LRsig"]])
    
  }
  
  
}





####13.7、按样本重新进行cellchat分析，计算P值，再按组对比p值####
{  
  #====读取数据，计算=========================================================

  #===
  # 完整流程：单样本CellChat -> 置换检验组间差异P值 -> 导出表格+绘图
  # 输入：RNA_harmony_annotation.rds（harmony校正+细胞注释完成seurat对象）
  # 适配meta：orig.ident=样本ID；group=control/carrier/patient；celltype=细胞类型
  #===
  library(Seurat)
  library(CellChat)
  library(tidyverse)
  library(writexl)
  library(future)
  
  #===【参数区，已按你的数据修改】===
  meta_col_sample   <- "orig.ident"   # 样本ID列 orig.ident
  meta_col_group    <- "group"        # 分组列：control / carrier / patient
  meta_col_celltype <- "celltype"     # 细胞类型注释列
  nboot_single      <- 100            # 单样本computeCommunProb nboot
  trim_val          <- 10
  min_cells_filter  <- 10
  n_perm_test       <- 200            # 置换检验次数，论文>=200，调试100
  future_workers    <- 2
  #===
  
  ##--- 步骤1：读取harmony注释后的seurat对象，内存精简 ---
  message("===== 步骤1：读取RNA_harmony_annotation.rds =====")
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  
  # V5 seurat 合并layers，精简对象
  RNA_harmony_annotation <- JoinLayers(RNA_harmony_annotation, assay = "RNA")
  RNA_harmony_annotation <- DietSeurat(
    RNA_harmony_annotation,
    assays = "RNA",
    dimreducs = NULL,
    graphs = FALSE
    #reductions = NULL,
    #features = NULL
  )
  
  # 提取样本‑分组对应关系，用于后续分组
  sample_group_map <- unique(RNA_harmony_annotation@meta.data[,c(meta_col_sample, meta_col_group)])
  print("===== 样本分组对应表 =====")
  print(sample_group_map)
  
  # 按orig.ident拆分为每个样本seurat子对象
  seurat_sample_list <- SplitObject(RNA_harmony_annotation, split.by = meta_col_sample)
  rm(RNA_harmony_annotation)
  gc(verbose = FALSE, reset = TRUE)
  
  saveRDS(sample_group_map, "sample_group_map.rds")
  
  
  

  ##--- 步骤2：循环每个样本，独立跑CellChat，merge前保存每个样本全部结果 ---
  message("\n===== 步骤2：逐个样本运行CellChat并保存结果 =====")
  CellChatDB <- CellChatDB.human
  lr_db <- CellChatDB$interaction
  lr_gene_list <- unique(c(lr_db$ligand, lr_db$receptor))
  
  cellchat_sample_list <- list()
  
  for(sample_name in names(seurat_sample_list)){
    message(sprintf("-------- 正在处理样本：%s --------", sample_name))
    seu_sub <- seurat_sample_list[[sample_name]]
    
    # 内存精简单个样本seurat，删除无效的reductions参数！！
    seu_sub <- DietSeurat(seu_sub, 
                          assays = "RNA", 
                          dimreducs = NULL, 
                          graphs = FALSE)
    
    # 校验meta列是否存在，方便debug
    if(!meta_col_celltype %in% colnames(seu_sub@meta.data)){
      stop(sprintf("样本 %s 的meta.data不存在列：%s！！", sample_name, meta_col_celltype))
    }
    
    data_input <- GetAssayData(seu_sub, assay = "RNA", layer = "data")
    
    # 直接提取meta，不需要手动改写vector！！ drop=FALSE保证data.frame格式
    meta_df <- seu_sub@meta.data[, meta_col_celltype, drop=FALSE]
    
    # 只保留配受体相关基因，减小矩阵
    overlap_gene <- intersect(rownames(data_input), lr_gene_list)
    data_input <- data_input[overlap_gene, ]
    
    # 构建cellchat对象
    cc <- createCellChat(object = data_input)
    cc <- addMeta(cc, meta = meta_df)
    cc <- setIdent(cc, ident.use = meta_col_celltype)
    cc@DB <- CellChatDB
    
    # ----关键缺失步骤补上---
    cc <- subsetData(cc)
    
    # 并行设置
    future::plan("multisession", workers = future_workers)
    suppressWarnings({ cc <- identifyOverExpressedGenes(cc) })
    gc(verbose = FALSE, reset = TRUE)
    
    cc <- identifyOverExpressedInteractions(cc)
    future::plan("sequential")
    gc(verbose = FALSE, reset = TRUE)
    
    # 通讯网络推断，沿用原有参数 raw.use=T
    cc <- computeCommunProb(cc, raw.use = TRUE, trim = trim_val, nboot = nboot_single)
    cc <- filterCommunication(cc, min.cells = min_cells_filter)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
    
    # === merge前：保存该样本的各类结果 ===
    # 1) 完整cellchat对象
    saveRDS(cc, paste0("cellchat_", sample_name, ".rds"))
    
    # 2) 该样本通讯结果表格（配受体对/通路/概率/p值），csv+rds双格式
    comm_tbl <- tryCatch(subsetCommunication(cc), error = function(e) NULL)
    if(!is.null(comm_tbl)){
      write.csv(comm_tbl, paste0("communication_", sample_name, ".csv"), row.names = FALSE)
      saveRDS(comm_tbl, paste0("communication_", sample_name, ".rds"))
    }
    
    # 3) 该样本net矩阵（count互作数、weight通讯强度）单独存rds
    saveRDS(cc@net$count,  paste0("net_count_",  sample_name, ".rds"))
    saveRDS(cc@net$weight, paste0("net_weight_", sample_name, ".rds"))
    
    message(sprintf("样本 %s 完成，已保存：cellchat / communication / net_count / net_weight", sample_name))
    
    cellchat_sample_list[[sample_name]] <- cc
    rm(seu_sub, cc, data_input, meta_df)
    gc(verbose = FALSE, reset = TRUE)
  }
  
  
  
  ##--- 步骤2.5：全部样本cellchat对象list打包保存，后续直接读取 ---
  saveRDS(cellchat_sample_list, "all_samples_cellchat_list.rds")
  message("已保存：all_samples_cellchat_list.rds（全部样本cellchat对象list）")
  
  #rm(seurat_sample_list)  可选清除以往缓存
  #gc(verbose = FALSE, reset = TRUE)  可选清除以往缓存
  

  
  #==读取数据，绘图===
  
  # 方式A：读取打包list（推荐，一步到位）
  cellchat_sample_list <- readRDS("all_samples_cellchat_list.rds")
  sample_group_map <- readRDS("sample_group_map.rds")
  
  ##---- 步骤3：按照分组，组装control / carrier / patient 的cellchat list ---
  message("\n===== 步骤3：按分组组装cellchat对象list =====")
  # ① 先定义向量！！
  ctrl_sample_vec    <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]] == "control"]
  carrier_sample_vec <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]] == "carrier"]
  stim_sample_vec    <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]] == "patient"]
  
  # ② 向量定义完成之后，再用来索引得到ctrl_list/carrier_list/stim_list
  ctrl_list    <- cellchat_sample_list[ctrl_sample_vec]
  carrier_list <- cellchat_sample_list[carrier_sample_vec]
  stim_list    <- cellchat_sample_list[stim_sample_vec]
  
  message(sprintf("control组样本数：%d", length(ctrl_list)))
  message(sprintf("carrier组样本数：%d", length(carrier_list)))
  message(sprintf("patient组样本数：%d", length(stim_list)))
  

  ##--- 步骤4：封装置换检验函数：输入两组cellchat list，输出p值表格 ---
  library(CellChat)
  library(tidyverse)
  library(writexl)
  
  #===【参数区】==
  n_perm_test <- 500
  logFC_cut   <- 0.3
  fdr_cut     <- 0.05

  
  ##---- 自实现样本水平置换检验函数（只做统计，无绘图） ---
  run_two_group_perm <- function(list_1, list_2, name_1, name_2, n.perm = 500){
    message(sprintf("\n======== 样本水平置换检验：%s VS %s ========", name_1, name_2))
    
    # 1. 提取每个样本的通讯强度矩阵
    mats1 <- lapply(list_1, function(cc) as.matrix(cc@net$weight))
    mats2 <- lapply(list_2, function(cc) as.matrix(cc@net$weight))
    
    # 2. 统一细胞类型（并集，缺失补0）
    all_types <- sort(unique(c(unlist(lapply(mats1, rownames)),
                               unlist(lapply(mats2, rownames)))))
    k <- length(all_types)
    
    standardize_mat <- function(mat){
      m <- matrix(0, nrow = k, ncol = k, dimnames = list(all_types, all_types))
      idx <- intersect(rownames(mat), all_types)
      m[idx, idx] <- mat[idx, idx, drop = FALSE]
      m
    }
    mats1 <- lapply(mats1, standardize_mat)
    mats2 <- lapply(mats2, standardize_mat)
    
    n1 <- length(mats1); n2 <- length(mats2)
    if(n1 == 0 || n2 == 0) stop("其中一组样本数量为0！")
    all_mats <- c(mats1, mats2)
    
    # 3. 真实观测：两组均值差
    mean1 <- Reduce(`+`, mats1) / n1
    mean2 <- Reduce(`+`, mats2) / n2
    obs_diff <- mean1 - mean2
    
    # 4. 置换检验（固定种子，可复现）
    set.seed(2024)
    n_total <- n1 + n2
    perm_count <- array(0, dim = c(k, k))
    
    for(b in 1:n.perm){
      perm_idx <- sample(n_total)
      g1_idx <- perm_idx[1:n1]
      g2_idx <- perm_idx[(n1+1):n_total]
      perm_mean1 <- Reduce(`+`, all_mats[g1_idx]) / n1
      perm_mean2 <- Reduce(`+`, all_mats[g2_idx]) / n2
      perm_diff <- perm_mean1 - perm_mean2
      perm_count <- perm_count + (abs(perm_diff) >= abs(obs_diff))
    }
    
    # 5. p值（+1修正避免0）+ BH校正
    pval <- (perm_count + 1) / (n.perm + 1)
    fdr_mat <- matrix(p.adjust(as.vector(pval), method = "BH"),
                      nrow = k, dimnames = list(all_types, all_types))
    
    # 6. 长表格 + log2FC
    tbl <- expand.grid(source = all_types, target = all_types, stringsAsFactors = FALSE) %>%
      mutate(
        p_perm   = as.vector(pval),
        fdr_BH   = as.vector(fdr_mat),
        !!sym(paste0("mean_",name_1,"_prob")) := as.vector(mean1),
        !!sym(paste0("mean_",name_2,"_prob")) := as.vector(mean2)
      ) %>%
      mutate(log2FC_prob = log2( (!!sym(paste0("mean_",name_2,"_prob")) + 1e-6) /
                                   (!!sym(paste0("mean_",name_1,"_prob")) + 1e-6) ))
    
    # 7. 导出Excel完整结果表
    out_file <- sprintf("CellChat_perm_%s_VS_%s.xlsx", name_1, name_2)
    write_xlsx(tbl, out_file)
    message("完整置换表格已输出：", out_file)
    
    # 返回统计对象，保存RDS，留存矩阵方便后续随时再分析/画图
    res <- list(pval = pval,
                fdr = fdr_mat,
                meanProb1 = mean1,
                meanProb2 = mean2,
                obs_diff = obs_diff,
                table = tbl,
                cell_types = all_types)
    saveRDS(res, sprintf("perm_stat_%s_VS_%s.rds", name_1, name_2))
    message("统计中间对象已保存：", sprintf("perm_stat_%s_VS_%s.rds", name_1, name_2))
    
    return(res)
  }
  
  ##--- 执行三组两两置换检验 ---
  message("\n===== 步骤4：执行三组两两置换检验（仅统计，无绘图） =====")
  res_ctrl_carrier  <- run_two_group_perm(ctrl_list, carrier_list, "control", "carrier", n.perm = n_perm_test)
  res_ctrl_stim     <- run_two_group_perm(ctrl_list, stim_list,    "control", "patient", n.perm = n_perm_test)
  res_carrier_stim  <- run_two_group_perm(carrier_list, stim_list, "carrier", "patient", n.perm = n_perm_test)
  
  ##--- 筛选显著差异通讯对并导出 ---
  message("\n===== 步骤5：筛选差异互作对 =====")
  sig_ctrl_carrier <- res_ctrl_carrier$table %>%
    filter(fdr_BH < fdr_cut, abs(log2FC_prob) > logFC_cut) %>%
    arrange(fdr_BH)
  sig_ctrl_stim    <- res_ctrl_stim$table %>%
    filter(fdr_BH < fdr_cut, abs(log2FC_prob) > logFC_cut) %>%
    arrange(fdr_BH)
  sig_carrier_stim <- res_carrier_stim$table %>%
    filter(fdr_BH < fdr_cut, abs(log2FC_prob) > logFC_cut) %>%
    arrange(fdr_BH)
  
  write_xlsx(sig_ctrl_carrier, "sig_interactions_control_VS_carrier.xlsx")
  write_xlsx(sig_ctrl_stim,    "sig_interactions_control_VS_patient.xlsx")
  write_xlsx(sig_carrier_stim, "sig_interactions_carrier_VS_patient.xlsx")
  
  message(sprintf("显著互作对：control_vs_carrier=%d 对", nrow(sig_ctrl_carrier)))
  message(sprintf("显著互作对：control_vs_patient=%d 对", nrow(sig_ctrl_stim)))
  message(sprintf("显著互作对：carrier_vs_patient=%d 对", nrow(sig_carrier_stim)))
  
  # 保存全部结果list，后续R会话直接readRDS恢复全部统计结果，不用重跑置换
  all_perm_results <- list(
    ctrl_carrier = res_ctrl_carrier,
    ctrl_stim = res_ctrl_stim,
    carrier_stim = res_carrier_stim
  )
  saveRDS(all_perm_results, "all_permutation_results.rds")
  message("\n✅ 全部置换统计结果保存完成：all_permutation_results.rds")
  message("后续读取：all_perm_results <- readRDS('all_permutation_results.rds')")
  
  message("\n######## 统计流程结束（跳过全部绘图） ########")
  
  
  
  
  
  
  
  #===
  # 配受体（LR）分子水平样本置换检验
  # 输入：ctrl_list / carrier_list / stim_list（每个元素是一个已跑完computeCommunProb的cellchat对象）
  # 输出：全量LR置换表格xlsx + 统计对象rds
  # 依赖：CellChat, tidyverse, writexl
  #===
  library(CellChat)
  library(tidyverse)
  library(writexl)
  
  run_LR_perm <- function(list_1, list_2, name_1, name_2, n.perm = 500){
    message(sprintf("\n======== LR分子水平样本置换检验：%s VS %s ========", name_1, name_2))
    
    ## --- 1. 提取每个样本的三维通讯概率数组 [source, target, LR] ---
    arrs1 <- lapply(list_1, function(cc) cc@net$prob)
    arrs2 <- lapply(list_2, function(cc) cc@net$prob)
    if(any(sapply(arrs1, is.null)) || any(sapply(arrs2, is.null))){
      stop("存在样本 @net$prob 为NULL，请先用完整性校验脚本确认！")
    }
    
    ## --- 2. 统一维度（并集，缺失补0） ---
    all_src <- sort(unique(c(unlist(lapply(arrs1, function(a) dimnames(a)[[1]])),
                             unlist(lapply(arrs2, function(a) dimnames(a)[[1]])))))
    all_tgt <- sort(unique(c(unlist(lapply(arrs1, function(a) dimnames(a)[[2]])),
                             unlist(lapply(arrs2, function(a) dimnames(a)[[2]])))))
    all_lr  <- sort(unique(c(unlist(lapply(arrs1, function(a) dimnames(a)[[3]])),
                             unlist(lapply(arrs2, function(a) dimnames(a)[[3]])))))
    nS <- length(all_src); nT <- length(all_tgt); nL <- length(all_lr)
    message(sprintf("细胞类型：%d 种source × %d 种target；LR分子总数：%d", nS, nT, nL))
    
    std_arr <- function(a){
      out <- array(0, dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
      s <- intersect(dimnames(a)[[1]], all_src)
      t <- intersect(dimnames(a)[[2]], all_tgt)
      l <- intersect(dimnames(a)[[3]], all_lr)
      if(length(s) > 0 && length(t) > 0 && length(l) > 0){
        out[s, t, l] <- a[s, t, l, drop = FALSE]
      }
      out
    }
    arrs1 <- lapply(arrs1, std_arr)
    arrs2 <- lapply(arrs2, std_arr)
    
    n1 <- length(arrs1); n2 <- length(arrs2)
    if(n1 == 0 || n2 == 0) stop("其中一组样本数量为0！")
    n_total <- n1 + n2
    all_arr <- c(arrs1, arrs2)
    
    ## --- 3. 压平成矩阵加速：每列一个样本，行 = nS×nT×nL ---
    big_mat <- sapply(all_arr, function(a) as.vector(a))   # 维度 [nS*nT*nL, n_total]
    message(sprintf("压平矩阵维度：%d 行 × %d 个样本，开始置换检验...", nrow(big_mat), ncol(big_mat)))
    
    ## --- 4. 真实观测：两组均值差 ---
    mean1_vec <- rowMeans(big_mat[, 1:n1, drop = FALSE])
    mean2_vec <- rowMeans(big_mat[, (n1+1):n_total, drop = FALSE])
    obs_diff_vec <- mean1_vec - mean2_vec
    
    ## --- 5. 样本水平置换 ---
    set.seed(2024)
    cnt <- numeric(nrow(big_mat))
    
    for(b in 1:n.perm){
      perm_idx <- sample(n_total)
      g1 <- perm_idx[1:n1]
      g2 <- perm_idx[(n1+1):n_total]
      p1 <- rowMeans(big_mat[, g1, drop = FALSE])
      p2 <- rowMeans(big_mat[, g2, drop = FALSE])
      pdiff <- p1 - p2
      cnt <- cnt + (abs(pdiff) >= abs(obs_diff_vec))
      if(b %% 100 == 0) message(sprintf("  置换进度：%d / %d", b, n.perm))
    }
    
    ## --- 6. p值（+1修正）+ BH校正 ---
    pval_vec <- (cnt + 1) / (n.perm + 1)
    fdr_vec  <- p.adjust(pval_vec, method = "BH")
    
    # 恢复成三维数组
    pval_arr  <- array(pval_vec,  dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
    fdr_arr   <- array(fdr_vec,   dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
    mean1_arr <- array(mean1_vec, dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
    mean2_arr <- array(mean2_vec, dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
    
    ## --- 7. 长表格 ---
    tbl <- expand.grid(source = all_src, target = all_tgt, LR = all_lr, stringsAsFactors = FALSE) %>%
      mutate(
        p_perm = pval_vec,
        fdr_BH = fdr_vec,
        !!sym(paste0("mean_", name_1, "_prob")) := mean1_vec,
        !!sym(paste0("mean_", name_2, "_prob")) := mean2_vec
      ) %>%
      mutate(log2FC_prob = log2( (!!sym(paste0("mean_", name_2, "_prob")) + 1e-6) /
                                   (!!sym(paste0("mean_", name_1, "_prob")) + 1e-6) ))
    
    ## --- 8. 保存输出 ---
    out_xlsx <- sprintf("CellChat_LRlevel_perm_%s_VS_%s.xlsx", name_1, name_2)
    write_xlsx(tbl, out_xlsx)
    message("✅ 全量LR置换表格输出：", out_xlsx)
    
    res <- list(
      pval = pval_arr, fdr = fdr_arr,
      meanProb1 = mean1_arr, meanProb2 = mean2_arr,
      table = tbl,
      celltype_source = all_src, celltype_target = all_tgt, LR_names = all_lr
    )
    saveRDS(res, sprintf("perm_LRstat_%s_VS_%s.rds", name_1, name_2))
    message("✅ 统计对象保存：", sprintf("perm_LRstat_%s_VS_%s.rds", name_1, name_2))
    return(res)
  }
  
  #===
  # 执行三组两两LR水平置换检验
  # ⚠️ 计算量大：n.perm=500，三组会跑较久，建议服务器执行
  #===
  message("\n===== 执行三组LR水平置换检验 =====")
  res_LR_ctrl_carrier <- run_LR_perm(ctrl_list, carrier_list, "control", "carrier", n.perm = 500)
  res_LR_ctrl_stim    <- run_LR_perm(ctrl_list, stim_list,    "control", "patient", n.perm = 500)
  res_LR_carrier_stim <- run_LR_perm(carrier_list, stim_list, "carrier", "patient", n.perm = 500)
  

  # （可选）快速查看结果概况
  message("\n===== 结果概况 =====")
  for(res in list(res_LR_ctrl_carrier, res_LR_ctrl_stim, res_LR_carrier_stim)){
    tb <- res$table
    message(sprintf("总LR检验数：%d；p<0.05：%d；FDR<0.1：%d；FDR<0.2：%d",
                    nrow(tb),
                    sum(tb$p_perm < 0.05),
                    sum(tb$fdr_BH < 0.1),
                    sum(tb$fdr_BH < 0.2)))
  }
  
  
  
  
  
  #==== 按样本单独分别提取配受体互作数量和强度、作图 ===
  library(CellChat)
  library(tidyverse)
  library(writexl)
  
  #==== 1.读取原始结果，拆分分组 ===
  all_samples_cellchat_list <- readRDS("all_samples_cellchat_list.rds")
  sample_group_map <- readRDS("sample_group_map.rds")
  
  meta_col_sample   <- "orig.ident"
  meta_col_group    <- "group"
  
  # 拆分三组list
  ctrl_sample_vec    <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]]=="control"]
  carrier_sample_vec <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]]=="carrier"]
  stim_sample_vec    <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]]=="patient"]
  
  ctrl_list    <- all_samples_cellchat_list[ctrl_sample_vec]
  carrier_list <- all_samples_cellchat_list[carrier_sample_vec]
  stim_list    <- all_samples_cellchat_list[stim_sample_vec]
  
  #=== 2.提取每个样本的count、strength ===
  extract_sample_comm_summary <- function(cc_list){
    res <- data.frame(
      sample_name = character(),
      total_LR_count = integer(),
      total_strength = numeric(),
      stringsAsFactors = FALSE
    )
    for(samp in names(cc_list)){
      cc <- cc_list[[samp]]
      cnt_sum <- sum(cc@net$count, na.rm = TRUE)
      wt_sum  <- sum(cc@net$weight, na.rm = TRUE)
      res <- rbind(res, data.frame(
        sample_name = samp,
        total_LR_count = as.integer(cnt_sum),
        total_strength = wt_sum
      ))
    }
    return(res)
  }
  
  # 提取全部18个样本
  sample_summary_df <- extract_sample_comm_summary(all_samples_cellchat_list)
  
  # ===修复点===
  # sample_group_map原始两列： orig.ident , group
  colnames(sample_group_map) <- c("sample_name", "group")
  
  # 拼接分组信息
  sample_summary_df <- left_join(sample_summary_df, sample_group_map, by = "sample_name")
  
  # 调整分组因子顺序，绘图顺序：control → carrier → patient
  sample_summary_df$group <- factor(sample_summary_df$group, levels = c("control","carrier","patient"))
  
  # 输出表格，保存每个样本数值
  write_xlsx(sample_summary_df, "per_sample_communication_summary.xlsx")
  message("✅已输出每个样本指标表格：per_sample_communication_summary.xlsx")
  print(sample_summary_df, row.names = FALSE)
  
  #===3.绘图：箱线图+散点（展示每个样本真实分布）===
  # ① Number of inferred interactions（LR总数量）
  p1 <- ggplot(sample_summary_df, aes(x = group, y = total_LR_count, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    labs(title = "Number of inferred interactions (per sample)",
         x = "Group", y = "Number of inferred LR interactions") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  # ② Interaction strength（总通讯强度）
  p2 <- ggplot(sample_summary_df, aes(x = group, y = total_strength, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    labs(title = "Interaction strength (per sample)",
         x = "Group", y = "Aggregated interaction strength") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  # 拼图输出pdf
  p_combine <- p1 + p2 + patchwork::plot_layout(ncol = 2)
  ggsave("per_sample_interaction_box_jitter.pdf", plot = p_combine, width = 12, height = 5, dpi = 300)
  message("✅已输出图：per_sample_interaction_box_jitter.pdf")
  
  #===（可选）做组间wilcoxon秩和检验，全局指标的简单统计 ===
  message("\n===== 全局指标组间Wilcoxon检验（仅全局聚合指标，不等于细胞‑细胞/LR置换检验） =====")
  
  # control vs carrier
  wc_ctrl_carrier_count <- wilcox.test(total_LR_count ~ group,
                                       data = filter(sample_summary_df, group %in% c("control","carrier")))
  wc_ctrl_carrier_strength <- wilcox.test(total_strength ~ group,
                                          data = filter(sample_summary_df, group %in% c("control","carrier")))
  
  # control vs patient
  wc_ctrl_patient_count <- wilcox.test(total_LR_count ~ group,
                                       data = filter(sample_summary_df, group %in% c("control","patient")))
  wc_ctrl_patient_strength <- wilcox.test(total_strength ~ group,
                                          data = filter(sample_summary_df, group %in% c("control","patient")))
  
  # carrier vs patient
  wc_carrier_patient_count <- wilcox.test(total_LR_count ~ group,
                                          data = filter(sample_summary_df, group %in% c("carrier","patient")))
  wc_carrier_patient_strength <- wilcox.test(total_strength ~ group,
                                             data = filter(sample_summary_df, group %in% c("carrier","patient")))
  
  cat("\ncontrol‑carrier count p‑value:", wc_ctrl_carrier_count$p.value)
  cat("\ncontrol‑carrier strength p‑value:", wc_ctrl_carrier_strength$p.value)
  cat("\ncontrol‑patient count p‑value:", wc_ctrl_patient_count$p.value)
  cat("\ncontrol‑patient strength p‑value:", wc_ctrl_patient_strength$p.value)
  cat("\ncarrier‑patient count p‑value:", wc_carrier_patient_count$p.value)
  cat("\ncarrier‑patient strength p‑value:", wc_carrier_patient_strength$p.value)
  
 
  
}
  
  



#####14.1、patient VS control组间基因表达差异分析####
{
  #====读取数据，计算=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  
  # 相同细胞类型在不同分组中的差异基因（FindMarkers）
  seu_harmony = RNA_harmony_annotation
  table(seu_harmony$celltype,seu_harmony$group)
  # 如果某一组少于3个细胞，不能做分组间差异分析，需要排除掉此细胞类型做
  # seu_harmony <- subset(seu_harmony, subset = celltype != "C6-lymphocyte")
  seu_harmony$celltype = as.character(seu_harmony$celltype)
  Idents(seu_harmony)="celltype"
  table(seu_harmony$celltype)
  ## 差异比较的是不同分组中相同细胞类型之间的差异基因
  type=unique(seu_harmony$celltype)
  table(seu_harmony$group)
  
  
  #寻找组间所有类型细胞的差异基因#
  r.deg=data.frame()
  for (i in 1:length(type)) {  # 遍历 type 中每一种细胞类型
    deg = FindMarkers(seu_harmony, ident.1 = "patient", ident.2 = "control",  # 使用 FindMarkers 进行差异分析
                      group.by = "group", subset.ident = type[i],min.pct = 0.25)  # 按疾病分组，并针对 type[i] 定义的细胞群体进行子集分析
    
    #write.csv(deg, file = paste0("Ctrl_patient_组间差异分析/", type[i], 'deg.csv'))  # 将每次差异分析的结果保存为单独的 CSV 文件
    
    ### 核心修正：将行名（基因名）转为显式列，避免行名重复导致后面使用rbind合并时使原基因名因重复而被增加尾缀！！！
    deg$gene <- rownames(deg)  # 新增gene列，存储原始标准基因名
    rownames(deg) <- NULL      # 清空行名，让R自动生成数字行名
    # 保存单个细胞类型的结果（gene列保留标准名）
    write.csv(deg, file = paste0("Ctrl_patient_组间差异分析/", type[i], 'deg.csv'), row.names = FALSE)
    
    deg$celltype = type[i]  # 给差异分析结果添加列，记录当前分析的细胞类型
    deg$unm = i - 1  # 给差异分析结果添加列，记录当前循环的索引
    r.deg = rbind(deg, r.deg)  # 将当前分析的结果与之前的结果合并，生成一个综合的差异基因数据表
  }
  table(r.deg$celltype)
  ###ident.1代表处理组
  ###ident.2代表对照组
  ###avg_log2FC = log2( ident.1组的基因平均表达量 / ident.2组的基因平均表达量 )，所以上调下调是患者组/对照组。一般红色是上调 即患者比对照上调。蓝色是下调。
  ###使用group.by指明patient，control在哪一列，即根据什么进行的分组
  ###subset.ident指定分析哪个cluster（原为celltype）差异基因，必须是active.ident的分组

  saveRDS(r.deg,"r.deg_filter_markers_patient_VS_control.rds")
  
  write.csv(r.deg,"Ctrl_patient_组间差异分析/r_deg_all_markers_patient_VS_control.csv",row.names = T)
  
  
  ###根据自己计算的marker基因数量确定log2FC的阈值，这里先定为1，筛选差异基因
  ###提取p<0.05&logFC>1
  s.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) #调整至0.5放宽，与第三个组间差异标准一致
  table(s.deg$celltype)
  ###标记出上调和下调并转化为因子
  s.deg$threshold <- as.factor(ifelse(s.deg$avg_log2FC > 1 , 'Up', 'Down'))
  table(s.deg$threshold)
  dim(s.deg)
  ###标记出显著和不显著并转化为因子
  s.deg$adj_p_signi <- as.factor(ifelse(s.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
  s.deg$thr_signi <- paste0(s.deg$threshold, "_", s.deg$adj_p_signi)
  
  saveRDS(s.deg,"s.deg_filter_markers_patient_VS_control.rds")
  
  write.csv(s.deg,"Ctrl_patient_组间差异分析/s_deg_filter_patient_VS_control.csv",row.names = T)
  
  
  #====读取数据，绘图=========================================================
  r.deg <- readRDS("r.deg_filter_markers_patient_VS_control.rds")
  s.deg <- readRDS("s.deg_filter_markers_patient_VS_control.rds")
  
  ###差异基因绘图
  ### 自定义显示想要展示的基因名，这里挑选log2FC为top5的基因进行展示
  ###每个celltype找前5个上调基因
  top_up_label <- s.deg %>% 
    subset(., threshold%in%"Up") %>% 
    group_by(celltype) %>% 
    top_n(n = 5, wt = avg_log2FC) %>% 
    as.data.frame()
  ###每个celltype找前5个下调基因
  top_down_label <- s.deg %>% 
    subset(., threshold %in% "Down") %>% 
    group_by(celltype) %>% 
    top_n(n = -5, wt = avg_log2FC) %>% 
    as.data.frame()
  ###将上调和下调进行合并
  top_label <- rbind(top_up_label,top_down_label)
  
  library(scRNAtoolVis)
  library(jjPlot)
  library(ggrepel)
  
  colors <- c("red", "blue", "green", "yellow", "purple", "orange", "pink", "cyan", "brown", "black")
  
  ### 画图展示jjVolcano
  r.deg$cluster <- r.deg$celltype #jjVolcano函数只能默认识别cluster列，故复制新增一列存放细胞类型
  jjVolcano(diffData =r.deg, 
            tile.col = colors[1:11],
            pSize = 0.4,###调节点的大小
            legend.position=c(0.1,0.9),###调节图注坐标位置
            celltypeSize=2,###调节细胞celltype文字大小
            topGeneN=5)+
    labs(title = "patient VS control")  # 添加标题
  ggsave(
    filename = "Ctrl_patient_组间差异分析/jjVolcano_patient_VS_control.pdf",
    plot = last_plot(),  # 关键：指定要保存的图（last_plot()获取最后生成的图）
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300  # 可选：增加分辨率，避免图模糊
  )
  
  #scRNAtoolVis::markerVolcano()
  
  # 绘制火山图2            
  markerVolcano(            
    markers = r.deg,            
    topn = 5,        
    labelCol = ggsci::pal_npg()(11)
  )
  ggsave(
    filename = "Ctrl_patient_组间差异分析/Volcano_patient_VS_control.pdf",
    plot = last_plot(),  # 关键：指定要保存的图（last_plot()获取最后生成的图）
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300  # 可选：增加分辨率，避免图模糊
  )
  
  
  sub_s.deg <- subset(s.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
  ##绘制韦恩图
  # 安装（仅首次需要）
  #install.packages("VennDiagram")
  # 加载包
  library(VennDiagram)
  # 安装（仅首次需要）
  #install.packages("ggvenn")
  # 加载包
  library(ggvenn)
  library(ggplot2)  # 依赖ggplot2
  
  # T cells 基因集：严格筛选（p_val_adj<0.05，|log2FC|>0.5）
  # 筛选特定细胞类型的行（逻辑索引：TRUE表示符合条件）
  filter_rows <- sub_s.deg$celltype == "T cell"
  # 提取基因列+去重+去NA
  degs_list1 <- unique(sub_s.deg$gene[filter_rows])
  degs_list1 <- degs_list1[!is.na(degs_list1)]
  # 查看结果
  print(degs_list1)
  
  # B cells 基因集：严格筛选（p_val_adj<0.05，|log2FC|>0.5）
  # 筛选特定细胞类型的行（逻辑索引：TRUE表示符合条件）
  filter_rows <- sub_s.deg$celltype == "B cell"
  # 提取基因列+去重+去NA
  degs_list2 <- unique(sub_s.deg$gene[filter_rows])
  degs_list2 <- degs_list2[!is.na(degs_list2)]
  # 查看结果
  print(degs_list2)
  
  # NK cells 基因集：严格筛选（p_val_adj<0.05，|log2FC|>0.5）
  # 筛选特定细胞类型的行（逻辑索引：TRUE表示符合条件）
  filter_rows <- sub_s.deg$celltype == "NK cell"
  # 提取基因列+去重+去NA
  degs_list3 <- unique(sub_s.deg$gene[filter_rows])
  degs_list3 <- degs_list3[!is.na(degs_list3)]
  # 查看结果
  print(degs_list3)
  
  # monocytes cells 基因集：严格筛选（p_val_adj<0.05，|log2FC|>0.5）
  # 筛选特定细胞类型的行（逻辑索引：TRUE表示符合条件）
  filter_rows <- sub_s.deg$celltype == "Monocyte"
  # 提取基因列+去重+去NA
  degs_list4 <- unique(sub_s.deg$gene[filter_rows])
  degs_list4 <- degs_list4[!is.na(degs_list4)]
  # 查看结果
  print(degs_list4)
  
  
  # 2. 整理为基因集列表（命名方便后续展示）
  gene_sets <- list(
    T_cell = degs_list1,  # 命名1
    B_cell = degs_list2,   # 命名2
    NK_cell = degs_list3,    # 命名3
    Monocytes = degs_list4  # 命名4
  )
  
  
  # 2. 提取你的4个基因集（无需修改，直接对应gene_sets）
  g1 <- gene_sets[["T_cell"]]
  g2 <- gene_sets[["B_cell"]]
  g3 <- gene_sets[["NK_cell"]]
  g4 <- gene_sets[["Monocytes"]]
  
  # 3. 自动计算4个基因集的所有交集区域数量（无需手动输入，避免出错）
  area1 <- length(g1)
  area2 <- length(g2)
  area3 <- length(g3)
  area4 <- length(g4)
  
  n12 <- length(intersect(g1, g2))
  n13 <- length(intersect(g1, g3))
  n14 <- length(intersect(g1, g4))
  n23 <- length(intersect(g2, g3))
  n24 <- length(intersect(g2, g4))
  n34 <- length(intersect(g3, g4))
  
  n123 <- length(intersect(intersect(g1, g2), g3))
  n124 <- length(intersect(intersect(g1, g2), g4))
  n134 <- length(intersect(intersect(g1, g3), g4))
  n234 <- length(intersect(intersect(g2, g3), g4))
  
  n1234 <- length(intersect(intersect(intersect(g1, g2), g3), g4))
  
  # 4. 绘制4细胞类型基因集韦恩图（标签自动使用你的基因集名称）
  venn_plot <- draw.quad.venn(
    area1 = area1,
    area2 = area2,
    area3 = area3,
    area4 = area4,
    n12 = n12,
    n13 = n13,
    n14 = n14,
    n23 = n23,
    n24 = n24,
    n34 = n34,
    n123 = n123,
    n124 = n124,
    n134 = n134,
    n234 = n234,
    n1234 = n1234,
    category = names(gene_sets),  # 自动显示T_cell/B_cell/NK_cell/monocytes_cell
    fill = c("#FF6B6B", "#4ECDC4", "#7927D1", "#4657D1"),  # 对应4个细胞类型颜色
    alpha = 0.6,                  # 透明度
    lwd = 1,                      # 边框线宽
    fontsize = 5,                 # 交集数量文字大小
    cat.fontsize = 6,             # 细胞类型标签文字大小
    cat.fontface = "bold"         # 标签加粗，更清晰
  )
  
  # 5. 显示韦恩图
  grid.draw(venn_plot)
  
  # 6. 保存高清图片（两种格式可选，满足论文投稿要求）
  # TIFF格式（无损高清）
  tiff("Ctrl_patient_组间差异分析/patient_VS_control_4_cell_types_venn_VennDiagram.tiff", width = 11, height = 10, units = "in", res = 300)
  grid.draw(venn_plot)
  dev.off()

  
}


#####14.2、carrier VS control组间基因表达差异分析####
{
  # 提前加载所有依赖包
  library(Seurat)
  library(dplyr)
  library(scRNAtoolVis)
  library(jjPlot)
  library(ggrepel)
  library(ggsci)
  library(VennDiagram)
  library(ggvenn)
  library(ggplot2)
  
  
#====读取数据，计算=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  seu_harmony = RNA_harmony_annotation
  table(seu_harmony$celltype, seu_harmony$group)
  
  # 1. 最先创建输出文件夹，解决循环写入报错
  out_dir <- "Ctrl_carrier_组间差异分析"
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  seu_harmony$celltype = as.character(seu_harmony$celltype)
  Idents(seu_harmony) = "celltype"
  type = unique(seu_harmony$celltype)
  table(seu_harmony$group)
  
  
  #寻找组间所有类型细胞的差异基因#
  r.deg=data.frame()
  for (i in 1:length(type)) {
    deg = FindMarkers(seu_harmony, ident.1 = "carrier", ident.2 = "control",
                      group.by = "group", subset.ident = type[i], min.pct = 0.25)
    
    deg$gene <- rownames(deg)
    rownames(deg) <- NULL
    # 单类细胞DEG输出
    write.csv(deg, file = paste0(out_dir, "/", type[i], 'deg.csv'), row.names = FALSE)
    
    deg$celltype = type[i]
    deg$unm = i - 1
    # 仅当当前细胞有差异基因才合并，避免空数据框
    if(nrow(deg) > 0){
      r.deg = rbind(deg, r.deg)
    }
  }
  table(r.deg$celltype)
  cat("carrier组总差异基因数量：", nrow(r.deg), "\n")
  
  ###ident.1=carrier组，ident.2=control组
  ###avg_log2FC>0：carrier中上调；<0：carrier下调
  
  # rds文件名加carrier区分，不会和患者组混淆
  saveRDS(r.deg, "r.deg_filter_markers_carrier_VS_control.rds")
  write.csv(r.deg, paste0(out_dir,"/r_deg_all_markers_carrier_VS_control.csv"), row.names = T)
  
  
  ###筛选显著差异基因 p_val_adj < 0.05 & |avg_log2FC| > 1
  s.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
  table(s.deg$celltype)
  s.deg$threshold <- as.factor(ifelse(s.deg$avg_log2FC > 1 , 'Up', 'Down'))
  table(s.deg$threshold)
  dim(s.deg)
  s.deg$adj_p_signi <- as.factor(ifelse(s.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
  s.deg$thr_signi <- paste0(s.deg$threshold, "_", s.deg$adj_p_signi)
  
  # 加carrier后缀区分文件
  saveRDS(s.deg, "s.deg_filter_markers_carrier_VS_control.rds")
  write.csv(s.deg, paste0(out_dir,"/s_deg_filter_carrier_VS_control.csv"), row.names = T)
  

  
#====读取数据，绘图=========================================================
  r.deg <- readRDS("r.deg_filter_markers_carrier_VS_control.rds")
  s.deg <- readRDS("s.deg_filter_markers_carrier_VS_control.rds")
  
  # 提取top5上下调标注基因
  top_up_label <- s.deg %>% 
    subset(., threshold%in%"Up") %>% 
    group_by(celltype) %>% 
    top_n(n = 5, wt = avg_log2FC) %>% 
    as.data.frame()
  top_down_label <- s.deg %>% 
    subset(., threshold %in% "Down") %>% 
    group_by(celltype) %>% 
    top_n(n = -5, wt = avg_log2FC) %>% 
    as.data.frame()
  top_label <- rbind(top_up_label,top_down_label)
  
  colors <- c("red", "blue", "green", "yellow", "purple", "orange", "pink", "cyan", "brown", "black")
  
  # 多细胞整合火山图
  r.deg$cluster <- r.deg$celltype
  jjVolcano(diffData =r.deg, 
            tile.col = colors[1:11],
            pSize = 0.4,
            legend.position=c(0.1,0.9),
            celltypeSize=2,
            topGeneN=5)+
    labs(title = "carrier VS control")
  ggsave(
    filename = paste0(out_dir,"/jjVolcano_carrier_VS_control.pdf"),
    plot = last_plot(),
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300
  )
  
  # 单图火山图
  markerVolcano(            
    markers = r.deg,            
    topn = 5,        
    labelCol = ggsci::pal_npg()(11)
  )
  ggsave(
    filename = paste0(out_dir,"/Volcano_carrier_VS_control.pdf"),
    plot = last_plot(),
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300
  )
  
  
  sub_s.deg <- subset(s.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
  ##绘制韦恩图
  # T cells
  filter_rows <- sub_s.deg$celltype == "T cell"
  degs_list1 <- unique(sub_s.deg$gene[filter_rows])
  degs_list1 <- degs_list1[!is.na(degs_list1)]
  print(degs_list1)
  
  # B cells
  filter_rows <- sub_s.deg$celltype == "B cell"
  degs_list2 <- unique(sub_s.deg$gene[filter_rows])
  degs_list2 <- degs_list2[!is.na(degs_list2)]
  print(degs_list2)
  
  # NK cells
  filter_rows <- sub_s.deg$celltype == "NK cell"
  degs_list3 <- unique(sub_s.deg$gene[filter_rows])
  degs_list3 <- degs_list3[!is.na(degs_list3)]
  print(degs_list3)
  
  # Monocytes
  filter_rows <- sub_s.deg$celltype == "Monocyte"
  degs_list4 <- unique(sub_s.deg$gene[filter_rows])
  degs_list4 <- degs_list4[!is.na(degs_list4)]
  print(degs_list4)
  
  gene_sets <- list(
    T_cell = degs_list1,
    B_cell = degs_list2,
    NK_cell = degs_list3,
    Monocytes = degs_list4
  )
  
  g1 <- gene_sets[["T_cell"]]
  g2 <- gene_sets[["B_cell"]]
  g3 <- gene_sets[["NK_cell"]]
  g4 <- gene_sets[["Monocytes"]]
  
  area1 <- length(g1)
  area2 <- length(g2)
  area3 <- length(g3)
  area4 <- length(g4)
  
  n12 <- length(intersect(g1, g2))
  n13 <- length(intersect(g1, g3))
  n14 <- length(intersect(g1, g4))
  n23 <- length(intersect(g2, g3))
  n24 <- length(intersect(g2, g4))
  n34 <- length(intersect(g3, g4))
  
  n123 <- length(intersect(intersect(g1, g2), g3))
  n124 <- length(intersect(intersect(g1, g2), g4))
  n134 <- length(intersect(intersect(g1, g3), g4))
  n234 <- length(intersect(intersect(g2, g3), g4))
  n1234 <- length(intersect(intersect(intersect(g1, g2), g3), g4))
  
  venn_plot <- draw.quad.venn(
    area1 = area1, area2 = area2, area3 = area3, area4 = area4,
    n12 = n12, n13 = n13, n14 = n14, n23 = n23, n24 = n24, n34 = n34,
    n123 = n123, n124 = n124, n134 = n134, n234 = n234, n1234 = n1234,
    category = names(gene_sets),
    fill = c("#FF6B6B", "#4ECDC4", "#7927D1", "#4657D1"),
    alpha = 0.6, lwd = 1, fontsize = 5, cat.fontsize = 6, cat.fontface = "bold"
  )
  
  grid.draw(venn_plot)
  
  # 韦恩图文件带carrier标识
  tiff(paste0(out_dir,"/carrier_VS_control_4_cell_types_venn_VennDiagram.tiff"), width = 11, height = 10, units = "in", res = 300)
  grid.draw(venn_plot)
  dev.off()
  
}

#####14.3、patient VS carrier组间基因表达差异分析####
{
  #====读取数据，计算=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  
  # 相同细胞类型在不同分组中的差异基因（FindMarkers）
  seu_harmony = RNA_harmony_annotation
  table(seu_harmony$celltype,seu_harmony$group)
  # 如果某一组少于3个细胞，不能做分组间差异分析，需要排除掉此细胞类型做
  # seu_harmony <- subset(seu_harmony, subset = celltype != "C6-lymphocyte")
  seu_harmony$celltype = as.character(seu_harmony$celltype)
  Idents(seu_harmony)="celltype"
  table(seu_harmony$celltype)
  ## 差异比较的是不同分组中相同细胞类型之间的差异基因
  type=unique(seu_harmony$celltype)
  table(seu_harmony$group)
  
  
  #寻找组间所有类型细胞的差异基因#
  r.deg=data.frame()
  for (i in 1:length(type)) {  # 遍历 type 中每一种细胞类型
    deg = FindMarkers(seu_harmony, ident.1 = "patient", ident.2 = "carrier",  # 使用 FindMarkers 进行差异分析
                      group.by = "group", subset.ident = type[i],min.pct = 0.25)  # 按疾病分组，并针对 type[i] 定义的细胞群体进行子集分析
    
    #write.csv(deg, file = paste0("carrier_patient_组间差异分析/", type[i], 'deg.csv'))  # 将每次差异分析的结果保存为单独的 CSV 文件
    
    ### 核心修正：将行名（基因名）转为显式列，避免行名重复导致后面使用rbind合并时使原基因名因重复而被增加尾缀！！！
    deg$gene <- rownames(deg)  # 新增gene列，存储原始标准基因名
    rownames(deg) <- NULL      # 清空行名，让R自动生成数字行名
    # 保存单个细胞类型的结果（gene列保留标准名）
    write.csv(deg, file = paste0("carrier_patient_组间差异分析/", type[i], 'deg.csv'), row.names = FALSE)
    
    deg$celltype = type[i]  # 给差异分析结果添加列，记录当前分析的细胞类型
    deg$unm = i - 1  # 给差异分析结果添加列，记录当前循环的索引
    r.deg = rbind(deg, r.deg)  # 将当前分析的结果与之前的结果合并，生成一个综合的差异基因数据表
  }
  table(r.deg$celltype)
  ###ident.1代表处理组
  ###ident.2代表携带者组
  ###avg_log2FC = log2( ident.1组的基因平均表达量 / ident.2组的基因平均表达量 )，所以上调下调是患者组/携带者组。一般红色是上调 即患者比携带者上调。蓝色是下调。
  ###使用group.by指明patient，carrier在哪一列，即根据什么进行的分组
  ###subset.ident指定分析哪个cluster（原为celltype）差异基因，必须是active.ident的分组

  saveRDS(r.deg,"r.deg_filter_markers_patient_VS_carrier.rds")
  
  write.csv(r.deg,"carrier_patient_组间差异分析/r_deg_all_markers_patient_VS_carrier.csv",row.names = T)
  
  
  ###根据自己计算的marker基因数量确定log2FC的阈值，这里先定为1，筛选差异基因
  ###提取p<0.05&logFC>1
  s.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)  #由于后续由1筛选的差异基因过少，富集不出来，故此组放宽至0.5
  table(s.deg$celltype)
  ###标记出上调和下调并转化为因子
  s.deg$threshold <- as.factor(ifelse(s.deg$avg_log2FC > 1 , 'Up', 'Down'))
  table(s.deg$threshold)
  dim(s.deg)
  ###标记出显著和不显著并转化为因子
  s.deg$adj_p_signi <- as.factor(ifelse(s.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
  s.deg$thr_signi <- paste0(s.deg$threshold, "_", s.deg$adj_p_signi)
  
  saveRDS(s.deg,"s.deg_filter_markers_patient_VS_carrier.rds")
  
  write.csv(s.deg,"carrier_patient_组间差异分析/s_deg_filter_patient_VS_carrier.csv",row.names = T)
  
  
  #====读取数据，绘图=========================================================
  r.deg <- readRDS("r.deg_filter_markers_patient_VS_carrier.rds")
  s.deg <- readRDS("s.deg_filter_markers_patient_VS_carrier.rds")
  
  ###差异基因绘图
  ### 自定义显示想要展示的基因名，这里挑选log2FC为top5的基因进行展示
  ###每个celltype找前5个上调基因
  top_up_label <- s.deg %>% 
    subset(., threshold%in%"Up") %>% 
    group_by(celltype) %>% 
    top_n(n = 5, wt = avg_log2FC) %>% 
    as.data.frame()
  ###每个celltype找前5个下调基因
  top_down_label <- s.deg %>% 
    subset(., threshold %in% "Down") %>% 
    group_by(celltype) %>% 
    top_n(n = -5, wt = avg_log2FC) %>% 
    as.data.frame()
  ###将上调和下调进行合并
  top_label <- rbind(top_up_label,top_down_label)
  
  library(scRNAtoolVis)
  library(jjPlot)
  library(ggrepel)
  
  colors <- c("red", "blue", "green", "yellow", "purple", "orange", "pink", "cyan", "brown", "black")
  
  ### 画图展示jjVolcano
  r.deg$cluster <- r.deg$celltype #jjVolcano函数只能默认识别cluster列，故复制新增一列存放细胞类型
  jjVolcano(diffData =r.deg, 
            tile.col = colors[1:11],
            pSize = 0.4,###调节点的大小
            legend.position=c(0.1,0.9),###调节图注坐标位置
            celltypeSize=2,###调节细胞celltype文字大小
            topGeneN=5)+
    labs(title = "patient VS carrier")  # 添加标题
  ggsave(
    filename = "carrier_patient_组间差异分析/jjVolcano_patient_VS_carrier.pdf",
    plot = last_plot(),  # 关键：指定要保存的图（last_plot()获取最后生成的图）
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300  # 可选：增加分辨率，避免图模糊
  )
  
  #scRNAtoolVis::markerVolcano()
  
  # 绘制火山图2            
  markerVolcano(            
    markers = r.deg,            
    topn = 5,        
    labelCol = ggsci::pal_npg()(11)
  )
  ggsave(
    filename = "carrier_patient_组间差异分析/Volcano_patient_VS_carrier.pdf",
    plot = last_plot(),  # 关键：指定要保存的图（last_plot()获取最后生成的图）
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300  # 可选：增加分辨率，避免图模糊
  )
  
  
  sub_s.deg <- subset(s.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
  ##绘制韦恩图
  library(VennDiagram)
  library(ggvenn)
  library(ggplot2)  # 依赖ggplot2
  
  # T cells 基因集：严格筛选（p_val_adj<0.05，|log2FC|>0.5）
  # 筛选特定细胞类型的行（逻辑索引：TRUE表示符合条件）
  filter_rows <- sub_s.deg$celltype == "T cell"
  # 提取基因列+去重+去NA
  degs_list1 <- unique(sub_s.deg$gene[filter_rows])
  degs_list1 <- degs_list1[!is.na(degs_list1)]
  # 查看结果
  print(degs_list1)
  
  # B cells 基因集：严格筛选（p_val_adj<0.05，|log2FC|>0.5）
  # 筛选特定细胞类型的行（逻辑索引：TRUE表示符合条件）
  filter_rows <- sub_s.deg$celltype == "B cell"
  # 提取基因列+去重+去NA
  degs_list2 <- unique(sub_s.deg$gene[filter_rows])
  degs_list2 <- degs_list2[!is.na(degs_list2)]
  # 查看结果
  print(degs_list2)
  
  # NK cells 基因集：严格筛选（p_val_adj<0.05，|log2FC|>0.5）
  # 筛选特定细胞类型的行（逻辑索引：TRUE表示符合条件）
  filter_rows <- sub_s.deg$celltype == "NK cell"
  # 提取基因列+去重+去NA
  degs_list3 <- unique(sub_s.deg$gene[filter_rows])
  degs_list3 <- degs_list3[!is.na(degs_list3)]
  # 查看结果
  print(degs_list3)
  
  # monocytes cells 基因集：严格筛选（p_val_adj<0.05，|log2FC|>0.5）
  # 筛选特定细胞类型的行（逻辑索引：TRUE表示符合条件）
  filter_rows <- sub_s.deg$celltype == "Monocyte"
  # 提取基因列+去重+去NA
  degs_list4 <- unique(sub_s.deg$gene[filter_rows])
  degs_list4 <- degs_list4[!is.na(degs_list4)]
  # 查看结果
  print(degs_list4)
  
  
  # 2. 整理为基因集列表（命名方便后续展示）
  gene_sets <- list(
    T_cell = degs_list1,  # 命名1
    B_cell = degs_list2,   # 命名2
    NK_cell = degs_list3,    # 命名3
    Monocytes = degs_list4  # 命名4
  )
  
  
  # 2. 提取你的4个基因集（无需修改，直接对应gene_sets）
  g1 <- gene_sets[["T_cell"]]
  g2 <- gene_sets[["B_cell"]]
  g3 <- gene_sets[["NK_cell"]]
  g4 <- gene_sets[["Monocytes"]]
  
  # 3. 自动计算4个基因集的所有交集区域数量（无需手动输入，避免出错）
  area1 <- length(g1)
  area2 <- length(g2)
  area3 <- length(g3)
  area4 <- length(g4)
  
  n12 <- length(intersect(g1, g2))
  n13 <- length(intersect(g1, g3))
  n14 <- length(intersect(g1, g4))
  n23 <- length(intersect(g2, g3))
  n24 <- length(intersect(g2, g4))
  n34 <- length(intersect(g3, g4))
  
  n123 <- length(intersect(intersect(g1, g2), g3))
  n124 <- length(intersect(intersect(g1, g2), g4))
  n134 <- length(intersect(intersect(g1, g3), g4))
  n234 <- length(intersect(intersect(g2, g3), g4))
  
  n1234 <- length(intersect(intersect(intersect(g1, g2), g3), g4))
  
  # 4. 绘制4细胞类型基因集韦恩图（标签自动使用你的基因集名称）
  venn_plot <- draw.quad.venn(
    area1 = area1,
    area2 = area2,
    area3 = area3,
    area4 = area4,
    n12 = n12,
    n13 = n13,
    n14 = n14,
    n23 = n23,
    n24 = n24,
    n34 = n34,
    n123 = n123,
    n124 = n124,
    n134 = n134,
    n234 = n234,
    n1234 = n1234,
    category = names(gene_sets),  # 自动显示T_cell/B_cell/NK_cell/monocytes_cell
    fill = c("#FF6B6B", "#4ECDC4", "#7927D1", "#4657D1"),  # 对应4个细胞类型颜色
    alpha = 0.6,                  # 透明度
    lwd = 1,                      # 边框线宽
    fontsize = 5,                 # 交集数量文字大小
    cat.fontsize = 6,             # 细胞类型标签文字大小
    cat.fontface = "bold"         # 标签加粗，更清晰
  )
  
  # 5. 显示韦恩图
  grid.draw(venn_plot)
  
  # 6. 保存高清图片（两种格式可选，满足论文投稿要求）
  # TIFF格式（无损高清）
  tiff("carrier_patient_组间差异分析/patient_VS_carrier_4_cell_types_venn_VennDiagram.tiff", width = 11, height = 10, units = "in", res = 300)
  grid.draw(venn_plot)
  dev.off()
  
  
  
}





##### 15.4、根据前述组间差异基因进行网站富集分析后的结果绘图（可加载数据）Metascape富集全套绘图（5类标准图）####
{
  # 首次运行安装包
  #install.packages(c("tidyverse","ggplot2","VennDiagram","grid","readxl"))
  library(tidyverse)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
  library(readxl)
  

  # 读取三组Metascape富集结果
  df_c <- read_excel("metascape_result_NK_carrier&control.xlsx", sheet = "Enrichment")
  df_c$group = "carrier_vs_control"
  
  df_p <- read_excel("metascape_result_NK_patient&control.xlsx", sheet = "Enrichment")
  df_p$group = "patient_vs_control"
  
  df_pvsc <- read_excel("metascape_result_NK_patient&carrier.xlsx", sheet = "Enrichment")
  df_pvsc$group = "patient_vs_carrier"
  
  all_enrich <- bind_rows(df_c, df_p, df_pvsc)
  
  # 清洗过滤：剔除Summary条目
  all_enrich_clean <- all_enrich %>%
    filter(!str_detect(GroupID, fixed("Summary")))
  
  keep_cat = c("Hallmark Gene Sets","KEGG Pathway","Reactome Gene Sets","GO Biological Processes")
  all_enrich_filter = all_enrich_clean %>%
    filter(Category %in% keep_cat) %>%
    mutate(
      neg_logq = -`Log(q-value)`,
      p.adjust = 10^(`Log(q-value)`),
      InTerm_InList = str_remove_all(InTerm_InList, " "),
      hit_n = as.numeric(str_extract(InTerm_InList, "^\\d+")),
      total_n = as.numeric(str_extract(InTerm_InList, "(?<=/)\\d+$")),
      GeneRatio = hit_n / total_n,
      Count = hit_n,
      # 简化通路名称
      short_term = case_when(
        str_starts(Term, "M") ~ str_remove(Description, "HALLMARK "),
        str_starts(Term, "hsa") ~ Description,
        str_starts(Term, "R-HSA") ~ Description,
        str_starts(Term, "GO") ~ Description
      ),
      # 数据库分组
      db = case_when(
        Category == "Hallmark Gene Sets" ~ "Hallmark",
        Category == "KEGG Pathway" ~ "KEGG",
        Category == "Reactome Gene Sets" ~ "Reactome",
        Category == "GO Biological Processes" ~ "GO_BP"
      ),
      db = factor(db, levels = c("Hallmark","KEGG","Reactome","GO_BP")),
      # 固定分组顺序，匹配配色
      group = factor(group, levels = c("carrier_vs_control","patient_vs_control","patient_vs_carrier"))
    ) %>%
    filter(`Log(q-value)` < -1.3) # q<0.05筛选显著通路
  
  # 统一配色 + 标签映射（绑定名称，彻底杜绝顺序错乱）
  color_list <- c(
    "carrier_vs_control" = "#2E86AB",
    "patient_vs_control" = "#A23B72",
    "patient_vs_carrier" = "#FF9533"
  )
  label_map <- c(
    "carrier_vs_control" = "carrier vs ctrl",
    "patient_vs_control" = "patient vs ctrl",
    "patient_vs_carrier" = "patient vs carrier"
  )
  
  # == 1、三组气泡图（横坐标组别）（5类图第1种） ===
  top_enrich = all_enrich_filter %>%
    group_by(group, db) %>%
    arrange(desc(neg_logq)) %>%
    slice_head(n = 12) %>%  #前12条通路，排序依据：desc(neg_logq)，-log10 (q-value) 降序，neg_logq = -log10 (校正后 q 值)，数值越大代表通路富集显著性越高
    ungroup()
  
  p_bubble = ggplot(top_enrich, aes(x = group, y = short_term)) +
    geom_point(aes(size = Count, fill = neg_logq), shape = 21, stroke = 0.3) +
    scale_fill_viridis_c(name = "-log10(q-value)") +
    scale_size(range = c(2, 8), name = "Gene count") +
    facet_wrap(~db, scales = "free_y", ncol = 2) +
    labs(x = "Comparison Group", y = "Pathway", title = "NK cell Enrichment Bubble Plot (Three Groups)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 7.5),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.background = element_rect(fill = "#336699"),
      strip.text = element_text(color = "white", face = "bold")
    )
  ggsave("NK_threeGroup_bubble_compare.pdf", p_bubble, width = 26, height = 14, device = "pdf", dpi = 300)
  print(p_bubble)
  
  # ==== 2、三组Top通路条形图（5类图第2种） ===
  top_all_sep = all_enrich_filter %>%
    group_by(group, db) %>%
    arrange(desc(neg_logq)) %>%
    slice_head(n = 12) %>%  #前12条通路，排序依据：desc(neg_logq)，-log10 (q-value) 降序，neg_logq = -log10 (校正后 q 值)，数值越大代表通路富集显著性越高
    ungroup()
  
  p_bar = ggplot(top_all_sep, aes(x = neg_logq, y = short_term, fill = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    facet_wrap(~db, scales = "free_y", ncol = 2) +
    scale_fill_manual(
      values = color_list,
      limits = names(color_list),
      labels = label_map
    ) +
    labs(x = "-log10(q-value)", y = "Pathway", title = "Top Significant Pathways Three Groups") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 7))
  ggsave("NK_threeGroup_bar_top.pdf", p_bar, width = 26, height = 12, device = "pdf", dpi = 300)
  print(p_bar)
  
  # ==== 3、三组三元韦恩图（5类图第3种）【优化：复用全局配色】 ===
  path_c <- unique(all_enrich_filter[all_enrich_filter$group=="carrier_vs_control",]$short_term)
  path_p <- unique(all_enrich_filter[all_enrich_filter$group=="patient_vs_control",]$short_term)
  path_pvsc <- unique(all_enrich_filter[all_enrich_filter$group=="patient_vs_carrier",]$short_term)
  
  venn3 <- draw.triple.venn(
    area1 = length(path_c),
    area2 = length(path_p),
    area3 = length(path_pvsc),
    n12 = length(intersect(path_c, path_p)),
    n13 = length(intersect(path_c, path_pvsc)),
    n23 = length(intersect(path_p, path_pvsc)),
    n123 = length(intersect(intersect(path_c, path_p), path_pvsc)),
    category = c("carrier_vs_ctrl","patient_vs_ctrl","patient_vs_carrier"),
    fill = c("#2E86AB","#A23B72","#FF9533"),
    alpha = 0.6, cat.cex = 1.1, cex = 1.5,
    cat.dist = 0.05, margin = 0.15
  )
  grid.draw(venn3)
  ggsave("NK_threeGroup_venn.pdf", plot = venn3, width = 8, height = 7, device = "pdf", dpi = 300)
  
  # === 4、各数据库通路数量柱状（5类图第4种） ===
  stat_df = all_enrich_filter %>%
    group_by(group, db) %>%
    tally(name = "count") %>%
    ungroup() %>%
    mutate(group = factor(group, levels = c("carrier_vs_control","patient_vs_control","patient_vs_carrier")))
  
  p_stat = ggplot(stat_df, aes(x = db, y = count, fill = group)) +
    geom_col(position = position_dodge(0.65), width = 0.65) +
    scale_fill_manual(
      values = color_list,
      labels = label_map,
      limits = names(color_list)
    ) +
    labs(x = "Database", y = "Number of significant pathways", title = "Pathway Count Comparison (Three Groups)") +
    theme_bw()
  ggsave("NK_threeGroup_count_bar.pdf", p_stat, width = 9, height = 6, device = "pdf", dpi = 300)
  print(p_stat)
  
  # === 5、三组分面标准气泡总图（横坐标基因比率）（5类图第5种） ===
  two_group_data <- all_enrich_filter %>%
    group_by(group, db) %>%
    arrange(desc(neg_logq)) %>%
    slice_head(n = 8) %>%  #前8条通路，排序依据：desc(neg_logq)，-log10 (q-value) 降序，neg_logq = -log10 (校正后 q 值)，数值越大代表通路富集显著性越高
    ungroup()
  
  p_two_facet <- ggplot(two_group_data, aes(x = GeneRatio, y = short_term)) +
    geom_point(aes(size = Count, fill = p.adjust), shape = 21, stroke = 0.3, color = "gray30") +
    # 【可选】交换高低颜色贴合文献习惯，按需启用
    # scale_fill_gradient(low = "#4178d8", high = "#de3838", trans = "log10", name = "p.adjust") +
    scale_fill_gradient(low = "#de3838", high = "#4178d8", trans = "log10", name = "p.adjust") +
    scale_size(range = c(2, 8), name = "Gene Count") +
    facet_grid(group ~ db, scales = "free_y") +
    labs(x = "GeneRatio", y = "") +
    theme_bw()
  ggsave("NK_threeGroup_all_facet_bubble.pdf", p_two_facet, width = 22, height = 10, device = "pdf", dpi = 300)
  print(p_two_facet)
  
  
  
}


