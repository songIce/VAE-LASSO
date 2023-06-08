library(GSVA)
library(GSEABase)
library(clusterProfiler)

Type_Group = "ECBS"
# FPKM表达量矩阵
geneExpr = fread("/home/bioinfo/FPKM_cli/fpkm/TCGA-LIHC_htseq_fpkm.tsv") %>% as.data.frame(.)
# ECBS分组模型
ECBS_Group = readRDS("/home/bioinfo/TCGA/pancancer/4.model_data/results/LIHC/ECBS_Model.rds")  
# #试一下EBG
# EBG_Group = readRDS("/home/bioinfo/TCGA/pancancer/data/EC_Model.rds") 
# ECBS_Group = readRDS("/home/bioinfo/TCGA/pancancer/data/EC_Model.rds") 
# rownames(ECBS_Group) = ECBS_Group$sample
# 基因名转换
eg <- bitr(geneExpr$gene_id, fromType = "ENSEMBL", toType = c("ENTREZID", "ENSEMBL",
                                                           "SYMBOL"), OrgDb = "org.Hs.eg.db")
# 基因表达量sample筛选
colnames(geneExpr) <- colnames(geneExpr) %>% gsub("[.]","-",.) 
geneExpr = geneExpr[,c(1,match(rownames(ECBS_Group),colnames(geneExpr)))]

mergedata <- merge(eg, geneExpr, by.x = "ENSEMBL", by.y = "gene_id")
mergedata <- mergedata[!duplicated(mergedata$SYMBOL), ]
expMatrix <- mergedata[, 4:ncol(mergedata)]
rownames(expMatrix) = mergedata[, 2]
#载入gmt文件
geneSets <- getGmt("/home/bioinfo/TCGA/pancancer/GSVA/pathway.gmt")  # 全部基因集
# KEGG
geneSets <- getGmt("/home/bioinfo/TCGA/pancancer/GSVA/msigdb_v7.0_files_to_download_locally/msigdb_v7.0_GMTs/c2.cp.kegg.v7.0.entrez.gmt") 
#geneSets <- getGmt("/home/bioinfo/TCGA/pancancer/GSVA/msigdb_v7.0_files_to_download_locally/msigdb_v7.0_GMTs/h.all.v7.0.entrez.gmt") 
GSVA_hall <- gsva(expr=as.matrix(expMatrix), 
                  gset.idx.list=geneSets, 
                  #   method="gsva", #c("gsva", "ssgsea", "zscore", "plage")
                  mx.diff=T, # 数据为正态分布则T，双峰则F
                  kcdf="Gaussian", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=20,# 并行线程数目
                  min.sz=2) 

library(limma)
# 设置或导入分组
ECBS_Group$group[ECBS_Group$group=="risk high"] = "high"
ECBS_Group$group[ECBS_Group$group=="risk low"] = "low"

design <- model.matrix(~0 + ECBS_Group$group)
colnames(design) = levels(factor(ECBS_Group$group))
rownames(design) = colnames(GSVA_hall)

compare <- makeContrasts(high - low, levels = design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef = 1, number = 1685)

FC_top_num = 40  #挑选前40个
pfiltGeneSets <- topTable(fit3, coef=1,p.value=0.05, number=Inf,adjust="fdr")
FC_threshold=ifelse(NROW(pfiltGeneSets) < FC_top_num,
                    sort(abs(pfiltGeneSets$logFC),decreasing = T)[NROW(pfiltGeneSets)]-1e-9,
                    sort(abs(pfiltGeneSets$logFC),decreasing = T)[FC_top_num]-1e-9)
print(paste0("FC threshold: ",FC_threshold))
DEgeneSets <- topTable(fit3, coef=1, number=Inf,
                       p.value=0.05,
                       lfc=FC_threshold,
                       adjust="fdr")

#################################################################################################
# 画条形图
set.seed(1234)
Diff_sample<-sample_n(Diff, 40, replace = FALSE)
Diff_sample = DEgeneSets
dat_plot <- data.frame(id = row.names(Diff_sample),
                       t = Diff_sample$t)
# 去掉"HALLMARK_"
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","") %>% 
  str_replace(. , "REACTOME_","") %>% 
  str_replace(. , "KEGG_","") %>% 
  gsub("_"," ",.)
# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t  >-2, ifelse(dat_plot$t >= 2 ,'Up','Stable'),'Down'),levels=c('Up','Down','Stable'))
# 排序
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)
library(ggthemes)
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#DC143C','NoSignifi'='#808080','Down'='#00008B')) +
  geom_hline(yintercept = c(-2,2),
             color = 'white',
             size = 0.5,
             lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, Risk high-VS-Risk Low') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
# 添加标签
# 小于-2的数量
low1 <- dat_plot %>% filter(t < -2) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于2总数量
high0 <- dat_plot %>% filter(t < 2) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)
# 依次从下到上添加标签
p = p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
              hjust = 0,color = 'black',size = 2.5) + # 小于-1的为黑色标签
  # geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
  #           hjust = 0,color = 'grey',size = 2.5) + # 灰色标签
  # geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
  #           hjust = 1,color = 'grey',size = 2.5) + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black',size = 2.5) # 大于1的为黑色标签
p
ggsave(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSVA/",
             Type_Group,"_gsva_bar.pdf",sep = ""),p,width = 10,height  = 6)
#################################################################################################





# 随机抽取40个
library(dplyr)
# GSVA_Sample <- sample_n(as.data.frame(GSVA_hall), 40, replace = FALSE)

# 更换热图和火山图的id，免得过长，详细信息可以附在supplementary materials中
T_ALLres = Diff
T_DEres_raw = DEgeneSets
T_DEres= T_DEres_raw %>%
  rownames_to_column(var="terms") %>%
  mutate(term_simp=gsub("^GO_","",terms) %>% gsub("^KEGG_","",.) %>% gsub("_PROCESS$","",.) %>% gsub("_PATHWAY$","",.)) %>%
  mutate(term_simp=tolower(term_simp)) %>%
  mutate(code=1:NROW(T_DEres_raw)) %>%
  mutate(term_simp_code=paste0(term_simp," (",code,") ")) %>%
  column_to_rownames(var="terms")

for (i in rownames(T_DEres)) {
  rownames(T_ALLres)[match(i,rownames(T_ALLres))]=T_DEres[match(i,rownames(T_DEres)),"code"]
}

order_ind = ECBS_Group %>% rownames_to_column(var="ind") %>% arrange(riskScore) %>% .$ind
# 画热图
heat_data=GSVA_hall %>% as.data.frame() %>% 
  .[rownames(T_DEres_raw),order_ind]%>% 
  rownames_to_column(var="terms") %>% 
  mutate(term_simp=gsub("^GO_","",terms) %>% 
           gsub("^KEGG_","",.) %>% 
           gsub("_PROCESS$","",.) %>% 
           gsub("_PATHWAY$","",.) %>% 
           gsub("REACTOME_","",.)) %>% 
  mutate(term_simp=tolower(term_simp)) %>% 
  mutate(code=1:NROW(T_DEres_raw)) %>% 
  mutate(term_simp_code=paste0(term_simp," (",code,") ")) %>% 
  mutate(term_simp_code=gsub("_"," ",term_simp_code)) %>% 
  dplyr::select(-c(term_simp,code,terms)) %>% 
  column_to_rownames(var="term_simp_code")

t1=pheatmap::pheatmap(heat_data,
                      # color = colorRampPalette(c("steelblue", "white","tomato"))(40),
                      cluster_rows = T, 
                      cluster_cols = F,
                      clustering_method="complete",
                      fontsize_row = 8,height = 11,
                      annotation_col = ECBS_Group[,"group",drop=F],
                      show_colnames = F)
T1 = ggplotify::as.ggplot(t1)
ggsave(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSVA/",Type_Group,
             "_gsva.heatmap.pdf",sep = ""),T1,width=12,height = 5)

# 画火山图
# rownames(Dat)= rownames(Dat) %>% gsub("^REACTOME_","",.) %>% tolower() %>% gsub("_"," ",.)
draw_volcano=function(Dat,adjPvalueCutoff,logFCcutoff,myxlim){
  #加载包
  library(ggplot2)
  library(ggrepel)
  library(ggThemeAssist)
  #读取数据
  #确定是上调还是下调，用于给图中点上色）
  Dat$threshold = factor(ifelse(Dat$adj.P.Val < adjPvalueCutoff & abs(Dat$logFC) >= logFCcutoff, 
                                ifelse(Dat$logFC>= logFCcutoff ,'Up','Down'),
                                'NoSignifi'),
                         levels=c('Up','Down','NoSignifi'))

  Dat2=Dat %>% rownames_to_column(var="id")
  Dat2$id = Dat2$id %>% gsub("KEGG_","",.) %>% gsub("REACTOME_","",.) %>% gsub("_"," ",.)
  
  p2=ggplot(Dat2,aes(x=logFC,y=-log10(adj.P.Val),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("#DC143C","#00008B","#808080"))+#确定点的颜色
    geom_text_repel(
      data = Dat2[Dat2$adj.P.Val<adjPvalueCutoff&abs(Dat2$logFC)>logFCcutoff,],
      aes(label = id),
      size = 3.5,
      segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
    theme_bw()+#修改图片背景
    theme(
      legend.title = element_blank(),
      panel.grid.major =element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank()#不显示图例标题
    )+
    ylab('-log10 (p-adj)')+#修改y轴名称
    xlab('log2 (FoldChange)')+#修改x轴名称
    geom_vline(xintercept=c(logFCcutoff*-1,logFCcutoff),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
    geom_hline(yintercept = -log10(adjPvalueCutoff),lty=3,col="black",lwd=0.5) + #添加竖线padj<0.05
    xlim(myxlim*-1,myxlim)
  return(p2)
}

#ECBS
logFCcutoff=0.6121460  #all geneset
logFCcutoff=0.3111918  #kegg
#EBG
logFCcutoff=0.4412650  #all geneset
logFCcutoff=0.2481833  #kegg

adjPvalueCutoff <- 0.05
T2=draw_volcano(T_ALLres,adjPvalueCutoff,logFCcutoff,0.71)
T2
ggsave(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSVA/",Type_Group,
             "_gsva.hotplot.pdf",sep = ""),T2,width=10,height = 5)


write.csv(T_DEres, paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSVA/",
                         Type_Group,"_gsva_hotplotData.csv",sep = ""))

