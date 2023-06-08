library("data.table")
library("dplyr")
library("TCGAutils")

# 脚本画图内容概要
# 1.癌和癌旁样本中的AEI值做一组对比的小提琴图
# 2.HCC功能性RNA编辑位点的鉴定部分



## 1.把癌和癌旁样本中的AEI值做一组对比的小提琴图，如下：
## 可以做两组，一组是癌和对应癌旁样本的，各50个；
## 另一组是所有癌和50个癌旁的。配色跟上面的图保持一致。两组图分两张。
## 还要对比一下两组中位数之间的差异显著性，用p值来标上面的※

# 先获取50对sample和所有的tumor
cancerData <- fread("/home/bioinfo/TCGA/pancancer/informative_RES/LIHC_EF_informative_no_SNP_2.csv")
colnames(exprADAR) <- gsub("[.]","-",colnames(exprADAR))
allSamples = colnames(cancerData)[-1:-26]
tumorSampleAll = c()
normalSample = c()
for(oneSamle in allSamples) {
  name_list <- oneSamle %>% substring(14,15) %>% as.numeric()
  if(name_list>=1&name_list<=9) {
    tumorSampleAll = c(tumorSampleAll,oneSamle)
  }else {
    normalSample = c(normalSample,oneSamle)
  }
}
match(substr(normalSample,1,12),substr(tumorSampleAll,1,12)) %>% tumorSampleAll[.] -> tumorSample

fread("/home/bioinfo/TCGA/pancancer/AEI/LIHC/LIHC_EditingIndex.csv") -> AEI_LIHC_Data
AEI_LIHC_Data[AEI_LIHC_Data$StrandDecidingMethod=="RefSeqThenMMSites"] -> AEI_LIHC_Data
AEI_LIHC_Data$SamplePath  %>% strsplit(.,"[/]") -> sampleName
allSamples = c()
for (oS in sampleName) {
  allSamples = c(allSamples,oS[8])
}
res <- TCGAutils::filenameToBarcode(allSamples)   ## 这个包需要联网，暂时做不了

# 先画50对tumor和normal的
boxplotAEI = data.frame()
groupTumorSample = res$file_name[match(tumorSample, res$aliquots.submitter_id)] %>% 
  match(.,allSamples) %>% AEI_LIHC_Data$A2GEditingIndex[.]
boxplotAEI = rbind(boxplotAEI,cbind( `AEI value` =  groupTumorSample, "group" = rep("Tumor",length(groupTumorSample))))
groupNormalSample = res$file_name[match(normalSample, res$aliquots.submitter_id)] %>% 
  match(.,allSamples) %>% AEI_LIHC_Data$A2GEditingIndex[.]
boxplotAEI = rbind( boxplotAEI, cbind( `AEI value` =  groupNormalSample, "group" = rep("Peritumor",length(groupNormalSample))))
boxplotAEI$`AEI value` <- boxplotAEI$`AEI value` %>% as.character() %>% as.numeric()
boxplotAEI$group <- boxplotAEI$group %>% as.character() %>% as.factor()
compaired <- list(c("Tumor", "Peritumor"))
p1 <- ggviolin(boxplotAEI, x="group", y="AEI value", 
              palette= "hue" , add = "boxplot", fill = "group",
              add.params = list(fill = "group"), xlab = "", ylab = "AEI (%)",
              ggtheme = theme_pubr(base_size = 15))
p1 <- p1 + stat_compare_means(comparisons = compaired,paired=T)
ggsave(paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,"/3.1.2.A.AEI_BoxPlot.pdf",sep = ''),p1)
wilcox.test(as.numeric(boxplotAEI$`AEI value`[boxplotAEI$group=="Tumor"]), 
            as.numeric(boxplotAEI$`AEI value`[boxplotAEI$group=="Peritumor"]),  paired = T)



# 再画所有tumor和normal的
boxplotAEI = data.frame()
groupTumorSample = res$file_name[match(tumorSampleAll, res$aliquots.submitter_id)] %>% 
  match(.,allSamples) %>% AEI_LIHC_Data$A2GEditingIndex[.]
boxplotAEI = rbind(boxplotAEI,cbind( `AEI value` =  groupTumorSample, "group" = rep("Tumor",length(groupTumorSample))))
boxplotAEI = rbind( boxplotAEI, cbind( `AEI value` =  groupNormalSample, "group" = rep("Peritumor",length(groupNormalSample))))
boxplotAEI$`AEI value` <- boxplotAEI$`AEI value` %>% as.character() %>% as.numeric()
boxplotAEI$group <- boxplotAEI$group %>% as.character() %>% as.factor()
compaired <- list(c("Tumor", "Peritumor"))
p2 <- ggviolin(boxplotAEI, x="group", y="AEI value", 
               palette= "hue" , add = "boxplot", fill = "group",
               add.params = list(fill = "group"), xlab = "", ylab = "AEI (%)",
               ggtheme = theme_pubr(base_size = 15))
p2 <- p2 + stat_compare_means(comparisons = compaired)
ggsave(paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,"/3.1.2.B.AEI_BoxPlot.pdf",sep = ''),p2)
# 手动算p值，自动算的不好改
wilcox.test(as.numeric(boxplotAEI$`AEI value`[boxplotAEI$group=="Tumor"]), 
            as.numeric(boxplotAEI$`AEI value`[boxplotAEI$group=="Peritumor"]),  paired = T)


# p = ggarrange(p1, p2, ncol = 2, nrow = 1)
# fileName = paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,"/3.1.2.C.AEI_BoxPlot.pdf",sep = '')
# ggsave(fileName, p, width=8, height = 6)


# colnames(boxplotAEI) <- c("value","group")
# boxplotAEI$value <- boxplotAEI$value %>% as.character() %>% as.numeric()
# boxplotAEI$group <- boxplotAEI$group %>% as.character() %>% as.factor()
# compaired <- list(c("A", "B"),c("B","C"),c("A","C"))
# p <- ggboxplot(boxplotAEI,x="group",y="value",color = "group",palette=c("#E41A1C","#4DAF4A","#1E90FF"))
# p <- p + stat_compare_means(comparisons = compaired)
# ggsave(paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,"/1.2.B.AEI_BoxPlot.pdf",sep = ''),p)





##########################################################################################################
# 2.HCC功能性RNA编辑位点的鉴定部分
# 1）editing level 和所在基因的表达量差异分析图（7个位点以及它们所在的5个基因）
# 的纵坐标不对，左图应该是”RNA editing level”；右图是“Gene expression level log2FPKM”.横坐标都写成“tumor”和“”，所有图中不再出现normal。

cancerName = "LIHC"
#先处理一下得到原始矩阵
cancerData = fread(paste("/home/bioinfo/TCGA/pancancer/informative_RES/",cancerName,"_EF_informative_no_SNP_2.csv",sep = ''))
#先导入stage和survival和peptide的数据
ostimeData = fread(paste('/home/bioinfo/TCGA/pancancer/CancerOstimepValue/','all_ostime_pvalue.csv',sep = ''))
peptideData = fread("/home/bioinfo/TCGA/pancancer/CancerPeptide/positionAndCancer.csv")
ostimeData = ostimeData[ostimeData$cancer==cancerName]
peptideData = peptideData[peptideData$cancer==cancerName]

paste(ostimeData$Chr,ostimeData$`Position(1base)`,sep = '_') -> ostimePosition
peptideData$position -> peptidePosition
match(peptidePosition,ostimePosition) %>% na.omit %>% ostimePosition[.]-> matchPosit
paste(cancerData$Chr,cancerData$`Position(1base)`,sep = '_') -> cancerPosit
match(matchPosit,cancerPosit) -> matchRow
cancerData$Gene.wgEncodeGencodeBasicV34[matchRow]

exprADAR = fread(paste("/home/bioinfo/FPKM_cli/fpkm/TCGA-",cancerName,"_htseq_fpkm.tsv",sep = ''))
allSamples = colnames(cancerData)[-1:-26]
tumorSample = c()
normalSample = c()
for(oneSamle in allSamples) {
  name_list <- oneSamle %>% substring(14,15) %>% as.numeric()
  if(name_list>=1&name_list<=9) {
    tumorSample = c(tumorSample,oneSamle)
  }else {
    normalSample = c(normalSample,oneSamle)
  }
}

match(substr(normalSample,1,12),substr(tumorSample,1,12)) %>% tumorSample[.] -> tumorSample
colNum = match(gsub("[.]","-",substr(tumorSample,1,28)),gsub("[.]","-",colnames(exprADAR))) %>% na.omit()
exprADARtumor = exprADAR[,..colNum]
colNum = match(gsub("[.]","-",substr(normalSample,1,28)),gsub("[.]","-",colnames(exprADAR))) %>% na.omit()
exprADARnormal = exprADAR[,..colNum]

geneName = c("ENSG00000188295","ENSG00000135451","ENSG00000188026","ENSG00000163453","ENSG00000129103")
exprADARtumorData = rbind( cbind(t(exprADARtumor[exprADAR$gene_id == "ENSG00000188295"]),rep("ZNF669",ncol(exprADARtumor))), 
                           cbind(t(exprADARtumor[exprADAR$gene_id == "ENSG00000135451"]),rep("TROAP",ncol(exprADARtumor))), 
                           cbind(t(exprADARtumor[exprADAR$gene_id == "ENSG00000188026"]),rep("RILPL1",ncol(exprADARtumor))), 
                           cbind(t(exprADARtumor[exprADAR$gene_id == "ENSG00000163453"]),rep("IGFBP7",ncol(exprADARtumor))), 
                           cbind(t(exprADARtumor[exprADAR$gene_id == "ENSG00000129103"]),rep("SUMF2",ncol(exprADARtumor)))) %>% as.data.frame(.)
exprADARtumorData = cbind(exprADARtumorData,rep("tumor",nrow(exprADARtumorData)))
colnames(exprADARtumorData) <- c("value","gene","group")
exprADARnormalData = rbind( cbind(t(exprADARnormal[exprADAR$gene_id == "ENSG00000188295"]),rep("ZNF669",ncol(exprADARnormal))), 
                            cbind(t(exprADARnormal[exprADAR$gene_id == "ENSG00000135451"]),rep("TROAP",ncol(exprADARnormal))), 
                            cbind(t(exprADARnormal[exprADAR$gene_id == "ENSG00000188026"]),rep("RILPL1",ncol(exprADARnormal))), 
                            cbind(t(exprADARnormal[exprADAR$gene_id == "ENSG00000163453"]),rep("IGFBP7",ncol(exprADARnormal))), 
                            cbind(t(exprADARnormal[exprADAR$gene_id == "ENSG00000129103"]),rep("SUMF2",ncol(exprADARnormal)))) %>% as.data.frame(.)
exprADARnormalData = cbind(exprADARnormalData,rep("peritumor",nrow(exprADARnormalData)))
colnames(exprADARnormalData) <- c("value","gene","group")

exprADARData = rbind(exprADARtumorData,exprADARnormalData)
exprADARData$value <- exprADARData$value %>% as.character(.) %>% as.numeric(.) 
exprADARData$value <- log2(exprADARData$value + 1)
exprADARData$gene <- as.factor(exprADARData$gene)
exprADARData$group <- as.factor(exprADARData$group)
remove_rownames(exprADARData) -> exprADARData

compare_means(value ~ group, data = exprADARData, 
              group.by = "gene", paired = TRUE)
# Box plot facetted by "dose"
p <- ggpaired(exprADARData, x = "group", y = "value",
              color = "group", palette = c("#EB6046","#397DB7"), 
              line.color = "gray", line.size = 0.3, 
              facet.by = "gene", short.panel.labs = T)
# Use only p.format as label. Remove method name.
p1 = p + stat_compare_means(label = "p.format", paired = TRUE, label.x = 1.2, label.y = 9.5) +
  labs(x="",y="Gene expression level log2FPKM") +
  theme(legend.position = "right") 
p1
ggsave(paste("/home/bioinfo/TCGA/pancancer/05.supplement/",cancerName,"hostGeneCompare.pdf",sep = ''),p1,width = 8,height = 7)


cancerData[matchRow,] -> matchCancerData
matchCol = match(tumorSample,colnames(matchCancerData))
tumorEditing = matchCancerData[,..matchCol] 
matchCol = match(normalSample,colnames(matchCancerData))
normalEditing = matchCancerData[,..matchCol]
tumorDataAll = data.frame()
geneName = cancerPosit[matchRow]
for (i in 1:nrow(normalEditing)) {
  tumorData = rbind(
    cbind(t(tumorEditing[i,]),rep("tumor",ncol(tumorEditing)),rep(geneName[i],ncol(tumorEditing))),
    cbind(t(normalEditing[i,]),rep("peritumor",ncol(normalEditing)),rep(geneName[i],ncol(normalEditing)))
  ) 
  tumorDataAll = rbind(tumorDataAll,tumorData)
}
colnames(tumorDataAll) <- c("value","group","position")
remove_rownames(tumorDataAll) -> tumorDataAll
tumorDataAll$position <- as.factor(tumorDataAll$position)
tumorDataAll$group <- as.factor(tumorDataAll$group)
tumorDataAll$value <- tumorDataAll$value %>% as.character(.) %>% as.numeric(.) 
compare_means(value ~ group, data = tumorDataAll, 
              group.by = "position", paired = TRUE)
# Box plot facetted by "dose"
p <- ggpaired(tumorDataAll, x = "group", y = "value",
              color = "group", palette = c("#397DB7","#EB6046"), 
              line.color = "gray", line.size = 0.4,
              facet.by = "position", short.panel.labs = T)
# Use only p.format as label. Remove method name.
p1 = p + stat_compare_means(label = "p.format", paired = TRUE, label.x = 1.3, label.y = 0.95) +
  labs(x="",y="RNA editing level") +
  theme(legend.position = "right") 
p1
ggsave(paste("/home/bioinfo/TCGA/pancancer/05.supplement/",cancerName,"editingLevelCompare.pdf",sep = ''),p1,height = 10,width = 10)


## 2）尝试增加一个热图
# ADAR表达量
exprADAR = fread("/home/bioinfo/FPKM_cli/fpkm/TCGA-LIHC_htseq_fpkm.tsv")
cancerData <- fread("/home/bioinfo/TCGA/pancancer/informative_RES/LIHC_EF_informative_no_SNP_2.csv")
colnames(exprADAR) <- gsub("[.]","-",colnames(exprADAR))
allSamples = colnames(cancerData)[-1:-26]
# 获取癌和癌旁的Sample
tumorSample = c()
normalSample = c()
for(oneSamle in allSamples) {
  name_list <- oneSamle %>% substring(14,15) %>% as.numeric()
  if(name_list>=1&name_list<=9) {
    tumorSample = c(tumorSample,oneSamle)
  }else {
    normalSample = c(normalSample,oneSamle)
  }
}
match(substr(normalSample,1,12),substr(tumorSample,1,12)) %>% tumorSample[.] -> tumorSample
# 获取癌和癌旁对应的ADAR1表达量
colNum = match(tumorSample,colnames(exprADAR)) %>% na.omit()
exprADARtumor = exprADAR[,..colNum]
colNum = match(normalSample,colnames(exprADAR)) %>% na.omit()
exprADARnormal = exprADAR[,..colNum]
exprADARtumorADAR1 = exprADARtumor[exprADAR$gene_id == "ENSG00000160710"]
exprADARnormalADAR1 = exprADARnormal[exprADAR$gene_id == "ENSG00000160710"]
# ADAR2,也就是ADARB1
exprADARtumorADAR2 = exprADARtumor[exprADAR$gene_id == "ENSG00000197381"]  
exprADARnormalADAR2 = exprADARnormal[exprADAR$gene_id == "ENSG00000197381"]
exprADARtumorADAR = rbind(exprADARtumorADAR1, exprADARtumorADAR2)
exprADARnormalADAR = rbind(exprADARnormalADAR1, exprADARnormalADAR2)

# 获取七个位点对应的编辑率
selectPosition = c("chr1_247100417","chr12_49324229","chr12_123518385","chr13_52643186","chr4_57110068","chr4_57110120","chr7_56078974")
matchCancerData = cancerData[match(selectPosition,paste(cancerData$Chr,cancerData$`Position(1base)`,sep = "_")),] %>% as.data.frame(.)
cancerDataTumor = matchCancerData[,c(1,2,match(tumorSample,colnames(matchCancerData)))]
cancerDataNormal = matchCancerData[,c(1,2,match(normalSample,colnames(matchCancerData)))]


# for (i in 1:5) {
#   wilcox.test(as.numeric(exprGeneTumor[i,]), as.numeric(exprGeneNormal[i,]),  paired = F) -> pvalues
#   cat(paste(pvalues$p.value,"\n"))
# }




# 获取五个基因的基因表达量
geneExpr = fread("/home/bioinfo/FPKM_cli/fpkm/TCGA-LIHC_htseq_fpkm.tsv")
colnames(geneExpr) <- gsub("[.]","-",colnames(geneExpr))
g2e=unique(toTable(org.Hs.egENSEMBL))
g2s=unique(toTable(org.Hs.egSYMBOL))
each_model_ori = c("ZNF669","TROAP","RILPL1","IGFBP7","SUMF2")
symbolIds = c()
for (i in 1:length(each_model_ori)) {
  symbolIds <- c(symbolIds, g2e$ensembl_id[g2e$gene_id == g2s$gene_id[g2s$symbol == each_model_ori[i]]])
}
exprGeneTumor = match(symbolIds,exprADAR$gene_id) %>% as.data.frame(exprADARtumor)[.,]
exprGeneNormal = match(symbolIds,exprADAR$gene_id) %>% as.data.frame(exprADARnormal)[.,]

# ADAR1数据标注
heatmapDataADAR = cbind(Tumor = exprADARtumorADAR, Normal = exprADARnormalADAR) %>% as.data.frame(.)
rownames(heatmapDataADAR) <- c("ARAR1","ADAR2")
annotation_col_ADAR <- data.frame( type = rep(c("Tumor", "peritumor"), each = 50) )
rownames(annotation_col_ADAR) <- colnames(heatmapDataADAR)
# 基因表达量数据标注
heatmapDataGeneExpr = cbind(Tumor = exprGeneTumor, Normal = exprGeneNormal)
rownames(heatmapDataGeneExpr) <- c("ZNF669","TROAP","RILPL1","IGFBP7","SUMF2")
# 位点对应的编辑率数据标注
heatmapDataEditRate = cbind(Tumor = cancerDataTumor[,c(-1:-2)], Normal = cancerDataNormal[,c(-1:-2)])
rownames(heatmapDataEditRate) <- c("chr1_247100417","chr12_49324229","chr12_123518385","chr13_52643186","chr4_57110068","chr4_57110120","chr7_56078974")

bk <- c(seq(-6,-2,by=0.1),seq(-1.9,-0.02,by=0.01),seq(0,1.9,by=0.01),seq(2,6,by=0.1))
pheatmap(log2(heatmapDataADAR), cluster_rows = F, show_colnames = F, cluster_cols = F, 
         annotation_col = annotation_col_ADAR, cellwidth = 5, cellheight = 20, scale = "row",
         colorRampPalette(colors = c("blue","white","red"))(50) ) -> pADAR1
pheatmap(log2(heatmapDataGeneExpr), cluster_rows = F, show_colnames = F, cluster_cols = F, 
        cellwidth = 5, cellheight = 20, scale = "row",
         colorRampPalette(colors = c("blue","white","red"))(50)) -> pGeneExpr
pheatmap(log2(heatmapDataEditRate+1), cluster_rows = F, show_colnames = F, cluster_cols = F, 
         cellwidth = 5, cellheight = 20, scale = "row",
         # color = c(colorRampPalette(colors = c("dodgerblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","tomato"))(length(bk)/2)),
         # breaks=bk,
         colorRampPalette(colors = c("dodgerblue","white","tomato"))(50)
         ) -> pEditRate


p = cowplot::plot_grid(pADAR1$gtable, pGeneExpr$gtable, pEditRate$gtable, ncol= 1, align = 'v', hjust=-10)
ggsave(paste("/home/bioinfo/TCGA/pancancer/05.supplement/heatmap.pdf",sep = ''),p,height = 7,width = 10)
ggdraw() +
  draw_plot(pADAR1$gtable, x = 0, y = 1, width = 1, height = .5) +
  draw_plot(pGeneExpr$gtable, x = 0, y = .5, width = 1, height = .5) +
  draw_plot(pEditRate$gtable, x = 0, y = 0, width = 1, height = .5) -> p2
ggsave(paste("/home/bioinfo/TCGA/pancancer/05.supplement/heatmap.pdf",sep = ''),p2,height = 5,width = 10)
