# 先取基因表达量数据的交集
geneExprLIRI = fread("/home/bioinfo/TCGA/pancancer/data/HCCDB18_mRNA_level3.txt")
geneExprGSE = fread("/home/bioinfo/TCGA/pancancer/data/GSE14520-GPL3921.gene.txt")
geneExprLIHC = fread("/home/bioinfo/TCGA/pancancer/LIHC_geneExpr_count.csv")

g2e=unique(toTable(org.Hs.egENSEMBL)) 
colnames(g2e) <- c("Entrez_ID","ensembl_id")
g2e$Entrez_ID = as.numeric(g2e$Entrez_ID)
selectGene = merge(geneExprGSE[,1], g2e, by="Entrez_ID")
selectGene = merge(selectGene, geneExprLIRI[,1], by="Entrez_ID")
geneLIHC = strsplit(geneExprLIHC$V1,"[.]")
allGene = c()
for (i in 1:length(geneLIHC)) {
  allGene = c(allGene,geneLIHC[[i]][1])
}
allGene = data.frame( ensembl_id = allGene )
selectGene = merge(selectGene, allGene, by="ensembl_id")
selectGene = unique(selectGene)
geneExprLIHC$ensembl_id = allGene
geneExprLIHCNew = merge(geneExprLIHC, unique(selectGene[,1]), by="ensembl_id")
write.csv(geneExprLIHCNew, "/home/bioinfo/TCGA/pancancer/LIHC_geneExpr_count_New.csv", row.names = F)







# 3.4.5 GES 预后模型的外部数据验证
# HCC的其他数据：选择：路径下的GSE14520表达量矩阵和LIRI 的表达量矩阵两组数据，以及对应的临床数据。
# 其他癌症数据：TCGA的其他癌症类型，包括NPC。
# 这两部分的额外数据都在iscience_data.zip
# 分别代入预后风险评估模型，
# a) 能否明确分成高、低两个风险组，以及两组间在预后、
# b) AUC、
# c) 单变量森林图、
# d)多因素cox回归森林图。HCC的全部都要，TCGA其他癌症的只把显著性的留下就可以。
#载入公共函数
{
  source("/home/bioinfo/TCGA/pancancer/Script/6.func.R")
  cancerNameComp = "LIRI"                      #给个名字
  # 先用lasso得到的基因和稀疏，建立新的模型
  setwd("/home/bioinfo/TCGA/pancancer/data")   #把目录设置到存放数据的位置
  geneExprCompFile = "HCCDB18_mRNA_level3.txt" #其他类型癌症的基因表达量
  immuneCompFile = "liri_immune.rds"           #其他类型癌症的免疫文件
  sampleCompFile = "HCCDB18.sample.txt"        #其他癌症的病人ID
  survCompFile = "HCCDB18.patient.txt"         #其他类型癌症的生存文件
  LIHCGeneFile = "/home/bioinfo/TCGA/pancancer/4.model_data/results/LIHC/3.4.1.A.Classifier.selectGeneData.csv"
  
  geneExprCompData = fread(geneExprCompFile)
  survCompData = fread(survCompFile) %>% t(.)
  colnames(survCompData) = survCompData[1,]
  survCompData = survCompData[-1,] %>% as.data.frame(.) %>% rownames_to_column(.)
  sampleCompData = fread(sampleCompFile)
  sampleCompData = t(sampleCompData) %>% as.data.frame(.) %>% rownames_to_column(.) %>% .[-1,]
  colnames(sampleCompData) <- c("sample", "TYPE", "SAMPLE_NAME", "SAMPLE_NAME1", "patient", "PATIENT")
  sampleCompData = sampleCompData[!sampleCompData$TYPE == "Adjacent",]   # 去掉癌旁
  # 载入LIHC得到的基因和系数
  selectGeneData = read.csv(LIHCGeneFile)
  colnames(selectGeneData) = c("gene", "coef", "Symbol")
  geneExprCompData = merge(selectGeneData,geneExprCompData,by="Symbol")
  # 用survminer包的函数surv_cutpoint来找到risk score的最佳cutoff值，
  # 不同cutoff时标准化的log-rank也会做一个统计量分布图，虚线为选择的最佳分割点;
  riskScoreComp = 0
  for (Ind in 5:ncol(geneExprCompData)) {
    riskScoreComp[Ind] = sum(geneExprCompData$coef * geneExprCompData[,Ind])
  }
  riskScoreComp = riskScoreComp[-1:-4]
  y_Comp = data.frame( sample = colnames(geneExprCompData)[-1:-4], riskScore = riskScoreComp ) 
  y_Comp = merge(y_Comp, sampleCompData[c("sample","patient")], by="sample")
  survCompDataPart = survCompData[c("rowname","SUR","STATUS")]
  colnames(survCompDataPart) <- c("patient", "time", "status")
  y_Comp = merge(y_Comp, survCompDataPart, by="patient")
  y_Comp$time = as.numeric(y_Comp$time)
  y_Comp$status = as.character(y_Comp$status)
  y_Comp$status[y_Comp$status == "Alive"] = 0
  y_Comp$status[y_Comp$status == "Dead"] = 1
  y_Comp$status = as.numeric(y_Comp$status)
  # 寻找cutpoint
  res.cut_Comp=survminer::surv_cutpoint(y_Comp, time = "time", event = "status", minprop = 0.004, variables=c("riskScore"), progressbar = F)
  cut_Comp= as.data.frame(summary(res.cut_Comp))["riskScore","cutpoint"]
  pdf("/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/3.4.1.d.comp.riskScore_cutPoint.pdf",height=6.5) 
  cut1_Comp = as.data.frame(summary(res.cut_Comp))["riskScore","cutpoint"]
  plot(res.cut_Comp$riskScore, pch=20, cex=1.4, xlab = "", ylab = "Standardized log-rank statistic") 
  text(x=cut1_Comp, y=res.cut_Comp$riskScore$stats %>% min(.) + 0.1, pos=4, labels = round(cut1_Comp,3))
  dev.off()
  y_Comp$group=ifelse(y_Comp$riskScore < cut1_Comp,"risk low","risk high")
  # e)预后模型评分在病人间的分布（每个病人的risk score是多少，其中高组和低组用不同的颜色表示）也要做一个图；
  library(ggsci); library(patchwork)
  plotRiskScore(y_Comp,"comp","LIRI")
  #############################
  #2）两个预后风险评估模型之间的对比分析
  # 画ROC曲线
  library(survivalROC); library(timeROC); library(timereg)
  colorName = c("#1E90FF","#FF4500")
  pdf("/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/3.4.2.A.comp.ROC.pdf") 
  each_yr_roc(y_Comp[c("time","status","riskScore","group")],colorName)
  dev.off()
  # KM生存分析
  survivalResultComp = survfit(Surv(time, status) ~ group, y_Comp)
  surv_By_RiskScore_Comp(survivalResultComp, "comp",cancerNameComp)
  cosResult = coxph(formula = Surv(time, status) ~ group, data = y_Comp)
  print(summary(cosResult))
  data.survdiff <- survdiff(Surv(time, status) ~ group,y_Comp)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
  cat(paste("HR:",HR,"up95:",up95,"low95",low95))
  KMdata = paste("HR:",HR,"up95:",up95,"low95",low95)
  write.table(KMdata, "/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/KMdata.txt",row.names = F)
  
  
  #热图不好画小提琴
  
  groupLabelData = merge(y_Comp,sampleCompData,by="sample")
  groupLabelData = groupLabelData[,c("group","PATIENT")]
  colnames(groupLabelData) <- c("group","rowname")

  immune_Data = readRDS(immuneCompFile)
  immune_Data = rownames_to_column(immune_Data)
  sigpara = c("tumor__Plasma.cells","tumor__T.cells.follicular.helper",
              "tumor__Macrophages.M0")
  violinPlotDataAll = c()
  for (i in 4:25) {
    violinPlotData = data.frame(rowname = immune_Data$rowname,value = immune_Data[i])
    violinPlotData = merge(violinPlotData , groupLabelData, by="rowname")
    violinPlotData$rowname = colnames(immune_Data)[i]
    colnames(violinPlotData) <- c("types","value","group")
    violinPlotDataAll = rbind(violinPlotDataAll ,  violinPlotData)
  }
  colnames(violinPlotDataAll) <- c("CellType","Composition","group")
  CellType = strsplit(violinPlotDataAll$CellType,"__")
  CellTypeShort = c()
  for (onetype in CellType) {
    CellTypeShort = c(CellTypeShort,onetype[2])
  }
  violinPlotDataAll$CellType <- CellTypeShort
  
  p<- ggboxplot(
    violinPlotDataAll,
    x = "CellType",
    y = "Composition",
    color = "black",
    fill = "group",
    xlab = "",
    palette = c("#015699", "#FAC00F"),
    ylab = "Cell composition",
    main = "TME Cell composition group by risk score" ) +
    stat_compare_means( aes(group = group),label = "p.signif",hide.ns =T) +
    theme(axis.text.x = element_text( angle = 90, hjust = 1, vjust = 0.5 )) 
  p
  
  ggsave("/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/3.4.2.C.immune.violin.pdf", p,width=10,heigh=7)
  
  
  
  #下面做cox
  survCompDataPart = survCompData[c("rowname", "AGE", "GENDER", "TNM_STAGE_T")]    # 取出需要的列
  colnames(survCompDataPart) <- c("patient", "age", "gender", "stage")             # 改一下列名
  survCompDataPart$age = as.numeric(survCompDataPart$age)                          # 把age改为numeric
  survCompDataPart$stage = as.character(survCompDataPart$stage)                    # 去掉factor
  survCompDataPart$stage[survCompDataPart$stage == 1] <- "Stage I"
  survCompDataPart$stage[survCompDataPart$stage == 2] <- "Stage II"
  survCompDataPart$stage[survCompDataPart$stage == 3] <- "Stage III"
  survCompDataPart$stage[survCompDataPart$stage == 4] <- "Stage IV"
  survCompDataPart = merge(y_Comp, survCompDataPart, by="patient")
  
  #A)先做批量单因素COX回归，看：risk score、age、gender、stage分别与生存之间是否相关。
  res1_Comp = unit_Cox_reg_Comp(survCompDataPart, c("group","age", "gender"))
  res2_Comp = unit_Cox_reg_Comp(survCompDataPart, c("stage"))
  result = rbind(res1_Comp, res2_Comp)
  result$name = str_remove(rownames(result), "stage")
  result = result[,c(6,1:5)]
  result = rbind(c("","","","","HR (95% CI)","P value"),result)
  result$name = c("","group (low vs. high)","age","gender (male vs.female)","TNM (II vs. I)","TNM (III vs. I)","TNM (IV vs. I)")
  pdf("/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/3.4.3.A.comp.forestPlot.pdf", onefile=F, height = 3)
  forestplot(result[,c(1,5,6)], mean=result[,2], lower=result[,3], upper=result[,4], zero=1, boxsize=0.4, graph.pos= "right",
             hrzl_lines=list("1" = gpar(lty=1,lwd=2), "2" = gpar(lty=1,columns=c(2:3)), "8"= gpar(lwd=2,lty=1,columns=c(1:3)) ),
             graphwidth = unit(.3,"npc"), xticks=c(0.25,0.5,1,1.5,2.25,4), is.summary=c(T,F,F,F,F,F,F,F),
             txt_gp=fpTxtGp( label=gpar(cex=1), ticks=gpar(cex=1), xlab=gpar(cex=1.5), title=gpar(cex=1)),
             lwd.zero=1, lwd.ci=1.5, lwd.xaxis=1.5, lty.ci=1.5, ci.vertices =T, ci.vertices.height=0.2,
             col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
             clip=c(0.25,2), lineheight = unit(8, 'mm'), line.margin=unit(8, 'mm'), colgap=unit(4, 'mm'),
             title=paste(cancerNameComp,"unicariate cox regression") )
  dev.off()
  #B)再做多变量Cox比例风险回归，看校正完三个临床信息后，risk score是否是独立影响生存的因素。主要是HR、Pvalue，统计这些数据，并做回归统计图。
  #1.前面已经载入多因素cox的函数，直接调用
  pdf("/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/3.4.3.B.comp.forestPlot.pdf", height = 3)
  multi_Cox_reg_Comp(survCompDataPart, "comp", cancerNameComp)
  dev.off()
 
  
  my_surv_df = y_Comp[c("time","status","riskScore")] %>% set_colnames(value = c("TIME_TO_EVENT","EVENT","risk"))  # 自己的模型
  rownames(my_surv_df) <- y_Comp$sample
  geneExprCompData = fread(geneExprCompFile)
  g2e=unique(toTable(org.Hs.egENSEMBL))
  colnames(g2e) <- c("Entrez_ID","ensembl_id")
  g2e$Entrez_ID = as.numeric(g2e$Entrez_ID)
  geneExprCompData = merge(geneExprCompData, g2e, by="Entrez_ID") %>% t(.)
  colnames(geneExprCompData) <- geneExprCompData[nrow(geneExprCompData),]
  geneExprCompData = geneExprCompData[-1:-2,]
  geneExprCompDataT = apply(geneExprCompData,2,as.numeric) %>% as.data.frame()
  rownames(geneExprCompDataT) = rownames(geneExprCompData)
  geneExprCompData = geneExprCompDataT[na.omit(match(y_Comp$sample,rownames(geneExprCompDataT))),]
  geneExprCompData$TIME_TO_EVENT = y_Comp$time
  geneExprCompData$EVENT = y_Comp$status
  filename = "/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/3.4.4.compareWithOthers.pdf"
  model7_list = compOtherModel(geneExprCompData, my_surv_df, cancerNameComp,filename)     # 第一个是基因表达矩阵，第二个是自己的模型
  dev.off()
  write.table(summary(model7_list,times=seq(12,36,12)),"/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/3.4.4.table.txt")
  saveRDS(geneExprCompData,"/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/geneExprCompData.rds")
  saveRDS(my_surv_df,"/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/my_surv_df.rds")
  # deciscion curves 
  aaa = survCompDataPart[c("time","status","riskScore","group","stage")]
  colnames(aaa) <- c("time","status","risk Score","ECBS.Group","TNM.Stage")
  for (i in c(1,2)) { aaa[,i] = aaa[,i] %>% as.numeric() }
  dd <- datadist(aaa)
  options(datadist="dd")
  # deciscion curves 
  nomo_score=predict(plot_nomogram(aaa)$mymodel) %>% as.data.frame() %>% set_colnames("nomo_score")
  nomodata.nomo_score=merge(aaa,nomo_score,by=0) %>% column_to_rownames(var="Row.names")
  nomodata.nomo_score = cbind(aaa,nomo_score)
  library(MASS)
  library(survival)
  nomodata.dca=nomodata.nomo_score
  nomodata.dca[nomodata.dca=="risk high"]=1
  nomodata.dca[nomodata.dca=="risk low"]=0
  nomodata.dca[nomodata.dca=="Stage I"]=1
  nomodata.dca[nomodata.dca=="Stage II"]=2
  nomodata.dca[nomodata.dca=="Stage III"]=3
  nomodata.dca[nomodata.dca=="Stage IV"]=4
  for (i in 1:NCOL(nomodata.dca)){
    nomodata.dca[,i]=as.numeric(nomodata.dca[,i])
  }
  # 画ROC曲线
  colorName = c("#4DAF4A","#E41A1C","#1E90FF")
  y_roc = nomodata.dca
  filename = "/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI"
  each_yr_roc_mul(y_roc,colorName,cancerNameComp,filename)
  
}








#################################################################################################################################################
{
  source("/home/bioinfo/TCGA/pancancer/Script/6.func.R")
  # 换个数据
  cancerNameComp = "GSE14520"                      #给个名字
  # 先用lasso得到的基因和稀疏，建立新的模型
  setwd("/home/bioinfo/TCGA/pancancer/data")   #把目录设置到存放数据的位置
  geneExprCompFile = "GSE14520-GPL3921.gene.txt" #其他类型癌症的基因表达量
  immuneCompFile = "GSE14520_immune.rds"           #其他类型癌症的免疫文件
  sampleCompFile = "HCCDB6.sample.txt"        #其他癌症的病人ID
  survCompFile = "GSE14520_clinical.txt"         #其他类型癌症的生存文件
  LIHCGeneFile = "/home/bioinfo/TCGA/pancancer/4.model_data/results/LIHC/3.4.1.A.Classifier.selectGeneData.csv"
  
  geneExprCompData = fread(geneExprCompFile)
  survCompData = fread(survCompFile)  %>% as.data.frame(.)
  survCompData = survCompData[survCompData$`Tissue Type` == "Tumor",]
  sampleCompData = fread(sampleCompFile)
  sampleCompData = t(sampleCompData) %>% as.data.frame(.) %>% rownames_to_column(.) %>% .[-1,]
  colnames(sampleCompData) <- c("sample", "TYPE", "SAMPLE_NAME", "GEO_ID", "SAMPLE_NAME_1", "PATIENT1", "PATIENT", "patient")
  sampleCompData = sampleCompData[!sampleCompData$TYPE == "Adjacent",]   # 去掉癌旁
  # 载入GSE14520得到的基因和系数
  selectGeneData = read.csv(LIHCGeneFile)
  colnames(selectGeneData) = c("gene", "coef", "Symbol")
  geneExprCompData = merge(selectGeneData,geneExprCompData,by="Symbol")
  # 用survminer包的函数surv_cutpoint来找到risk score的最佳cutoff值，
  # 不同cutoff时标准化的log-rank也会做一个统计量分布图，虚线为选择的最佳分割点;
  riskScoreComp = 0
  for (Ind in 5:ncol(geneExprCompData)) {
    riskScoreComp[Ind] = sum(geneExprCompData$coef * geneExprCompData[,Ind])
  }
  riskScoreComp = riskScoreComp[-1:-4]
  y_Comp = data.frame( sample = colnames(geneExprCompData)[-1:-4], riskScore = riskScoreComp ) 
  y_Comp = merge(y_Comp, sampleCompData[c("sample","patient","PATIENT")], by="sample")
  survCompDataPart = survCompData[c("LCS ID","Survival months","Survival status",
                                    "Age","Gender","TNM staging","HBV viral status",
                                    "ALT(>/<=50U/L)","Multinodular","AFP (>/<=300ng/ml)","Cirrhosis")]
  colnames(survCompDataPart) <- c("PATIENT", "time", "status","age","gender","stage","HBV","ALT","Multinodular","AFP","Cirrhosis")
  survCompDataPart$PATIENT = substr(survCompDataPart$PATIENT,1,7) %>% gsub("_","-",.)
  y_Comp = merge(y_Comp, survCompDataPart, by="PATIENT")
  y_Comp = y_Comp[-which(is.na(y_Comp$age)),]
  # 寻找cutpoint
  res.cut_Comp=survminer::surv_cutpoint(y_Comp, time = "time", event = "status", minprop = 0.1, variables=c("riskScore"), progressbar = F)
  cut_Comp= as.data.frame(summary(res.cut_Comp))["riskScore","cutpoint"]
  pdf("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/3.4.1.d.comp.riskScore_cutPoint.pdf",height=6.5) 
  cut1_Comp = as.data.frame(summary(res.cut_Comp))["riskScore","cutpoint"]
  plot(res.cut_Comp$riskScore, pch=20, cex=1.4, xlab = "", ylab = "Standardized log-rank statistic") 
  text(x=cut1_Comp, y=res.cut_Comp$riskScore$stats %>% min(.) + 0.1, pos=4, labels = round(cut1_Comp,3))
  dev.off()
  y_Comp$group=ifelse(y_Comp$riskScore < cut1_Comp,"risk low","risk high")
  # e)预后模型评分在病人间的分布（每个病人的risk score是多少，其中高组和低组用不同的颜色表示）也要做一个图；
  library(ggsci); library(patchwork)
  plotRiskScore(y_Comp,"comp",cancerNameComp)
  #############################
  #2）两个预后风险评估模型之间的对比分析
  # 画ROC曲线
  library(survivalROC); library(timeROC); library(timereg)
  colorName = c("#1E90FF","#FF4500")
  pdf("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/3.4.2.A.comp.ROC.pdf") 
  each_yr_roc(y_Comp[c("time","status","riskScore","group")],colorName)
  dev.off()
  # KM生存分析
  survivalResultComp = survfit(Surv(time, status) ~ group, y_Comp)
  surv_By_RiskScore_Comp(survivalResultComp, "comp", cancerNameComp)
  cosResult = coxph(formula = Surv(time, status) ~ group, data = y_Comp)
  print(summary(cosResult))
  data.survdiff <- survdiff(Surv(time, status) ~ group,y_Comp)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
  cat(paste("HR:",HR,"up95:",up95,"low95",low95))
  KMdata = paste("HR:",HR,"up95:",up95,"low95",low95)
  write.table(KMdata, "/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/KMdata.txt",row.names = F)
  
  
  #热图不好画小提琴
  
  groupLabelData = y_Comp[,c("group","PATIENT")]
  colnames(groupLabelData) <- c("group","rowname")
  
  immune_Data = readRDS(immuneCompFile)
  immune_Data = rownames_to_column(immune_Data)
  sigpara = c("tumor__T.cells.CD4.naive","tumor__T.cells.gamma.delta",
              "tumor__Monocytes","tumor__Macrophages.M0",
              "tumor__Dendritic.cells.activated","tumor__Mast.cells.resting",
              "tumor__Mast.cells.activated")
  violinPlotDataAll = c()
  for (i in 4:25) {
    violinPlotData = data.frame(rowname = immune_Data$rowname,value = immune_Data[i])
    violinPlotData = merge(violinPlotData , groupLabelData, by="rowname")
    violinPlotData$rowname = colnames(immune_Data)[i]
    colnames(violinPlotData) <- c("types","value","group")
    violinPlotDataAll = rbind(violinPlotDataAll ,  violinPlotData)
  }
  colnames(violinPlotDataAll) <- c("CellType","Composition","group")
  CellType = strsplit(violinPlotDataAll$CellType,"__")
  CellTypeShort = c()
  for (onetype in CellType) {
    CellTypeShort = c(CellTypeShort,onetype[2])
  }
  CellTypeShort = gsub("[.]"," ",CellTypeShort)
  violinPlotDataAll$CellType <- CellTypeShort
  
  
  p<- ggboxplot(
    violinPlotDataAll,
    x = "CellType",
    y = "Composition",
    color = "black",
    fill = "group",
    xlab = "",
    palette = c("#015699", "#FAC00F"),
    ylab = "Cell composition",
    main = "TME Cell composition group by risk score" ) +
    stat_compare_means( aes(group = group),label = "p.signif",hide.ns =T) +
    theme(axis.text.x = element_text( angle = 90, vjust = 0.5,hjust = 1 )) 
  p
  
  ggsave("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/3.4.2.C.immune.violin.pdf", p,width = 10,height = 7)
  
  
  # cox
  survCompDataPart = y_Comp
  survCompDataPart$stage[survCompDataPart$stage == "I"   ] <- "Stage I"
  survCompDataPart$stage[survCompDataPart$stage == "II"  ] <- "Stage II"
  survCompDataPart$stage[survCompDataPart$stage == "III" ] <- "Stage III"
  survCompDataPart$stage[survCompDataPart$stage == "IIIA"] <- "Stage III"
  survCompDataPart$stage[survCompDataPart$stage == "IIIB"] <- "Stage III"
  survCompDataPart$stage[survCompDataPart$stage == "IIIC"] <- "Stage III"
  survCompDataPartStage = survCompDataPart[!survCompDataPart$stage == ".",]

  survCompDataPart$HBV[survCompDataPart$HBV == "N"] = "."
  survCompDataPartHBV =  survCompDataPart[!survCompDataPart$HBV == ".",]
  survCompDataPartAFP =  survCompDataPart[!survCompDataPart$AFP == ".",]
  #A)先做批量单因素COX回归，看：risk score、age、gender、stage分别与生存之间是否相关。
  res1_Comp = unit_Cox_reg_Comp(survCompDataPart, c("group","age", "gender","ALT","Multinodular","Cirrhosis"))
  resStage_Comp = unit_Cox_reg_Comp(survCompDataPartStage, c("stage"))
  resHBV_Comp = unit_Cox_reg_Comp(survCompDataPartHBV, c("HBV"))
  resAPF_Comp = unit_Cox_reg_Comp(survCompDataPartAFP, c("AFP"))
  result = rbind(res1_Comp, resHBV_Comp, resAPF_Comp, resStage_Comp)
  result$name = str_remove(rownames(result), "stage")
  result = result[,c(6,1:5)]
  result = rbind(c("","","","","HR (95% CI)","P value"),result)
  result$name = c("","group (low vs. high)","age","gender (male vs.female)",
                  "ALT (low vs. high)","Multinodular (Y vs. N)",
                  "Cirrhosis (Y vs. N)","HBV (CC vs. AVR-CC)",
                  "AFP (low vs. high)","TNM (II vs. I)","TNM (III vs. I)")
  pdf("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/3.4.3.A.comp.forestPlot.pdf", onefile=F, height = 5)
  forestplot(result[,c(1,5,6)], mean=result[,2], lower=result[,3], upper=result[,4],
             zero=1, boxsize=0.4, graph.pos= "right",
             hrzl_lines=list("1" = gpar(lty=1,lwd=2), 
                             "2" = gpar(lty=1,columns=c(2:3)), 
                             "12"= gpar(lwd=2,lty=1,columns=c(1:3)) ),
             graphwidth = unit(.3,"npc"), xticks=c(0.25,0.5,1,1.5,2.25,4), 
             is.summary=c(T,F,F,F,F,F,F,F,F,F,F),
             txt_gp=fpTxtGp( label=gpar(cex=1), ticks=gpar(cex=1), 
                             xlab=gpar(cex=1.5), title=gpar(cex=1)),
             lwd.zero=1, lwd.ci=1.5, lwd.xaxis=1.5, lty.ci=1.5, 
             ci.vertices =T, ci.vertices.height=0.2,
             col=fpColors(box="royalblue",line="darkblue", summary="royalblue", 
                          hrz_lines = "#444444"),
             clip=c(0.25,2), lineheight = unit(8, 'mm'), line.margin=unit(8, 'mm'), 
             colgap=unit(4, 'mm'),
             title=paste(cancerNameComp,"unicariate cox regression") )
  dev.off()



  
  #B)再做多变量Cox比例风险回归，看校正完三个临床信息后，risk score是否是独立影响生存的因素。主要是HR、Pvalue，统计这些数据，并做回归统计图。
  #1.前面已经载入多因素cox的函数，直接调用
  # 先把无用数据剔除
  survCompDataPart = survCompDataPart[!survCompDataPart$stage == ".",]
  survCompDataPart =  survCompDataPart[!survCompDataPart$HBV == ".",]
  survCompDataPart =  survCompDataPart[!survCompDataPart$AFP == ".",]
  pdf("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/3.4.3.B.comp.forestPlot.pdf", height = 5)
  multi_Cox_reg_CompM(survCompDataPart, "comp", cancerNameComp)
  dev.off()

  
  my_surv_df = y_Comp[c("time","status","riskScore")] %>% set_colnames(value = c("TIME_TO_EVENT","EVENT","risk"))  # 自己的模型
  rownames(my_surv_df) <- y_Comp$sample
  geneExprCompData = fread(geneExprCompFile)
  g2e=unique(toTable(org.Hs.egENSEMBL))
  colnames(g2e) <- c("Entrez_ID","ensembl_id")
  g2e$Entrez_ID = as.numeric(g2e$Entrez_ID)
  geneExprCompData = merge(geneExprCompData, g2e, by="Entrez_ID") %>% t(.)
  colnames(geneExprCompData) <- geneExprCompData[nrow(geneExprCompData),]
  geneExprCompData = geneExprCompData[-1:-2,]
  geneExprCompDataT = apply(geneExprCompData,2,as.numeric) %>% as.data.frame()
  rownames(geneExprCompDataT) = rownames(geneExprCompData)
  geneExprCompData = geneExprCompDataT[na.omit(match(y_Comp$sample,rownames(geneExprCompDataT))),]
  geneExprCompData$TIME_TO_EVENT = y_Comp$time
  geneExprCompData$EVENT = y_Comp$status
  filename = filename = paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/",cancerNameComp,"/3.4.4.compareWithOthers.pdf",sep = "")
  model7_list = compOtherModel(geneExprCompData, my_surv_df, cancerNameComp, filename)     # 第一个是基因表达矩阵，第二个是自己的模型
  dev.off()
  write.table(summary(model7_list,times=seq(12,36,12)),"/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/3.4.4.table.txt")
  saveRDS(geneExprCompData,"/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/geneExprCompData.rds")
  saveRDS(my_surv_df,"/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/my_surv_df.rds")
  # deciscion curves 
  aaa = survCompDataPartStage[c("time","status","riskScore","group","stage")]
  colnames(aaa) <- c("time","status","risk Score","ECBS.Group","TNM.Stage")
  for (i in c(1,2)) { aaa[,i] = aaa[,i] %>% as.numeric() }
  dd <- datadist(aaa)
  options(datadist="dd")
  # deciscion curves 
  nomo_score=predict(plot_nomogram(aaa)$mymodel) %>% as.data.frame() %>% set_colnames("nomo_score")
  nomodata.nomo_score=merge(aaa,nomo_score,by=0) %>% column_to_rownames(var="Row.names")
  nomodata.nomo_score = cbind(aaa,nomo_score)
  library(MASS)
  library(survival)
  nomodata.dca=nomodata.nomo_score
  nomodata.dca[nomodata.dca=="risk high"]=1
  nomodata.dca[nomodata.dca=="risk low"]=0
  nomodata.dca[nomodata.dca=="Stage I"]=1
  nomodata.dca[nomodata.dca=="Stage II"]=2
  nomodata.dca[nomodata.dca=="Stage III"]=3
  for (i in 1:NCOL(nomodata.dca)){
    nomodata.dca[,i]=as.numeric(nomodata.dca[,i])
  }
  # 画ROC曲线
  colorName = c("#4DAF4A","#E41A1C","#1E90FF")
  y_roc = nomodata.dca
  filename = "/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520"
  each_yr_roc_mul(y_roc,colorName,cancerNameComp,filename)
  
  
}



#############################
# 两个模型的验证
geneExprGSE = readRDS("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/geneExprCompData.rds")
my_surf_GES = readRDS("/home/bioinfo/TCGA/pancancer/4.model_data/results/GSE14520/my_surv_df.rds")
geneExprLIRI = readRDS("/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/geneExprCompData.rds")
my_surf_LIRI = readRDS("/home/bioinfo/TCGA/pancancer/4.model_data/results/LIRI/my_surv_df.rds")
geneExprLIRI = geneExprLIRI[,unique(na.omit(match(colnames(geneExprGSE),colnames(geneExprLIRI))))]
geneExprGSE = geneExprGSE[,na.omit(match(colnames(geneExprLIRI),colnames(geneExprGSE)))]
geneExpr = rbind(geneExprLIRI, geneExprGSE)
my_surv_df = rbind(my_surf_LIRI, my_surf_GES)
filename = paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/3.4.4.compareWithOthers.pdf",sep = "")
model7_list = compOtherModel(geneExpr, my_surv_df, cancerNameComp,filename)
dev.off()


##########################################################################################################################
{
  source("/home/bioinfo/TCGA/pancancer/Script/6.func.R")
  KMdata = ""
  splots = list()
  TCGACancers = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
                  "KICH","KIRC","KIRP","LAML","LGG","LUAD","LUSC","MESO","OV","PAAD",
                  "PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
  # # 选出显著性的数据
  TCGACancers = c("ACC","BRCA","COAD","ESCA","HNSC","KICH","KIRC","KIRP","LAML","LGG",
                  "LUAD","MESO","OV","PAAD","PRAD","SARC","SKCM","UCEC","UCS","UVM")
  TCGACancers = c("ACC","BRCA","HNSC","KICH","KIRC","KIRP","LAML","LGG","LUAD","LUSC",
                  "MESO","OV","PAAD","PCPG","SARC","SKCM","UCEC","UVM")
  for (cancerName in TCGACancers) {
    
    cat(paste0(cancerName,"\n"))
    # 换TCGA的数据
    cancerNameComp = cancerName                      #给个名字
    # 先用lasso得到的基因和稀疏，建立新的模型
    setwd("/home/bioinfo/TCGA/pancancer/data")   #把目录设置到存放数据的位置
    geneExprCompFile = paste("/home/bioinfo/FPKM_cli/TCGA_express_matrix/TCGA-",cancerNameComp,"_raw_exp.FPKM.tsv", sep = "") #其他类型癌症的基因表达量
    survCompFile = '/home/bioinfo/Music/REIA/multi_omics/tcga_patient_clinical_info.csv'         #其他类型癌症的生存文件
    
    geneExprCompData = fread(geneExprCompFile) %>% as.data.frame(.)
    colnames(geneExprCompData) <- gsub("[.]","-",colnames(geneExprCompData))
    allSamples = colnames(geneExprCompData)[-1:-26]
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
    geneExprCompData = geneExprCompData[c("gene_id",tumorSample)]          # 去掉癌旁
    survCompData = fread(survCompFile) 
    survCompData <- survCompData[survCompData$cancer == cancerNameComp,]
    survCompData = survCompData %>% as.data.frame(.) 
    # 载入LIHC得到的基因和系数
    selectGeneData = read.csv(LIHCGeneFile)
    colnames(selectGeneData) = c("gene_id", "coef", "Symbol")
    geneExprCompData = merge(selectGeneData,geneExprCompData,by="gene_id")
    # 用survminer包的函数surv_cutpoint来找到risk score的最佳cutoff值，
    # 不同cutoff时标准化的log-rank也会做一个统计量分布图，虚线为选择的最佳分割点;
    riskScoreComp = 0
    for (Ind in 4:ncol(geneExprCompData)) {
      riskScoreComp[Ind] = sum(geneExprCompData$coef * geneExprCompData[,Ind])
    }
    riskScoreComp = riskScoreComp[-1:-3]
    y_Comp = data.frame( sample = colnames(geneExprCompData)[-1:-3], riskScore = riskScoreComp ) 
    survCompDataPart = survCompData[c("patient_id","os","os_time")]
    colnames(survCompDataPart) <- c("patient", "status", "time")
    y_Comp$sample = substr(y_Comp$sample,1,12)
    colnames(y_Comp) <- c("patient","riskScore")
    y_Comp = merge(y_Comp, survCompDataPart, by="patient")
    y_Comp$time = y_Comp$time/30
    y_Comp = y_Comp[!is.na(y_Comp$time),]
    # 寻找cutpoint
    res.cut_Comp=survminer::surv_cutpoint(y_Comp, time = "time", event = "status", minprop = 0.1, variables=c("riskScore"), progressbar = F)
    cut_Comp= as.data.frame(summary(res.cut_Comp))["riskScore","cutpoint"]
    y_Comp$group=ifelse(y_Comp$riskScore < cut_Comp,"risk low","risk high")
    #############################
    # KM生存分析
    survivalResultComp = survfit(Surv(time, status) ~ group, y_Comp)
    pp = surv_By_RiskScore_TCGA(survivalResultComp, cancerNameComp)
    
    splots[[which(TCGACancers == cancerName)]] <- pp
    
    data.survdiff <- survdiff(Surv(time, status) ~ group,y_Comp)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
    cat(paste("HR:",HR,"up95:",up95,"low95",low95))
    cat("\n")
    KMdata = paste(KMdata,"cancer:",cancerName,"HR:",HR,"up95:",up95,"low95",low95,"\n")
    
  }
  write.table(KMdata, "/home/bioinfo/TCGA/pancancer/4.model_data/results/TCGA/KMdata.txt",row.names = F)
  
  # 上面的尺寸适合放论文
  p <- arrange_ggsurvplots(splots,ncol = 4, nrow = 5) #定义行数和列数
  fileName = paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/TCGA/TCGA_ALL_KM_survival.pdf",sep = '')
  ggsave(fileName, p, width=21, height = 32)
  # 下面的尺寸适合放PPT
  p1 <- arrange_ggsurvplots(splots,ncol = 10, nrow = 2) #定义行数和列数
  fileName = paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/TCGA/TCGA_ALL_KM_survival_wide.pdf",sep = '')
  ggsave(fileName, p1, width=30, height = 14)
  
}