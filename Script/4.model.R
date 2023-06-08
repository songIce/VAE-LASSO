#1. VAE压缩+无监督聚类
#1)这里改成用VAE压缩后的矩阵再聚类
#设置一些公共参数
{
  kCluster = 3
  cancerName = "LIHC"
  #先处理一下得到原始矩阵
  aData = fread(paste("/home/bioinfo/TCGA/pancancer/02_difference_of_tumor_normal/04_/04_",cancerName,"_tumor_RESs.csv",sep = ''))
  #先导入stage和survival和peptide的数据
  stageData = fread(paste('/home/bioinfo/TCGA/pancancer/CancerWithpValue/','allStagepValue.csv',sep = '')) 
  ostimeData = fread(paste('/home/bioinfo/TCGA/pancancer/CancerOstimepValue/','all_ostime_pvalue.csv',sep = ''))
  peptideData = fread("/home/bioinfo/TCGA/pancancer/CancerPeptide/positionAndCancer.csv")
  stageData = stageData[stageData$cancer==cancerName] #太少了没必要去掉了吧
  ostimeData = ostimeData[ostimeData$cancer==cancerName] 
  peptideData = peptideData[peptideData$cancer==cancerName] #这个也是太少了
  #去重(要对两个数据简单处理一下，结构不同)
  paste(stageData$Chr,stageData$`Position(1base)`,sep = '') -> stagePosition
  paste(ostimeData$Chr,ostimeData$`Position(1base)`,sep = '') -> ostimePosition
  paste(aData$Chr,aData$`Position(1base)`,sep = '') -> cancerPosition
  peptideData$position -> peptidePosition
  c(stagePosition,ostimePosition,cancerPosition,peptidePosition) %>% unique() -> allPositUniq
  #读入并处理cancer数据
  #1. 根据EF——info数据中的列名，区分出tumor和normal，两者的区别在sampleID的第14、15位，小于等于9 的是tumor，大于9的则是normal。
  fread(paste("/home/bioinfo/TCGA/pancancer/informative_RES/",cancerName,"_EF_informative_no_SNP_2.csv",sep = '')) -> cancerData 
  if (cancerData$Chr[nrow(cancerData)] %>% is.na()) {
    cancerData = cancerData[-nrow(cancerData),]   # 删掉最后一列统计的
  }
  aDataPosition = paste(cancerData$Chr,cancerData$`Position(1base)`,sep = '') 
  match(allPositUniq,aDataPosition) %>% na.omit(.) %>% unique(.)  %>% cancerData[.,] -> matchData
  tumorSample = colnames(matchData)[1:26]
  for(oneColname in colnames(matchData)[-1:-26]) {
    name_list <- oneColname %>% substring(14,15) %>% as.numeric()
    if(name_list>=1&name_list<=9) {
      if (substr(oneColname,16,16)=="A") {
        tumorSample = c(tumorSample,oneColname)
      }
    }
  }
  substr(tumorSample,1,12)[which(duplicated(substr(tumorSample,1,12)))] -> dupSample
  for (oneSample in dupSample) {
    which(substr(tumorSample,1,12)==oneSample) %>%  tumorSample[.] -> dupS
    substr(dupS,15,15) -> sa
    if (sa[1]!=sa[2]) {
      tumorSample = tumorSample[-match(dupS[1],tumorSample)]
      next
    }
    substr(dupS,22,25) %>% max(.) ->maxV
    if (!is.na(maxV)) {
      dupS[which(substr(dupS,22,25)!=maxV)] -> delS
      tumorSample = tumorSample[-match(delS,tumorSample)]
    }
  }
  
  
  #tumorData
  matchData = matchData[,..tumorSample]
  #好多值大于1？最大的还超过一万，必须置1，好像只能一列一列做？
  for (i in 27:ncol(matchData)) {
    matchData[which(matchData[,..i]>1),i] = 1
  }
  saveRowName = paste(matchData$Chr,matchData$`Position(1base)`,sep = '')
  write.csv(matchData,file = paste("/home/bioinfo/TCGA/pancancer/VAE_data/",cancerName,"/",cancerName,"all_matchData.csv",sep = ''),row.names = F)
  matchData = matchData[,-1:-26] %>% as.matrix(.) %>% t(.)
  colnames(matchData) <- saveRowName
  #保存一下处理完的LIHC数据
  # matchData = scale(matchData)
  write.csv(matchData,file = paste("/home/bioinfo/TCGA/pancancer/VAE_data/",cancerName,"/",cancerName,"_matchData.csv",sep = ''))
  
}
####################################
#去做VAE









#设置一些公共参数
{
  kCluster = 3           # 聚类数量
  cancerName = "LIHC"
  set.seed(1)
  
  # 载入需要使用的文件路径
  setwd(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/",cancerName,sep = ''))
  vae_data_file_name = paste("/home/bioinfo/TCGA/pancancer/4.model_data/VAE_data/",cancerName,"/1.",cancerName,".latent.tsv",sep = '')
  clinical_data_file_name = '/home/bioinfo/TCGA/pancancer/4.model_data/tcga_patient_clinical_info.csv'
  AEI_data_file_name = paste("/home/bioinfo/TCGA/pancancer/4.model_data/AEI/",cancerName,"/",cancerName,"_EditingIndex.csv",sep = '')
  # Immune_data_file_name = "/home/bioinfo/TCGA/pancancer/4.model_data/CIBERSORTx_Job2_Results.csv"
  Immune_data_file_name = "/home/bioinfo/TCGA/pancancer/4.model_data/immune_data.csv"
  # 免疫的源文件做了预处理，主要是修改列名和patiend id
  GeneExpr_data_filename = "/home/bioinfo/TCGA/pancancer/4.model_data/LIHC_geneExpr_count_New.csv"
  GeneExpr_fpkm_filename = "/home/bioinfo/TCGA/pancancer/4.model_data/fpkm/TCGA-LIHC_htseq_fpkm.tsv"
  
  # 载入公共函数
  source("/home/bioinfo/TCGA/pancancer/Script/6.func.R")
  source('/home/bioinfo/TCGA/pancancer/Script/dca.R')
  
  # 载入需要的包
  library(plyr); library(dplyr); library(data.table); 
  library(tidyr); library(marray); library(magrittr)
  library(tibble); library(ConsensusClusterPlus); 
  library(survival); library(survminer); library(glmnet)
  #主成分分析！
  library(psych); library(reshape2); library(ggplot2); library(factoextra)
  #免疫热图
  library(pheatmap)
  #构建预后模型
  library(DESeq2); library(org.Hs.eg.db); library(clusterProfiler)
  library(ggsci); library(patchwork); library(survivalROC); 
  library(timeROC); library(timereg); library(tableone)  
  ##画森林图的包
  library(forestplot); library(stringr); library(Hmisc)
  library(grid); library(lattice); library(Formula);
  library(rms); library(MASS);
  library(TCGAutils)
  
}




#下面使用VAE压缩后的数据，需要把第一列作为列名
fread( vae_data_file_name, header = T ) -> vae_data
vae_data_Names = vae_data[,1]  %>% as.data.frame(.)
vae_data = vae_data[,-1]  %>% t(.)
colnames(vae_data) <- vae_data_Names[,1]
rm(vae_data_Names)

# 聚类
{
  # browseVignettes("ConsensusClusterPlus")#查看帮助文档
  cluster_LIHC = ConsensusClusterPlus(vae_data,maxK=3,reps=1000,pItem=0.8,pFeature=1,title="1.1.clusterResults",
                                      clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="pdf",writeTable=TRUE)
  cluster_LIHC[[kCluster]]["consensusClass"] -> cluster_Results
  rm(cluster_LIHC)
}
# 聚类结束

# 开始画PCA
{
  vae_data_t = t(vae_data)
  vae_data_t = vae_data_t[ , which(apply(vae_data_t, 2, var) != 0)]
  pca_result <- prcomp(vae_data_t, scale=T)
  Letter = c("EC1","EC2","EC3","EC4","EC5")
  groupPCA = rep(1,nrow(vae_data_t))
  for (i in 1:kCluster) {
    groupPCA[match(names(cluster_Results[[1]])[cluster_Results[[1]] == i],rownames(vae_data_t))] = Letter[i]
  }
  
  
  summ<-summary(pca_result)
  xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
  ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
  
  
  pca_result$x %>% .[,1:2] %>% as.data.frame() %>% 
    cbind(Groups = groupPCA) %>%
    ggplot(aes(x = PC1, y = PC2),
           ellipse = TRUE, ellipse.prob=0.95, circle = T,var.axes=T) + 
    geom_point(aes(color = Groups),size = 2.5,alpha = 10) + 
    labs(x=xlab,y=ylab,color="") +
    stat_ellipse(aes(fill=Groups),
                 type = "norm", geom ="polygon",alpha=0.2,color=NA)+
    theme_bw(base_size = 30) -> p
  ggsave("3.3.1.A.PCA.pdf", p, width=9, height=6.5);
}
# PCA结束



#2)构建分类器：用聚出来的亚型，以及原本未压缩的editing level矩阵来做分类器。
# 做好后，看亚型之间在生存、免疫评分、AEI上的区别。
#A.生存用KM曲线分析（3组的生存曲线，如果生存时间太长，可以把生存时间超过1200天的数据过滤掉。

{
  #匹配一下临床数据，取出相同的sample
  clinical_data_full <- fread(clinical_data_file_name) %>% unique(.)
  clinical_selected <- clinical_data_full[clinical_data_full$cancer == cancerName,]
  clinical_selected <- clinical_selected[match(unique(substr(colnames(vae_data),1,12)),clinical_selected$patient_id)]
  
  #分组
  Group = rep(1,nrow(clinical_selected))
  for (i in 1:kCluster) {
    Group[match(substr(names(cluster_Results[[1]])[cluster_Results[[1]] == i] , 1 , 12),clinical_selected$patient_id)] = Letter[i]
  }
  clinical_selected = cbind(clinical_selected,Group) %>% na.omit(.)
  
  # 生存分析，画KM
  # 先画EC1,EC2,EC3三组的
  {
    survival_Result = survfit(Surv(os_time, os) ~ Group, data = clinical_selected)
    res <- pairwise_survdiff(Surv(os_time, os) ~ Group, data = clinical_selected, p.adjust.method = "BH"); res
    write.list(res, "1.2.LogRankTest.txt")
    p <- ggsurvplot(survival_Result, pval = T, conf.int = FALSE, risk.table = TRUE, 
                    risk.table.col = "strata", linetype = "strata", surv.median.line = "none", 
                    xlim = c(0,2000), ggtheme = theme_bw(base_size = 20), legend.labs = c("EC1", "EC2", "EC3"),) 
    pdf("3.3.1.B.LIHC_KM_survival.pdf", onefile=F, width=7, height = 9)
    print(p)
    dev.off()
    
    data.survdiff <- survdiff(Surv(os_time, os) ~ Group, clinical_selected)
    HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
    cat(paste("EC2/EC1","HR:",HR,"up95:",up95,"low95",low95))
    cat("\n")
    HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[3]/data.survdiff$exp[3])
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[3]))
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[3])) 
    cat(paste("EC2/EC3","HR:",HR,"up95:",up95,"low95",low95))
    cat("\n")
    
  }
  # 再画EC-high和EC-low两组的
  {
    cli_EC = clinical_selected
    cli_EC$Group[cli_EC$Group=="EC1"] = "Non-EC2"
    cli_EC$Group[cli_EC$Group=="EC3"] = "Non-EC2"
    cli_EC$Group[cli_EC$Group=="EC2"] = "EC2"
    survival_Result = survfit(Surv(os_time, os) ~ Group, data = cli_EC)
    res <- pairwise_survdiff(Surv(os_time, os) ~ Group, data = clinical_selected, p.adjust.method = "BH"); res
    write.list(res, "1.2.LogRankTest.EC.txt")
    p <- ggsurvplot(survival_Result, pval = T, conf.int = FALSE, risk.table = TRUE, 
                    palette = c("#015699", "#FAC00F"),
                    risk.table.col = "strata", linetype = "strata", surv.median.line = "none", 
                    xlim = c(0,2000), ggtheme = theme_bw(base_size = 20), 
                    legend.labs = c("EC2", "Non-EC2"),) 
    pdf("3.3.1.B.LIHC_KM_survival.EC.pdf", onefile=F, width=7.5, height = 9)
    print(p)
    dev.off()
  }
  
  # 下面是输出一些具体的数据
  data.survdiff <- survdiff(Surv(os_time, os) ~ Group,clinical_selected)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
  cat(paste("HR:",HR,"up95:",up95,"low95",low95))
  
  #B.AEI上的区别，柱形图，看整体编辑率的差异性。
  AEI_data = fread(AEI_data_file_name)
  AEI_data[AEI_data$StrandDecidingMethod=="RefSeqThenMMSites"] -> AEI_data
  AEI_data$SamplePath  %>% strsplit(.,"[/]") -> sampleName
  allSamples = data.frame()
  for (oS in sampleName) {
    UUIDs = UUIDhistory(oS[7]) 
    newUUID = UUIDs$uuid[UUIDs$file_change=="released"]
    oneSample = data.frame( Sample = oS[7], UUID = newUUID )
    allSamples = rbind( allSamples, oneSample )
  }
  
  UUID_set = UUIDtoBarcode(allSamples$UUID, from_type = "file_id")
  colnames(UUID_set) <- c("UUID","Barcode")
  UUID_set$UUID <- as.character(UUID_set$UUID)
  
  AEI_data$UUID_old = allSamples$Sample
  AEI_data$UUID_new = allSamples$UUID

  write.csv(merge(UUID_set,allSamples,by="UUID"),"UUID_Sample_Barcode.csv")

  # 画EC1 EC2 EC3三组的
  {
    boxplotAEI = data.frame()
    for (i in 1:kCluster) {
      groupUUID = UUID_set$UUID[match(names(cluster_Results[[1]])[cluster_Results[[1]] == i],UUID_set$Barcode)] 
      groupData = AEI_data$A2GEditingIndex[match(groupUUID,AEI_data$UUID_new)]
      boxplotAEI = rbind(boxplotAEI,cbind(groupData,rep(Letter[i],length(groupData))))
    }
    colnames(boxplotAEI) <- c("AEI value","group")
    boxplotAEI$`AEI value` <- boxplotAEI$`AEI value` %>% as.numeric()
    boxplotAEI$group <- boxplotAEI$group %>% as.character() %>% as.factor()
    compaired <- list(c("EC1", "EC2"),c("EC2","EC3"),c("EC1","EC3"))
    p <- ggboxplot(boxplotAEI,x="group",y="AEI value",color = "group", 
                   palette=c("#f8766d","#00ba38","#619cff"), add = "jitter", xlab = "", ylab = "AEI (%)",
                   ggtheme = theme_pubr(base_size = 15))
    p <- p + stat_compare_means(comparisons = compaired)
    ggsave("3.3.1.B.AEI_BoxPlot.pdf",p,height = 6,width = 7)
  }
  
  
  # 画normal EC1 EC2 EC3四组的
  {
    boxplotAEI = data.frame()
    for (i in 1:kCluster) {
      groupUUID = UUID_set$UUID[match(names(cluster_Results[[1]])[cluster_Results[[1]] == i],UUID_set$Barcode)] 
      groupData = AEI_data$A2GEditingIndex[match(groupUUID,AEI_data$UUID_new)]
      boxplotAEI = rbind(boxplotAEI,cbind(groupData,rep(Letter[i],length(groupData))))
    }
    groupUUIDNor = UUID_set$UUID[match(names(cluster_Results[[1]]),UUID_set$Barcode)] 
    groupData = AEI_data$A2GEditingIndex[- match(groupUUIDNor,AEI_data$UUID_new)]
    boxplotAEI = rbind(boxplotAEI,cbind(groupData,rep("Normal",length(groupData))))
    colnames(boxplotAEI) <- c("AEI value","group")
    boxplotAEI$`AEI value` <- boxplotAEI$`AEI value` %>% as.numeric()
    boxplotAEI$group <- boxplotAEI$group %>% as.character() %>% as.factor()
    compaired <- list(c("Normal","EC3"),c("EC1", "EC2"),c("EC2","EC3"),
                      c("EC1","EC3"),c("Normal","EC2"),c("Normal", "EC1"))
    p <- ggboxplot(boxplotAEI,x="group",y="AEI value",color = "group", 
                   palette=c("#f8766d","#00ba38","#619cff","Gray"), add = "jitter", xlab = "", ylab = "AEI (%)",
                   ggtheme = theme_pubr(base_size = 15))
    p <- p + stat_compare_means(comparisons = compaired)
    ggsave("3.3.1.B.AEI_BoxPlot.normal.pdf",p,height = 6,width = 7)
  }
  
  # 画EC2和Non-EC2两组的
  {
    
    boxplotAEI = data.frame()
    for (i in 1:kCluster) {
      groupUUID = UUID_set$UUID[match(names(cluster_Results[[1]])[cluster_Results[[1]] == i],UUID_set$Barcode)] 
      groupData = AEI_data$A2GEditingIndex[match(groupUUID,AEI_data$UUID_new)]
      boxplotAEI = rbind(boxplotAEI,cbind(groupData,rep(Letter[i],length(groupData))))
    }
    colnames(boxplotAEI) <- c("AEI value","group")
    boxplotAEI$group[boxplotAEI$group=="EC1"] = "Non-EC2"
    boxplotAEI$group[boxplotAEI$group=="EC3"] = "Non-EC2"
    boxplotAEI$`AEI value` <- boxplotAEI$`AEI value` %>% as.numeric()
    boxplotAEI$group <- boxplotAEI$group %>% as.character() %>% as.factor()
    p <- ggboxplot(boxplotAEI,x="group",y="AEI value",color = "group", 
                   palette=c("#015699", "#FAC00F"), add = "jitter", xlab = "", ylab = "AEI (%)",
                   ggtheme = theme_pubr(base_size = 15))
    p <- p + stat_compare_means(comparisons = list(c("EC2","Non-EC2")))
    ggsave("3.3.1.B.AEI_BoxPlot.EC.pdf",p,height = 6,width = 7)
    
  }
  
  # 画normal EC2和Non-EC2三组的
  {
    
    boxplotAEI = data.frame()
    for (i in 1:kCluster) {
      groupUUID = UUID_set$UUID[match(names(cluster_Results[[1]])[cluster_Results[[1]] == i],UUID_set$Barcode)] 
      groupData = AEI_data$A2GEditingIndex[match(groupUUID,AEI_data$UUID_new)]
      boxplotAEI = rbind(boxplotAEI,cbind(groupData,rep(Letter[i],length(groupData))))
    }
    groupUUIDNor = UUID_set$UUID[match(names(cluster_Results[[1]]),UUID_set$Barcode)] 
    groupData = AEI_data$A2GEditingIndex[- match(groupUUIDNor,AEI_data$UUID_new)]
    boxplotAEI = rbind(boxplotAEI,cbind(groupData,rep("Normal",length(groupData))))
    colnames(boxplotAEI) <- c("AEI value","group")
    boxplotAEI$group[boxplotAEI$group=="EC1"] = "Non-EC2"
    boxplotAEI$group[boxplotAEI$group=="EC3"] = "Non-EC2"
    boxplotAEI$`AEI value` <- boxplotAEI$`AEI value` %>% as.numeric()
    boxplotAEI$group <- boxplotAEI$group %>% as.character() %>% as.factor()
    compaired <- list(c("EC2","Non-EC2"),c("Normal", "Non-EC2"),c("EC2","Normal"))
    p <- ggboxplot(boxplotAEI,x="group",y="AEI value",color = "group", 
                   palette=c("#015699", "#FAC00F","Gray"), add = "jitter", xlab = "", ylab = "AEI (%)",
                   ggtheme = theme_pubr(base_size = 15))
    p <- p + stat_compare_means(comparisons = compaired)
    ggsave("3.3.1.B.AEI_BoxPlot.EC.normal.pdf",p,height = 6,width = 7)
    
  }
  
  
  #C.免疫评分：数据：/home/bioinfo/Music/REIA/multi_omics/tcga_immune_cell_infiltration.csv
  # 根据sampleID找出数据中对应样本的免疫细胞组成成分值，然后再看不同亚型间的区别
  # (参考REIA上的分析，得到每组有3个柱的柱形图)
  
  # 先画EC1,EC2,EC3三组的
  {
    immune_Data = fread(Immune_data_file_name) 
    groupPosi = c()
    allSam = cluster_Results[[1]] %>% sort()
    groupLabel = data.frame(sample = substr(names(allSam),1,12), group =  as.vector(allSam) 
                            %>% gsub("1","EC1",.) %>% gsub("2","EC2",.) %>% gsub("3","EC3",.) %>% as.factor(.) ) %>% unique(.)
    remove_rownames(groupLabel) %>% column_to_rownames(.,var = "sample")  -> groupLabel
    for (i in 1:kCluster) {
      groupPosi = c(groupPosi,match(unique(substr(names(cluster_Results[[1]])[cluster_Results[[1]] == i] , 1 , 12)),immune_Data$patient_id))
    }
    immune_Data = immune_Data[groupPosi]
    immune_Data = column_to_rownames(immune_Data,"patient_id") %>% t(.)
    
    # immune_Data = scale(immune_Data)
    fileName = "3.3.1.C.immune.pdf"
    pheatmap(immune_Data,show_colnames = F, filename = fileName, annotation_col = groupLabel, 
             scale = "row",
             fontsize_row = 8, cluster_cols = F, 
             colorRampPalette(colors = c("dodgerblue4","white","red"))(50) )
  }
  # 再画EC-high和EC-low两组的
  {
    immune_Data = fread(Immune_data_file_name) 
    groupPosi = c()
    allSam = cluster_Results[[1]] %>% sort()
    groupL_EC = data.frame(sample = substr(names(allSam),1,12), group =  as.vector(allSam) 
                            %>% gsub("2","EC2",.) %>% gsub("1","Non-EC2",.) %>% gsub("3","Non-EC2",.) %>% as.factor(.) ) %>% unique(.)
    remove_rownames(groupL_EC) %>% column_to_rownames(.,var = "sample")  -> groupL_EC
    for (i in c(1,3,2)) {
      groupPosi = c(groupPosi,match(unique(substr(names(cluster_Results[[1]])[cluster_Results[[1]] == i] , 1 , 12)),immune_Data$patient_id))
    }
    immune_Data = immune_Data[groupPosi]
    immune_Data = column_to_rownames(immune_Data,"patient_id") %>% t(.)
    
    # immune_Data = scale(immune_Data)
    fileName = "3.3.1.C.immune.EC.pdf"
    pheatmap(immune_Data,show_colnames = F, filename = fileName, annotation_col = groupL_EC, 
             scale = "row",
             fontsize_row = 8, cluster_cols = F, 
             colorRampPalette(colors = c("dodgerblue4","white","red"))(50) )
  }
  
  # 3）如果上面的3组间免疫浸润热图效果不好的话，就用这部分分组数据，增加一个组间免疫浸润水平差距分析的小提琴图。
  # 之前做过两组间的热图，效果不理想，现在可以将26种免疫浸润数据一一对比，其中显著性差异的挑出来做成小提琴图，
  # 每组2个提琴代表一个免疫水平，有几个差异就做几组，放在一个大图上。
  
  #先做一个显著的图，再做一个全部的图
  immune_Data = fread(Immune_data_file_name)
  violinPlotDataAll = c()
  for (i in c(3,15)) {
    violinPlotData = data.frame(rowname = immune_Data$patient_id,value = immune_Data[,..i])
    groupLabelData = rownames_to_column(groupLabel)
    violinPlotData = merge(violinPlotData , groupLabelData, by="rowname")
    violinPlotData$rowname = colnames(immune_Data)[i]
    colnames(violinPlotData) <- c("CellType","Composition","group")
    violinPlotDataAll = rbind(violinPlotDataAll ,  violinPlotData)
  }
  compaired <- list(c("EC1", "EC2"),c("EC2","EC3"),c("EC1","EC3"))
  
  p <- ggboxplot(violinPlotDataAll, x="group", y="Composition", color = "group", facet.by = "CellType",
                 palette=c("#f8766d","#00ba38","#619cff")) +
    stat_compare_means(comparisons = compaired)
  ggsave("3.3.1.D.violin.sig.pdf", p)
  
  #做全部的图
  violinPlotDataAll = c()
  for (i in 2:23) {
    violinPlotData = data.frame(rowname = immune_Data$patient_id,value = immune_Data[,..i])
    groupLabelData = rownames_to_column(groupLabel)
    violinPlotData = merge(violinPlotData , groupLabelData, by="rowname")
    violinPlotData$rowname = colnames(immune_Data)[i]
    colnames(violinPlotData) <- c("CellType","Composition","group")
    violinPlotDataAll = rbind(violinPlotDataAll ,  violinPlotData)
  }
  compaired <- list(c("EC1", "EC2"),c("EC2","EC3"),c("EC1","EC3"))
  
  p <- ggboxplot(violinPlotDataAll, x="group", y="Composition", color = "group", facet.by = "CellType",
                 palette=c("#f8766d","#00ba38","#619cff"), ncol = 11, nrow = 2) +
    stat_compare_means(comparisons = compaired)
  ggsave("3.3.1.D.violin.all.pdf", p,height = 10,width = 22) 
  
}

{
  # 再做一个只有EC2和EC1+EC3比较的
  immune_Data = fread(Immune_data_file_name)
  groupLabelData = rownames_to_column(groupLabel)
  groupp = rep("EC2",length(groupLabelData$group))
  groupp[which(groupLabelData$group=="EC3")] = "Non-EC2"
  groupp[which(groupLabelData$group=="EC1")] = "Non-EC2"
  groupLabelData$group <- as.factor(groupp)
  violinPlotDataAll = c()
  for (i in 2:23) {
    violinPlotData = data.frame(rowname = immune_Data$patient_id,value = immune_Data[,..i])
    violinPlotData = merge(violinPlotData , groupLabelData, by="rowname")
    violinPlotData$rowname = colnames(immune_Data)[i]
    colnames(violinPlotData) <- c("CellType","Composition","group")
    violinPlotDataAll = rbind(violinPlotDataAll ,  violinPlotData)
  }
  
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
  
  ggsave("3.3.1.D.violin.EC2.NonEC2.pdf", p,width = 10)
  
}



# 画桑基图
library(plotly)
library(htmlwidgets)

dataType = fread("/home/bioinfo/TCGA/pancancer/data/tcga-subtype.txt")
dataLIHC = dataType[grep("LIHC",dataType$subtype)]
# LIHCiCluster1 = dataLIHC$barcode[dataLIHC$subtype == "LIHC.iCluster:1"]
# LIHCiCluster2 = dataLIHC$barcode[dataLIHC$subtype == "LIHC.iCluster:2"]
# LIHCiCluster3 = dataLIHC$barcode[dataLIHC$subtype == "LIHC.iCluster:3"]
# 
# LIHCEC2 = names(cluster_Results[[1]])[cluster_Results[[1]] == 2] %>% substr(.,1,16)
# LIHCNon-EC2 = names(cluster_Results[[1]])[cluster_Results[[1]] == 1] %>% substr(.,1,16)
#
## value的顺序
## EC2-iCluster:1, EC2-iCluster:2, EC2-iCluster:3,
## Non-EC2-iCluster:1, Non-EC2-iCluster:2, Non-EC2-iCluster:3, 

# sAsc$sample[sAsc$group == "risk high"] EC2
# sAsc$sample[sAsc$group == "risk high"] Non-EC2
sankeyValues = c()
for (i in c("risk high","risk low")) {
  LIHCEC = sAsc$sample[sAsc$group == i] %>% substr(.,1,16)
  for (j in 1:3) {
     LIHCiCluster = dataLIHC$barcode[dataLIHC$subtype == paste("LIHC.iCluster:",j,sep = "")]
     sankeyValue = match(LIHCEC,LIHCiCluster) %>% na.omit(.) %>% length(.)
    sankeyValues = c(sankeyValues,sankeyValue)
  }
}





dataType = fread("/home/bioinfo/TCGA/pancancer/data/tcga-subtype.txt")
 
dataLIHC = dataType[grep("LIHC",dataType$subtype)]
clinical <- fread('/home/bioinfo/Music/REIA/multi_omics/tcga_patient_clinical_info.csv') %>% unique(.)
clinicalLIHC <- clinical[clinical$cancer == cancerName,]

library("survminer")
EC_Group = c("EC1","EC2","EC3")
splots <- list()
iClusterGroup = c("LIHC.iCluster:1","LIHC.iCluster:2","LIHC.iCluster:3")
for (j in 1:3) {
  LIHCiCluster = dataLIHC$barcode[dataLIHC$subtype == iClusterGroup[j]]
  GroupData = data.frame(Sample = LIHCiCluster, Group = rep(1,length(LIHCiCluster)))
  for (i in EC_Group) {
      GroupData$Group[match(substr(names(cluster_Results[[1]])[cluster_Results[[1]] == which(EC_Group==i)], 1, 16), GroupData$Sample)] = i
  }
  if (is.element(1,GroupData$Group)) { GroupData = GroupData[!GroupData$Group==1,] }
  clinicalLIHCData = cbind(clinicalLIHC[match(substr(GroupData$Sample,1,12),clinicalLIHC$patient_id),],Group = GroupData$Group) %>% na.omit(.)
  survivalResult = survfit(Surv(os_time, os) ~ Group, data = clinicalLIHCData)
  res <- pairwise_survdiff(Surv(os_time, os) ~ Group, data = clinicalLIHCData, p.adjust.method = "BH"); res
  pp <- ggsurvplot(survivalResult, pval = T, conf.int = TRUE, risk.table = TRUE,
                                        risk.table.col = "strata", legend=c(0.80,0.85),
                                        linetype = "strata", surv.median.line = "none",
                                        legend.labs = c("EC1", "EC2", "EC3"), title = paste("iCluster:",j),
                                        xlim = c(0,1800), ggtheme = theme_survminer(base_size = 16))
  splots[[j]] <- pp
  print(paste("LIHC.iCluster:",j))
  print(res)
  data.survdiff <- survdiff(Surv(os_time, os) ~ Group,clinicalLIHCData)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
  cat(paste("EC2/EC1","HR:",HR,"up95:",up95,"low95",low95))
  cat("\n")
  HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2])) 
  cat(paste("EC1/EC2","HR:",HR,"up95:",up95,"low95",low95))
  cat("\n")
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[3]/data.survdiff$exp[3])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[3]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[3])) 
  cat(paste("EC2/EC3","HR:",HR,"up95:",up95,"low95",low95))
  cat("\n")
  HR = (data.survdiff$obs[3]/data.survdiff$exp[3])/(data.survdiff$obs[2]/data.survdiff$exp[2])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[3]+1/data.survdiff$exp[2]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[3]+1/data.survdiff$exp[2])) 
  cat(paste("EC3/EC2","HR:",HR,"up95:",up95,"low95",low95))
  cat("\n")
}
p <- arrange_ggsurvplots(splots,ncol = 3, nrow = 1) #定义行数和列数
ggsave("3.3.4.B.EC_in_iCluster_LIHC_KM_survival.pdf", p, width=12, height = 6.5)

ggsave("3.3.4.B.EC_in_iCluster_LIHC_KM_survival.pdf", arrange_ggsurvplots(splots,ncol = 3, nrow = 1), width=12, height = 6.5)


#保存一下所有数据，不用重新跑
save.image(file = paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,".Step1.RData",sep = ''))


cancerName = "LIHC"
load(file = paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,".Step1.RData",sep = ''))
source("/home/bioinfo/TCGA/pancancer/Script/6.func.R")
source('/home/bioinfo/TCGA/pancancer/Script/dca.R')
setwd(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/",cancerName,sep = ''))




##########################################################################################################################################
#3.构建预后模型
# 全部数据
# 取了交集的数据
geneExprNor = fread("/home/bioinfo/FPKM_cli/fpkm/TCGA-LIHC_htseq_fpkm.tsv")
geneExprNor = column_to_rownames(geneExprNor,var = "gene_id")


g2e=unique(toTable(org.Hs.egENSEMBL))
g2s=unique(toTable(org.Hs.egSYMBOL))
g2all = merge(g2e,g2s,by="gene_id")
RNA <- intersect(rownames(geneExprNor),g2all$ensembl_id)
geneExprNor = geneExprNor[RNA,]


colnames(geneExprNor) <- colnames(geneExprNor) %>% gsub("[.]","-",.)
tumorSample = c()
for(oneColname in colnames(geneExprNor)) {
  name_list <- oneColname %>% substring(14,15) %>% as.numeric()
  if(name_list>=1&name_list<=9) {
    if (substr(oneColname,16,16)=="A") {
      tumorSample = c(tumorSample,oneColname)
    }
  }
}
geneExprNor = geneExprNor[,tumorSample]


# 挑选差异明显的的基因
#方法1
{
  geneExpr = fread(GeneExpr_data_filename)
  geneExpr = geneExpr[,-1]
  geneExpr = column_to_rownames(geneExpr,var = "V1")
  
  #取出高低风险组的样本名，EC2分为高风险组，EC1和EC3为低风险组
  highRiskSample = names(cluster_Results[[1]])[cluster_Results[[1]] == 2] %>% as.vector(.)
  lowRiskSample = names(cluster_Results[[1]])[cluster_Results[[1]] != 2] %>% as.vector(.)
  
  match(substr(highRiskSample,1,16),colnames(geneExpr)) %>% na.omit(.) %>% geneExpr[,.] -> highRiskGeneExpr
  match(substr(lowRiskSample,1,16),colnames(geneExpr)) %>% na.omit(.) %>% geneExpr[,.] -> lowRiskGeneExpr
  highRiskSample <- colnames(highRiskGeneExpr)
  lowRiskSample <- colnames(lowRiskGeneExpr)
  
  
  sampleType <- data.frame( type = c( rep("lowRisk",length(lowRiskSample)),rep("highRisk",length(highRiskSample)) ) )
  rownames(sampleType) <- c(lowRiskSample,highRiskSample)
  dds <- DESeqDataSetFromMatrix(countData = cbind(lowRiskGeneExpr,highRiskGeneExpr), colData = sampleType, design = ~ type)
  # 过滤这个效果没有提升，暂时不要了
  # keep <- rowSums(counts(dds) >= 10) >= 3   #过滤低表达基因，至少在3个样品中都有10个reads 
  # dds <- dds[keep, ] 
  dds<-DESeq(dds)
  res<-results(dds)
  res<-as.data.frame(res)
  resSig<-res[which(res$pvalue<0.05),] %>% rownames(.)
  
  match(substr(resSig,1,15),rownames(geneExprNor)) %>% geneExprNor[.,] %>% na.omit(.) -> geneExprNor
}



#方法2
{
  #匹配一下临床数据，取出相同的sample
  clinical_data_full <- fread(clinical_data_file_name) %>% unique(.)
  clinical_selected <- clinical_data_full[clinical_data_full$cancer == cancerName,]
  clinical_selected <- clinical_selected[match(unique(substr(colnames(geneExprNor),1,12)),clinical_selected$patient_id)]
  clinical_selected$os_time = round(clinical_selected$os_time/365,3)
  geneExprNor <- geneExprNor[match(clinical_selected$patient_id,unique(substr(colnames(geneExprNor),1,12)))]

  geneExprNor = log(geneExprNor+1) 
  cox <- apply(
    geneExprNor,1,function(x){
      clinical_selected$gene <- as.numeric(x)
      cox_genes <- coxph(Surv(os_time, os) ~ gene, data = clinical_selected)
      coef <- coef(cox_genes) #回归系数
      SE <- sqrt(diag(vcov(cox_genes))) #标准误
      HR <- exp(coef) #风险比
      cox_need <- cbind(HR = HR,
                        HR.95L = exp(coef - qnorm(.975, 0, 1) * SE),
                        HR.95H = exp(coef + qnorm(.975, 0, 1) * SE),
                        pvalue = 1 - pchisq((coef/SE)^2, 1))
      return(cox_need['gene',])
    }
  )
  unicox <- t(cox)
  
  select_gene_name = rownames(unicox)[ unicox[,4]<0.05 ] %>% na.omit(.)
}




#a）LASSO模型中基于交叉验证进行变量选择; (在运行时自动生成)

#匹配一下临床数据，取出相同的sample
clinical_data_full <- fread(clinical_data_file_name) %>% unique(.)
clinical_selected <- clinical_data_full[clinical_data_full$cancer == cancerName,]
clinical_selected <- clinical_selected[match(unique(substr(colnames(geneExprNor),1,12)),clinical_selected$patient_id)]

#取出y的值，即os和ostime
matchCol = match(unique(clinical_selected$patient_id),substr(colnames(geneExprNor),1,12)) %>% na.omit(.)
survival_data = geneExprNor[,matchCol]  %>% t(.)
dataSurvivalRowName <- rownames(survival_data)
survival_data = cbind(dataSurvivalRowName,survival_data,clinical_selected[match(substr(dataSurvivalRowName,1,12),clinical_selected$patient_id),8:9]) %>% na.omit(.)

#要删除存活时间为0的
survival_data = survival_data[-which(survival_data$os_time == 0),]
survival_data_selt = survival_data[,(ncol(survival_data)-1):ncol(survival_data)] %>% set_colnames(c("status","time")) %>% remove_rownames() %>% as.matrix()
survival_data_selt = survival_data_selt[,c(2,1)]

#方法1
{
  xMatrix = geneExprNor[,match(survival_data$dataSurvivalRowName, colnames(geneExprNor))] %>% t(.) %>% as.matrix(.)
  xMatrix = log(xMatrix+1)
}

#方法2
{
  xMatrix = geneExprNor[select_gene_name,match(survival_data$dataSurvivalRowName, colnames(geneExprNor))] %>% t(.) %>% as.matrix(.)
}


#保存一下所有数据，不用重新跑
save.image(file = paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,".Step_lasso.RData",sep = ''))

cancerName = "LIHC"
load(file = paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,".Step_lasso.RData",sep = ''))
source("/home/bioinfo/TCGA/pancancer/Script/6.func.R")
source('/home/bioinfo/TCGA/pancancer/Script/dca.R')
setwd(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/",cancerName,sep = ''))




########################################################
set.seed(600)
cvfit <<- cv.glmnet(xMatrix, survival_data_selt, family = "cox")
pdf("3.4.1.A.Classifier.cvfit.pdf", onefile=F)
plot(cvfit)
dev.off()
## 随着lambda值的变化，每个观察值对应的系数的变化趋势。
fitModel <- glmnet(xMatrix, as.matrix(survival_data_selt), family = "cox", alpha = 1)
pdf("3.4.1.B.Classifier.fit.pdf", onefile=F)
plot(fitModel)
dev.off()
##
coef.min = coef(cvfit, s = "lambda.min")
active.min = which(as.matrix(coef.min) != 0)
index.min = coef.min[active.min]
each_model = data.frame( gene=colnames(xMatrix)[active.min], coef = index.min)
singifPosiData = xMatrix[,each_model$gene] %>% as.data.frame(.) 
rownames(singifPosiData) <- survival_data$dataSurvivalRowName
singifPosiData = rbind(singifPosiData, each_model$coef)
write.csv(singifPosiData, paste("3.4.1.A.Classifier.",cvfit["lambda.min"][[1]],".csv",sep = ''), row.names = F)
print(each_model)

## 保存基因和稀疏，先进行基因名转换
# geneName=bitr(each_model$ENSEMBL,fromType="ENSEMBL" ,toType="SYMBOL",OrgDb="org.Hs.eg.db")
g2e=unique(toTable(org.Hs.egENSEMBL))
g2s=unique(toTable(org.Hs.egSYMBOL))

symbolIds = c()
record_row = c()
for (i in 1:nrow(each_model)) {
  if (length(g2s$symbol[g2s$gene_id == g2e$gene_id[g2e$ensembl_id == each_model$gene[i]]]) != 0) {
    symbolIds <- c(symbolIds, g2s$symbol[g2s$gene_id == g2e$gene_id[g2e$ensembl_id == each_model$gene[i]]])
    record_row = c(record_row,i)
  }
  cat(paste(g2s$symbol[g2s$gene_id == g2e$gene_id[g2e$ensembl_id == each_model$gene[i]]],each_model$gene[i],"\n"))
}
selectGeneData = cbind(each_model[record_row,], symbol = symbolIds)
write.csv(selectGeneData, "3.4.1.A.Classifier.selectGeneData.csv",row.names = F)




# d)对应的是下图中的d。用survminer包的函数surv_cutpoint来找到risk score的最佳cutoff值，
# 不同cutoff时标准化的log-rank也会做一个统计量分布图，虚线为选择的最佳分割点;
riskScore = 0
for (Ind in 1:nrow(xMatrix)) {
  riskScore[Ind] = sum(index.min * xMatrix[Ind,active.min])
}

survival_data_selt = as.data.frame(survival_data_selt) %>% set_rownames(survival_data$dataSurvivalRowName)
y_ori = survival_data_selt
y_ori$riskScore = riskScore

res.cut = survminer::surv_cutpoint(y_ori, time = "time", event = "status",minprop = 0.1, variables=c("riskScore"),progressbar = F)
cut = as.data.frame(summary(res.cut))["riskScore","cutpoint"]
pdf("3.4.1.d.ori.riskScore_cutPoint.pdf", height=6.5) 
cut1 = as.data.frame(summary(res.cut))["riskScore","cutpoint"]
plot(res.cut$riskScore,pch=20,cex=1.4,xlab = "",ylab = "Standardized log-rank statistic") 
text(x=cut1,y=res.cut$riskScore$stats %>% min(.) + 0.1,pos=4, labels = round(cut1,3))
dev.off()
y_ori$group=ifelse(y_ori$riskScore < cut1,"risk low","risk high")
# e)预后模型评分在病人间的分布（每个病人的risk score是多少，其中高组和低组用不同的颜色表示）也要做一个图；

plotRiskScore(y_ori,"ori",cancerName)

#根据高低组画一个基因表达量热图
geneExpr = fread( GeneExpr_fpkm_filename ) %>% as.data.frame(.)
colnames(geneExpr) <- gsub("[.]","-",colnames(geneExpr))
geneExpr = geneExpr[,c(1,match(rownames(y_ori)[order(y_ori$riskScore)], colnames(geneExpr)))]
each_model_heat = selectGeneData
colnames(each_model_heat) <- c("gene_id","coef", "symbol")
each_model_heat = data.frame( gene_id = as.character(each_model_heat$gene_id), symbol = as.character(each_model_heat$symbol) )
geneExprHeatPlot = merge(each_model_heat, geneExpr, by = "gene_id")
annotationCol = data.frame( ECBS = y_ori$riskScore , patient = rownames(y_ori))
annotationCol = annotationCol[order(annotationCol$ECBS),]
rownames(annotationCol) = annotationCol$patient
annotationCol = annotationCol["ECBS"]

geneExprHeatPlot = column_to_rownames(geneExprHeatPlot, var = "symbol")
geneExprHeatPlot = geneExprHeatPlot[,-1]


bk <- c(seq(-6,-1,by=0.1),seq(-0.9,-0.1,by=0.01),seq(0,0.9,by=0.01),seq(1,6,by=0.1))
pheatmap(log10(geneExprHeatPlot+1), 
         color = c(colorRampPalette(colors = c("steelblue","white"))(length(bk)/2),colorRampPalette(colors = c("white","tomato"))(length(bk)/2)),
         breaks=bk,
         # color = colorRampPalette(c("steelblue", "white","tomato"))(50),
         cluster_rows = F, cluster_cols = F, 
         clustering_method="ward.D2",
         # gaps_row = c(25),
         # clustering_distance_rows = "correlation",
         # clustering_distance_rows = drows,
         scale = "row",
         annotation_col = annotationCol,
         #         annotation_colors=ann_colors,
         show_rownames = T,
         show_colnames = F,
         # main=title_,
         filename = "3.4.1.G.lasso.gene.heatmap.pdf",
         # treeheight_row=15,
         fontsize_row=10,
         cellwidth = 1,
         # cellwidth = 2,cellheight=6,
         # cutree_row = 2,
         cutree_cols=50,
         border=T, legend = T)




#############################
#2）两个预后风险评估模型之间的对比分析


# 画ROC曲线
colorName = c("#1E90FF","#FF4500")
pdf("3.4.2.A.ori.ROC.pdf") 
y_new = y_ori
y_new$time = y_new$time/30
each_yr_roc(y_new,colorName)
dev.off()

# KM生存分析
survival_Result = survfit(Surv(time, status) ~ group, y_ori)
surv_By_RiskScore(survival_Result)


data.survdiff <- survdiff(Surv(time, status) ~ group,y_ori)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1])) 
cat(paste("HR:",HR,"up95:",up95,"low95",low95))

#3)免疫浸润评分的热图：高低风险组之间在免疫评分上有什么整体区别。这里跟第1-2)-C步骤相似，做高低风险组间免疫浸润评分的热图。

immune_Data = fread(Immune_data_file_name)
immune_Score_Heatmap(y_ori,immune_Data,"3.4.2.C.immune.pdf",cancerName)

#热图不好画小提琴

y_violin = y_ori
y_violin$group[which(y_violin$group =='risk high')] <- 'Risk-high'
y_violin$group[which(y_violin$group =='risk low')] <- 'Risk-low'


groupLabel = y_violin['group']
groupLabelData = rownames_to_column(groupLabel)
groupLabelData$rowname = substr(groupLabelData$rowname,1,12)

immune_Data = fread(Immune_data_file_name)
# sigpara = c("B cells memory","T cells CD4 memory resting",
#             "T cells CD4 memory activated","T cells CD4 memory activated",
#             "NK cells resting","Monocytes","Monocytes",
#             "Mast cells resting","Eosinophils")
violinPlotDataAll = c()
for (i in 2:23) {
  
  violinPlotData = data.frame(rowname = immune_Data$patient_id,value = immune_Data[,..i])
  violinPlotData = merge(violinPlotData , groupLabelData, by="rowname")
  violinPlotData$rowname = colnames(immune_Data)[i]
  colnames(violinPlotData) <- c("types","value","group")
  violinPlotDataAll = rbind(violinPlotDataAll ,  violinPlotData)
  
}

colnames(violinPlotDataAll) <- c("CellType","Composition","group")

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
  stat_compare_means( aes(group = group),label = "p.format",hide.ns =T) +
  theme(axis.text.x = element_text( angle = 90, hjust = 1, vjust = 0.5 )) 

ggsave("3.4.2.C.immune.violin.pdf", p,width = 20,heigh=7)





#4)单因素Cox+多因素Cox 验证预后/风险评估模型独立性

#再匹配一下临床数据，取出相同的sample
clinical_data_full <- fread(clinical_data_file_name)  %>% unique(.)
clinical <- clinical_data_full[clinical_data_full$cancer == cancerName,]
clinical <- clinical[match(unique(substr(rownames(survival_data_selt),1,12)),clinical$patient_id)]

#处理一下stage部分的数据
clinical$ajcc_pathologic_tumor_stage[grep("Stage IA",clinical$ajcc_pathologic_tumor_stage)] <- "Stage I"
clinical$ajcc_pathologic_tumor_stage[grep("Stage IB",clinical$ajcc_pathologic_tumor_stage)] <- "Stage I"
clinical$ajcc_pathologic_tumor_stage[grep("Stage IIA",clinical$ajcc_pathologic_tumor_stage)] <- "Stage II"
clinical$ajcc_pathologic_tumor_stage[grep("Stage IIB",clinical$ajcc_pathologic_tumor_stage)] <- "Stage I"
clinical$ajcc_pathologic_tumor_stage[grep("III",clinical$ajcc_pathologic_tumor_stage)] <- "Stage III"
clinical$ajcc_pathologic_tumor_stage[grep("IV",clinical$ajcc_pathologic_tumor_stage)] <- "Stage IV"

#A)先做批量单因素COX回归，看：risk score、age、gender、stage分别与生存之间是否相关。
res1 = unit_Cox_reg(y_ori,clinical,c("group","age", "gender"))
res2 = unit_Cox_reg(y_ori,clinical,c("stage"))
result = rbind(res1,res2)
result$name = str_remove(rownames(result),"stage")
result = result[,c(6,1:5)]
result = rbind(c("","","","","HR (95% CI)","P value"),result) 
result$name = c("","group (low vs. high)","age","gender (male vs.female)","TNM (II vs. I)","TNM (III vs. I)","TNM (IV vs. I)")
pdf("3.4.3.A.forestPlot.pdf", onefile=F, height = 3)
forestplot(result[,c(1,5,6)], mean=result[,2], lower=result[,3], upper=result[,4], zero=1, boxsize=0.4, graph.pos= "right",
           hrzl_lines=list("1" = gpar(lty=1,lwd=2), "2" = gpar(lty=1,columns=c(2:3)), "8"= gpar(lwd=2,lty=1,columns=c(1:3)) ),
           graphwidth = unit(.3,"npc"), xticks=c(0.25,0.5,1,1.5,2.25,4), is.summary=c(T,F,F,F,F,F,F,F),
           txt_gp=fpTxtGp( label=gpar(cex=1), ticks=gpar(cex=1), xlab=gpar(cex=1.5), title=gpar(cex=1)),
           lwd.zero=1, lwd.ci=1.5, lwd.xaxis=1.5, lty.ci=1.5, ci.vertices =T, ci.vertices.height=0.2, 
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
           clip=c(0.1,3), lineheight = unit(8, 'mm'), line.margin=unit(8, 'mm'), colgap=unit(4, 'mm'),
           title=paste(cancerName,"unicariate cox regression") )
dev.off()

#B)再做多变量Cox比例风险回归，看校正完三个临床信息后，risk score是否是独立影响生存的因素。主要是HR、Pvalue，统计这些数据，并做回归统计图。
#1.前面已经载入多因素cox的函数，直接调用
pdf("3.4.3.B.forestPlot.pdf", height = 3)
clinicalMul = clinical[!(clinical$ajcc_pathologic_tumor_stage == "[Discrepancy]" | clinical$ajcc_pathologic_tumor_stage=="[Not Available]")]
multi_Cox_reg(y_ori,clinicalMul)
dev.off()


saveRDS(y_ori, file = "/home/bioinfo/TCGA/pancancer/4.model_data/results/LIHC/ECBS_Model.rds")
#取出高低风险组的样本名，EC2分为高风险组，EC1和EC3为低风险组
highRiskSample = names(cluster_Results[[1]])[cluster_Results[[1]] == 2] %>% as.vector(.)
lowRiskSample = names(cluster_Results[[1]])[cluster_Results[[1]] != 2] %>% as.vector(.)

EC_model = rbind(data.frame(sample = highRiskSample, group = "risk high"),
                 data.frame(sample = lowRiskSample, group = "risk low"))
saveRDS(EC_model, file = "/home/bioinfo/TCGA/pancancer/4.model_data/results/LIHC/EC_Model.rds")

######
#保存一下所有数据，不用重新跑
save.image(file = paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,".Step2.RData",sep = ''))
cancerName = "LIHC"
load(file = paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,".Step2.RData",sep = ''))



#5)与其他一些常见的预后模型作对比，评估优势。只需看ROC曲线的线下面积（AUC）就行，常见预后风险模型可以参考黄璐用的那几个，另外再加上
geneExpr = fread(GeneExpr_fpkm_filename)
colnames(geneExpr) <- gsub("[.]","-",colnames(geneExpr))
geneExpr = t(geneExpr)
colnames(geneExpr) = geneExpr[1,]
geneExpr = geneExpr[-1,]
match(rownames(y_ori),rownames(geneExpr)) %>% na.omit(.) -> matchCol
geneExprLIHC = geneExpr[matchCol,]
rownames(geneExprLIHC) <- rownames(geneExpr)[matchCol]
y_ori = y_ori[match(rownames(geneExprLIHC),rownames(y_ori)) %>% na.omit(.),]
geneExprLIHC = apply(geneExprLIHC,2,as.numeric)
rownames(geneExprLIHC) <- rownames(geneExpr)[matchCol]
geneExprLIHC = as.data.frame(geneExprLIHC,stringsAsFactors = F)
y_ori_1 = y_ori
y_ori_1$time = y_ori_1$time/30
geneExprLIHC$TIME_TO_EVENT = y_ori_1$time
geneExprLIHC$EVENT = y_ori_1$status
my_surv_df = y_ori_1[,1:3] %>% set_colnames(value = c("TIME_TO_EVENT","EVENT","risk"))

# compOtherModel函数包括计算和画图，比较其他模型
# 第一个输入是基因表达矩阵，第二个输入是自己的模型
model7_list = compOtherModel(geneExprLIHC, my_surv_df, cancerName, "3.4.4.compareWithOthers.pdf")     
summary(model7_list,times=seq(12,36,12))


plot(model7_list, 
     col=c("#E41A1C","#4DAF4A","#1E90FF","#FF8C00","#984EA3","#EB6046","#AFEEEE","#A9A9A9"), 
     xlim=c(0,60),lagend=F)
write.table(summary(model7_list,times=seq(12,36,12)),"3.4.4.table.txt")
legend(x = "bottomright",          # Position
       legend = c("Reference",
                  "ECBS", 
                  "Yuan et al", 
                  "Kim et al", 
                  "Jiang et al",
                  "Kong et al",
                  "Liu et al",
                  "Huang et al"),  # Legend texts
       lty = c(1,1,1,1,1,1,1),           # Line types
       col = c("#E41A1C","#4DAF4A","#1E90FF","#FF8C00","#984EA3","#EB6046","#AFEEEE","#A9A9A9"),           # Line colors
       lwd = 2,
       cex = 1.2)                 # Line width
dev.off()


#5预后列线图的构建与评估
# deciscion curves 
aaa = cbind(y_ori,clinical[,6])
aaa$time = aaa$time/30
colnames(aaa) <- c("time","status","risk Score","ECBS.Group","TNM.Stage")
for (i in c(1,2)) { aaa[,i] = aaa[,i] %>% as.numeric() }
aaa = aaa[!(aaa$TNM.Stage == "[Discrepancy]" | aaa$TNM.Stage =="[Not Available]" | aaa$TNM.Stage =="[Discrepancy]"),]
dd <- datadist(aaa)
options(datadist="dd")

coxm <-cph(Surv(time,status)~ECBS.Group+TNM.Stage,x=T,y=T,data=aaa,surv=T)
surv<- Survival(coxm) # 建立生存函数
surv1<- function(x)surv(12,lp=x) # 定义time.inc,1年OS
surv2<- function(x)surv(24,lp=x) # 定义time.inc,2年OS
surv3<- function(x)surv(36,lp=x) # 定义time.inc,3年OS
pdf("3.5.1.A.liexiantu.pdf", width = 9, height = 6)
plot(nomogram(coxm,fun=list(surv1,surv2,surv3),lp=F,
              funlabel=c('12 month OS','24 month OS','36 month OS'),
              maxscale=100,
              est.all=F,
              fun.at=c('0.95','0.9','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1','0.05')),xfrac=.35,
     lplabel="Linear Predictor", varname.label = T , varname.label.sep= "=", 
     ia.space=.2, tck = NA, tcl = -0.20, lmgp = 0.3,
     points.label = 'Points', total.points.label = 'Total Points', 
     total.sep.page = F, cap.labels = F, cex.var = 1.6, 
     cex.axis = 1.05, lwd = 5, label.every = 1, col.grid = gray(c(0.8, 0.95)))
dev.off()

# calibration curves
colName = c("#E41A1C","#4DAF4A","#1E90FF")
my_formula="Surv(time,status) ~ ECBS.Group + TNM.Stage"
pdf("3.5.1.B.deciscionCurves.pdf")
# plot_nomogram_cali_curve(my_formula,aaa,60,12,colName[1]); par(new=TRUE)
plot_nomogram_cali_curve(my_formula,aaa,60,24,colName[2]); par(new=TRUE)
plot_nomogram_cali_curve_re(my_formula,aaa,60,32,colName[3])
legend(x = "bottomright",          # Position
       legend = c("24 Months","36 Months"),  # Legend texts
       lty = c(1,1),           # Line types
       pch = c(16,16),
       col = colName[2:3],           # Line colors
       lwd = 3, cex = 1.4)                 # Line width
dev.off()
coxm = coxph(as.formula(my_formula),data=aaa)
sum.surv<-summary(coxm)



# deciscion curves 
nomo_score=predict(plot_nomogram(aaa)$mymodel) %>% as.data.frame() %>% set_colnames("nomo_score")
nomodata.nomo_score=merge(aaa,nomo_score,by=0) %>% column_to_rownames(var="Row.names")
nomodata.nomo_score = cbind(aaa,nomo_score)



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

colnames(nomodata.dca) <- c("TIME_TO_EVENT","EVENT","risk Score","ECBS Group","TNM","nomo_score")

pdf("3.5.1.C.12month.decisionCurces.pdf")
result1=stdca(data=nomodata.dca, outcome="EVENT", ttoutcome="TIME_TO_EVENT",
              timepoint=12, predictors=c("TNM","ECBS Group","nomo_score"),
              probability=c(F,F,F), xstop=0.7)
dev.off()
pdf("3.5.1.C.24month.decisionCurces.pdf")
result2=stdca(data=nomodata.dca, outcome="EVENT", ttoutcome="TIME_TO_EVENT",
              timepoint=24, predictors=c("TNM","ECBS Group","nomo_score"), 
              probability=c(F,F,F), xstop=0.8)
dev.off()
pdf("3.5.1.C.36month.decisionCurces.pdf")
result3=stdca(data=nomodata.dca, outcome="EVENT", ttoutcome="TIME_TO_EVENT",
              timepoint=36, predictors=c("TNM","ECBS Group","nomo_score"), 
              probability=c(F,F,F), xstop=0.5)
dev.off()

#2）另外，这里还要再加一组ROC曲线，分别用TCGA、GSE14520、和LIRI去做，对比24个月和36个月的：列线图、EGBS模型、TNM分期 的性能。。
# 画ROC曲线
colorName = c("#4DAF4A","#E41A1C","#1E90FF")
y_roc = nomodata.dca
filename = "/home/bioinfo/TCGA/pancancer/4.model_data/results/LIHC"
each_yr_roc_mul(y_roc,colorName,cancerName,filename)






#6. 肝细胞癌Driver RNA编辑位点的功能分析
{
  
  
  
  #1)对于未压缩模型中做出贡献的位点，看这些位点对mRNA、miRNA之类的影响，以及其host gene的功能。
  
  #A)对mRNA的影响：看这些位点与peptide那一部分的overlap，看有哪些位点是可以导致peptide产生非同义突变的；
  
  paste(matchData$Chr,matchData$`Position(1base)`,sep = '')[active.min_ori] -> nonCompPosit
  peptideData = fread("/home/bioinfo/TCGA/pancancer/CancerPeptide/positionAndCancer.csv")
  LIHCPeptide = peptideData[peptideData$cancers==cancerName]
  overlapPosi <- match(nonCompPosit,LIHCPeptide$position)  %>% na.omit(.) %>% LIHCPeptide$position[.]
  
  
  #B)对miRNA的影响：对miRNA本身或者其靶点的影响，用miRanda软件，预测miRNA与mRNA结合力来预测编辑区域和对照区域的miRNA结合情况，统计editing loss/gain/change。
  
  
  
  #C)根据A中位点的坐标（chr+end/position）在EF_informative.csv找到这些位点，每个位点对应的第11列就是它的host gene name，也就是这个位点所在的基因名。
  library(org.Hs.eg.db)
  edInfo <- fread(paste("/home/bioinfo/TCGA/pancancer/informative_RES/",cancerName,"_EF_informative_no_SNP_2.csv",sep = ""))
  
  geneExpr = fread(paste("/home/bioinfo/FPKM_cli/fpkm/","TCGA-",cancerName,"_htseq_fpkm.tsv",sep = ''))
  hostGene <- match(nonCompPosit,paste(edInfo$Chr,edInfo$`Position(1base)`,sep = '')) %>% edInfo$Gene.wgEncodeGencodeBasicV34[.]
  positions <- match(nonCompPosit,paste(edInfo$Chr,edInfo$`Position(1base)`,sep = '')) %>% paste(edInfo$Chr,edInfo$`Position(1base)`,sep = '')[.]
  geneAndPost = cbind(hostGene,positions) %>% na.omit(.) %>% as.data.frame(.)
  colnames(geneAndPost) <- c("gene","position")
  write.csv(geneAndPost,file = paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,"/6.1.C.ori.hostGene.csv",sep = ''),row.names = F)
  
  #D)host gene 的功能分析：根据C中找出的host gene,看高低风险组之间这些基因的基因表达量有无显著性差异。找出高低组中对应基因的表达量，用热图来对照差异。
  
  
  g2e=unique(toTable(org.Hs.egENSEMBL)) 
  g2s=unique(toTable(org.Hs.egSYMBOL))
  
  exprResult = data.frame()
  for (i in 1:nrow(geneAndPost)) {
    geneName = geneAndPost$gene[i] %>% as.character()
    ensemblId1 <- g2e$ensembl_id[g2e$gene_id == g2s$gene_id[g2s$symbol == geneName]] 
    if(isEmpty(ensemblId1)) {
      print(geneName)
      next
    }else if(length(ensemblId1)>1){
      for (eId in ensemblId1) {
        if (!isEmpty(which(geneExpr$gene_id == eId))) {
          ensemblId1 = eId
          break()
        }
      }
    }
    if(isEmpty(which(geneExpr$gene_id == ensemblId1))) {
      print(geneName)
      next
    }
    exprResult = rbind(exprResult,cbind(geneAndPost[i,], log2(as.numeric(geneExpr[geneExpr$gene_id == ensemblId1,-1]) + 1) %>% t(.) ))
  }
  
  colnames(exprResult) <- c("gene","position",gsub("[.]","-",colnames(geneExpr)[-1]))
  
  xData = exprResult[,c(match(rownames(y_ori)[y_ori$group=="risk high"],colnames(exprResult)),
                        match(rownames(y_ori)[y_ori$group=="risk low"],colnames(exprResult)))] %>% as.matrix(.)
  rownames(xData) <- paste(exprResult$position,": ",exprResult$gene,sep = '')
  
  labX  =  rep("risk high",match(rownames(y_ori)[y_ori$group=="risk high"],colnames(exprResult)) %>% length()) %>% 
    c(.,rep("risk low",match(rownames(y_ori)[y_ori$group=="risk low"],colnames(exprResult)) %>% length())) %>%
    t(.) %>% as.factor(.) 
  sampleClass = data.frame(group = labX)
  rownames(sampleClass) <- colnames(xData)
  
  fileName = paste("/home/bioinfo/TCGA/pancancer/04_/",cancerName,"/6.1.D.","ori",".heatmap.pdf",sep = '')
  pheatmap(xData,show_colnames = F,filename = fileName,fontsize_row=5, annotation_col =sampleClass,cluster_cols = F)
  
  
}














##########################################################
#GSVA
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(enrichplot)
gsvaData = fread("/home/bioinfo/TCGA/pancancer/04_/GSVA/rawdata/HTSeq-counts/TCGA-LIHC_counts.tsv")
colnames(gsvaData) <- gsub("[.]", "-", colnames(gsvaData))

match(rownames(y_ori[y_ori$group=="risk high",]),colnames(gsvaData)) -> matchColHigh
match(rownames(y_ori[y_ori$group=="risk low",]),colnames(gsvaData)) -> matchColLow
gsvaDataHigh = gsvaData[,..matchColHigh]
gsvaDataLow = gsvaData[,..matchColLow]

apply(gsvaDataHigh,1,mean) -> sumHigh
apply(gsvaDataLow,1,mean) -> sumLow
foldChange = log2(sumHigh/(sumLow))

gene = gsvaData$gene_id %>% as.data.frame() 
colnames(gene) <- "gene_id"
gene$logFC = foldChange

geneName<-str_trim(gene$gene_id,"both") #定义gene
#开始ID转换
geneName=bitr(geneName,fromType="ENSEMBL" ,toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
## 去重
geneName <- dplyr::distinct(geneName,ENSEMBL,.keep_all=TRUE)
gene = gene[match(geneName$ENSEMBL,gene$gene_id),]
geneList <- gene$logFC
names(geneList) <- as.character(geneName$ENTREZID)
geneList <- sort(geneList,decreasing = T)




gsvaGeneSet = read.gmt("E:/bioinfo/GSVA/rawdata/HTSeq-counts/h.all.v7.4.entrez.gmt")

KEGG<-GSEA(geneList,TERM2GENE = gsvaGeneSet)
dotplot(KEGG) #出点图 
dotplot(KEGG,color="pvalue")  #按p值出点图 
dotplot(KEGG,split=".sign")+facet_grid(~.sign) #出点图，并且分面激活和抑制

#特定通路作图
gseaplot2(KEGG,1,color="red",pvalue_table = T) # 按第一个做二维码图，并显示p值
gseaplot2(KEGG,1:10,color="red") #按第一到第十个出图，不显示p值
########################################################################################################################
