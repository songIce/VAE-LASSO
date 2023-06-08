
#单因素cos回归函数
#下面的代码是抄 https://zhuanlan.zhihu.com/p/339754573 的
unit_Cox_reg <- function(yData,cliData,covariates){
  #先把riskScore放进临床数据中(好像可以不用做因为输入数据已经匹配好了但是保险起见还是做一下)
  cliData$riskScore = yData$riskScore[match(cliData$patient_id,substr(rownames(yData),1,12))]
  cliData$group = yData$group[match(cliData$patient_id,substr(rownames(yData),1,12))]
  colnames(cliData) <- c("patient_id", "cancer", "age", "gender", "race", "stage", "status", "os", "os_time", "pfs", "pfs_time", "riskScore","group")
  if (covariates[1] == "stage") {
    cliData = cliData[!(cliData$stage == "[Discrepancy]" | cliData$stage=="[Not Available]")]
    univ_model = coxph( Surv(os_time, os)~stage, data = cliData )
    mul_cox1 <- summary(univ_model)
    multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
    #一-3、multi2：提取：HR(95%CI)和P
    multi2<-ShowRegTable(univ_model, exp=TRUE, digits=2, pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)
    res <-cbind(multi1,multi2); #一-4.将两次提取结果合并成表；取名result
  }else{
    #分别对每一个变量，构建生存分析的公式
    univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(os_time, os)~', x)))
    #循环对每一个特征做cox回归分析
    univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cliData)})
    univ_results <- lapply(univ_models,
                           function(x){ 
                             mul_cox1 <- summary(x) # multi1：提取：变量+HR+95%CI+95%CI multi2：提取：HR(95%CI)和P
                             multi1<-round(mul_cox1$conf.int[, c(1, 3, 4)], 2)
                             multi2<-ShowRegTable(x, exp=TRUE, digits=2, pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)
                             result <-cbind(t(multi1),multi2) #一-4.将两次提取结果合并成表；取名result 
                             return(result)
                           })
    #转换成数据框，并转置
    res = data.frame()
    for (univ in univ_results) { res = rbind(res,univ) }
  }
  for(i in c(1:3)) {res[, i] = as.character(res[,i]) %>% as.numeric(.)}
  for(i in c(4,5)) {res[, i] = as.character(res[,i])}
  return(res)
}

# 用在compare癌症版本的
unit_Cox_reg_Comp <- function(cliData,covariates){
  #先把riskScore放进临床数据中(好像可以不用做因为输入数据已经匹配好了但是保险起见还是做一下)
  if (covariates[1] == "stage") {
    univ_model = coxph( Surv(time, status)~stage, data = cliData )
    mul_cox1 <- summary(univ_model)
    multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
    #一-3、multi2：提取：HR(95%CI)和P
    multi2<-ShowRegTable(univ_model, exp=TRUE, digits=2, pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)
    res <-cbind(multi1,multi2); #一-4.将两次提取结果合并成表；取名result
  }else{
    #分别对每一个变量，构建生存分析的公式
    univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(time, status)~', x)))
    #循环对每一个特征做cox回归分析
    univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cliData)})
    univ_results <- lapply(univ_models,
                           function(x){ 
                             mul_cox1 <- summary(x) # multi1：提取：变量+HR+95%CI+95%CI multi2：提取：HR(95%CI)和P
                             multi1<-round(mul_cox1$conf.int[, c(1, 3, 4)], 2)
                             multi2<-ShowRegTable(x, exp=TRUE, digits=2, pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)
                             result <-cbind(t(multi1),multi2) #一-4.将两次提取结果合并成表；取名result 
                             return(result)
                           })
    #转换成数据框，并转置
    res = data.frame()
    for (univ in univ_results) { res = rbind(res,univ) }
  }
  for(i in c(1:3)) {res[, i] = as.character(res[,i]) %>% as.numeric(.)}
  for(i in c(4,5)) {res[, i] = as.character(res[,i])}
  return(res)
}


#多因素cos回归函数
multi_Cox_reg <- function(yData,cliData){
  #先把riskScore放进临床数据中(好像可以不用做因为输入数据已经匹配好了但是保险起见还是做一下)
  cliData$riskScore = yData$riskScore[match(cliData$patient_id,substr(rownames(yData),1,12))]
  cliData$group = yData$group[match(cliData$patient_id,substr(rownames(yData),1,12))]
  colnames(cliData) <- c("patient_id", "cancer", "age", "gender", "race", "stage", "status", "os", "os_time", "pfs", "pfs_time", "riskScore", "group")
  
  #一-1.cox多因素回归分析
  mul_cox<-coxph(Surv(os_time,os)~group+age+gender+stage,data=cliData)
  #一-2 multi1：提取：变量+HR+95%CI+95%CI
  mul_cox1 <- summary(mul_cox)
  multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
  #一-3、multi2：提取：HR(95%CI)和P
  multi2<-ShowRegTable(mul_cox, exp=TRUE, digits=2, pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)
  result <-cbind(multi1,multi2); #一-4.将两次提取结果合并成表；取名result
  #一-5.行名转为表格第一列，并给予命名"Characteristics"
  result<-tibble::rownames_to_column(result, var = "Characteristics");result
  fig1<- forestplot(result[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                    mean=result[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                    lower=result[,3],  #告诉函数表格第3列为5%CI，
                    upper=result[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                    zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                    boxsize=0.3,       #设置小黑块的大小
                    graph.pos=2)       #森林图应插在图形第2列
  #1.删除部分变量名，只保留亚变量即："ER_Positive"变为"Positive"
  result$Characteristics<-str_remove(result$Characteristics,"stage")
  for(i in c(2:4)) {result[, i] = as.character(result[,i]) %>% as.numeric(.)}
  for(i in c(5,6)) {result[, i] = as.character(result[,i])}
  
  result = rbind(c("","","","","HR (95% CI)","P value"),result)
  result$Characteristics = c("","group (low vs. high)","age","gender (male vs.female)","TNM (II vs. I)","TNM (III vs. I)","TNM (IV vs. I)")
  forestplot(result[,c(1,5,6)], mean=result[,2], lower=result[,3], upper=result[,4], zero=1,            
             boxsize=0.4, graph.pos= "right" ,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=1,columns=c(2:3)),
                             "8"= gpar(lwd=2,lty=1,columns=c(1:3)) ),
             graphwidth = unit(.3,"npc"), xticks=c(0.25,0.5,1,1.5,2.25,4), is.summary=c(T,F,F,F,F,F,F,F),
             txt_gp=fpTxtGp( label=gpar(cex=1), ticks=gpar(cex=1), xlab=gpar(cex=1.5), title=gpar(cex=1)),
             lwd.zero=1, lwd.ci=1.5, lwd.xaxis=1.5,  lty.ci=1.5, ci.vertices =T, ci.vertices.height=0.2, 
             col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
             clip=c(0.1,2), lineheight = unit(8, 'mm'), line.margin=unit(8, 'mm'), colgap=unit(4, 'mm'),
             title=paste(cancerName,"multivariable cox regression"))
  
}

#compare版本的多因素cos回归函数
multi_Cox_reg_Comp <- function(cliData,ki,cancerName){
  #一-1.cox多因素回归分析
  mul_cox<-coxph(Surv(time,status)~group+age+gender+stage,data=cliData)
  #一-2 multi1：提取：变量+HR+95%CI+95%CI
  mul_cox1 <- summary(mul_cox)
  multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
  #一-3、multi2：提取：HR(95%CI)和P
  multi2<-ShowRegTable(mul_cox, exp=TRUE, digits=2, pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)
  result <-cbind(multi1,multi2); #一-4.将两次提取结果合并成表；取名result
  #一-5.行名转为表格第一列，并给予命名"Characteristics"
  result<-tibble::rownames_to_column(result, var = "Characteristics");result
  fig1<- forestplot(result[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                    mean=result[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                    lower=result[,3],  #告诉函数表格第3列为5%CI，
                    upper=result[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                    zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                    boxsize=0.3,       #设置小黑块的大小
                    graph.pos=2)       #森林图应插在图形第2列
  #1.删除部分变量名，只保留亚变量即："ER_Positive"变为"Positive"
  result$Characteristics<-str_remove(result$Characteristics,"stage")
  for(i in c(2:4)) {result[, i] = as.character(result[,i]) %>% as.numeric(.)}
  for(i in c(5,6)) {result[, i] = as.character(result[,i])}
  
  result = rbind(c("","","","","HR (95% CI)","P value"),result)
  result$Characteristics = c("","group (low vs. high)","age","gender (male vs.female)","TNM (II vs. I)","TNM (III vs. I)","TNM (IV vs. I)")
  forestplot(result[,c(1,5,6)], mean=result[,2], lower=result[,3], upper=result[,4], zero=1,            
             boxsize=0.4, graph.pos= "right" ,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=1,columns=c(2:3)),
                             "8"= gpar(lwd=2,lty=1,columns=c(1:3)) ),
             graphwidth = unit(.3,"npc"), xticks=c(0.25,0.5,1,1.5,2.25,4), is.summary=c(T,F,F,F,F,F,F,F),
             txt_gp=fpTxtGp( label=gpar(cex=1), ticks=gpar(cex=1), xlab=gpar(cex=1.5), title=gpar(cex=1)),
             lwd.zero=1, lwd.ci=1.5, lwd.xaxis=1.5,  lty.ci=1.5, ci.vertices =T, ci.vertices.height=0.2, 
             col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
             clip=c(0.25,2), lineheight = unit(8, 'mm'), line.margin=unit(8, 'mm'), colgap=unit(4, 'mm'),
             title=paste(cancerName,"multivariable cox regression"))
  
}



#compare版本的多因素cos回归函数
multi_Cox_reg_CompM <- function(cliData,ki,cancerName){
  #一-1.cox多因素回归分析
  mul_cox<-coxph(Surv(time,status)~group+age+gender+ALT+Multinodular+Cirrhosis+HBV+AFP+stage,data=cliData)
  #一-2 multi1：提取：变量+HR+95%CI+95%CI
  mul_cox1 <- summary(mul_cox)
  multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
  #一-3、multi2：提取：HR(95%CI)和P
  multi2<-ShowRegTable(mul_cox, exp=TRUE, digits=2, pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)
  result <-cbind(multi1,multi2); #一-4.将两次提取结果合并成表；取名result
  #一-5.行名转为表格第一列，并给予命名"Characteristics"
  result<-tibble::rownames_to_column(result, var = "Characteristics");result
  fig1<- forestplot(result[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                    mean=result[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                    lower=result[,3],  #告诉函数表格第3列为5%CI，
                    upper=result[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                    zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                    boxsize=0.3,       #设置小黑块的大小
                    graph.pos=2)       #森林图应插在图形第2列
  #1.删除部分变量名，只保留亚变量即："ER_Positive"变为"Positive"
  result$Characteristics<-str_remove(result$Characteristics,"stage")
  for(i in c(2:4)) {result[, i] = as.character(result[,i]) %>% as.numeric(.)}
  for(i in c(5,6)) {result[, i] = as.character(result[,i])}
  
  result = rbind(c("","","","","HR (95% CI)","P value"),result)
  result$Characteristics = c("","group (low vs. high)","age","gender (male vs.female)",
                             "ALT (low vs. high)","Multinodular (Y vs. N)",
                             "Cirrhosis (Y vs. N)","HBV (CC vs. AVR-CC)",
                             "AFP (low vs. high)","TNM (II vs. I)","TNM (III vs. I)")
  forestplot(result[,c(1,5,6)], mean=result[,2], lower=result[,3], upper=result[,4], zero=1,            
             boxsize=0.4, graph.pos= "right" ,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=1,columns=c(2:3)),
                             "12"= gpar(lwd=2,lty=1,columns=c(1:3)) ),
             graphwidth = unit(.3,"npc"), xticks=c(0.25,0.5,1,1.5,2.25,4), is.summary=c(T,F,F,F,F,F,F,F,F,F,F),
             txt_gp=fpTxtGp( label=gpar(cex=1), ticks=gpar(cex=1), xlab=gpar(cex=1.5), title=gpar(cex=1)),
             lwd.zero=1, lwd.ci=1.5, lwd.xaxis=1.5,  lty.ci=1.5, ci.vertices =T, ci.vertices.height=0.2, 
             col=fpColors(box="royalblue",line="darkblue", summary="royalblue", hrz_lines = "#444444"),
             clip=c(0.25,2), lineheight = unit(8, 'mm'), line.margin=unit(8, 'mm'), colgap=unit(4, 'mm'),
             title=paste(cancerName,"multivariable cox regression"))
  
}


#多因素cos回归函数 用在EC亚型的
multi_Cox_reg_ec <- function(yData,cliData,ki){
  #先把riskScore放进临床数据中(好像可以不用做因为输入数据已经匹配好了但是保险起见还是做一下)
  cliData$riskScore = yData$riskScore[match(cliData$patient_id,substr(yData$sample,1,12))]
  colnames(cliData) <- c("patient_id", "cancer", "age", "gender", "race", "stage", "status", "os", "os_time", "pfs", "pfs_time", "AEI","ADAR1", "riskScore")
  #一-1.cox多因素回归分析
  mul_cox<-coxph(Surv(os_time,os)~riskScore+age+gender+stage+AEI+ADAR1,data=cliData)
  #一-2 multi1：提取：变量+HR+95%CI+95%CI
  mul_cox1 <- summary(mul_cox)
  multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
  #一-3、multi2：提取：HR(95%CI)和P
  multi2<-ShowRegTable(mul_cox, exp=TRUE, digits=2, pDigits =3,printToggle = TRUE, quote=FALSE, ciFun=confint)
  result <-cbind(multi1,multi2); #一-4.将两次提取结果合并成表；取名result
  #一-5.行名转为表格第一列，并给予命名"Characteristics"
  result<-tibble::rownames_to_column(result, var = "Characteristics");result
  fig1<- forestplot(result[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                    mean=result[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                    lower=result[,3],  #告诉函数表格第3列为5%CI，
                    upper=result[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                    zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                    boxsize=0.3,       #设置小黑块的大小
                    graph.pos=2)       #森林图应插在图形第2列
  #1.删除部分变量名，只保留亚变量即："ER_Positive"变为"Positive"
  result$Characteristics<-str_remove(result$Characteristics,"stage")
  for(i in c(2:4)) {result[, i] = as.character(result[,i]) %>% as.numeric(.)}
  for(i in c(5,6)) {result[, i] = as.character(result[,i])}
  
  result = rbind(c("","","","","HR (95% CI)","P value"),result)
  result$Characteristics = c("","group (low vs. high)","age","gender (male vs.female)",
                             "TNM (II vs. I)","TNM (III vs. I)","TNM (IV vs. I)","AEI","ADAR1")
  forestplot(result[,c(1,5,6)], 
             mean=result[,2], 
             lower=result[,3], 
             upper=result[,4], 
             zero=1,            
             boxsize=0.4, 
             graph.pos= "right" ,
             hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                             "2" = gpar(lty=1,columns=c(2:3)),
                             "10"= gpar(lwd=2,lty=1,columns=c(1:3)) ),
             graphwidth = unit(.3,"npc"),
             xticks=c(0.75,1,1.25,1.5,2,2.5), 
             is.summary=c(T,F,F,F,F,F,F,F,F),
             txt_gp=fpTxtGp( label=gpar(cex=1), 
                             ticks=gpar(cex=1), 
                             xlab=gpar(cex=1.5), 
                             title=gpar(cex=1)),
             lwd.zero=1, 
             lwd.ci=1.5, 
             lwd.xaxis=1.5,  
             lty.ci=1.5, 
             ci.vertices =T, 
             ci.vertices.height=0.2, 
             col=fpColors(box="royalblue",
                          line="darkblue", 
                          summary="royalblue", 
                          hrz_lines = "#444444"),
             clip=c(0.75,2.5), 
             lineheight = unit(8, 'mm'), 
             line.margin=unit(8, 'mm'), 
             colgap=unit(4, 'mm'),
             title=paste(cancerName,"multivariable cox regression"))
  
}


plot_nomogram_cali_curve=function(my_formula,mydata,sampling_num,months,mycolor){
  # my_formula="Surv(TIME_TO_EVENT,EVENT)~risk_group+TNM+tumor__Plasma.cells+nat__Macrophages.M2+nat__Mast.cells.activated"
  coxm <-cph(as.formula(my_formula),x=T,y=T,data=mydata,surv=T)
  cal<- calibrate(coxm, cmethod='KM', method='boot', u=months, m=sampling_num,B=NROW(mydata),conf.int=T)
  plot(cal,lwd=2,lty=1, legend =T,
       errbar.col=c(rgb(0,118,192,maxColorValue=255)),
       xlim=c(0,1),ylim=c(0,1),
       xlab="Nomogram-Predicted Probability of OS",
       ylab="Actual OS", 
       cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4,
       col=c(rgb(192,98,83,maxColorValue=255)),plot=F)
  lines(cal[,c("mean.predicted","KM")],type="b",lwd=3,col=mycolor, pch=16)
  abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))
}

plot_nomogram_cali_curve_re=function(my_formula,mydata,sampling_num,months,mycolor){
  # my_formula="Surv(TIME_TO_EVENT,EVENT)~risk_group+TNM+tumor__Plasma.cells+nat__Macrophages.M2+nat__Mast.cells.activated"
  coxm <-cph(as.formula(my_formula),x=T,y=T,data=mydata,surv=T)
  cal<- calibrate(coxm, cmethod='KM', method='boot', u=months, m=sampling_num,B=NROW(mydata),conf.int=T)
  plot(cal,lwd=2,lty=1,
       errbar.col=c(rgb(0,118,192,maxColorValue=255)),
       xlim=c(0,1),ylim=c(0,1),
       xlab="",
       ylab="", 
       cex.lab=1.4, cex.axis=1.4, cex.main=1.4, cex.sub=1.4,
       col=c(rgb(192,98,83,maxColorValue=255)),plot=F)
  lines(cal[,c("mean.predicted","KM")],type="b",lwd=2,col=mycolor, pch=16)
  abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue=255)))
}

plot_nomogram=function(mydata){
  library(Hmisc); library(grid); library(lattice);library(Formula); library(ggplot2);library(rms);library(survival)
  dd=datadist(mydata) 
  options(datadist="dd") 
  
  coxm <-cph(Surv(time,status)~ECBS.Group+TNM.Stage,x=T,y=T,data=mydata,surv=T)
  surv<- Survival(coxm) # 建立生存函数
  surv1<- function(x)surv(12,lp=x) # 定义time.inc,1年OS
  surv2<- function(x)surv(24,lp=x) # 定义time.inc,2年OS
  surv3<- function(x)surv(36,lp=x) # 定义time.inc,3年OS
  # plot(nomogram(coxm,fun=list(surv2,surv3),lp=F,
  #               funlabel=c('24 month OS','36 month OS'),
  #               maxscale=100,
  #               est.all=F,
  #               # maxscale=75,
  #               fun.at=c('0.95','0.9','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1','0.05')),xfrac=.3)
  #maxscale 参数指定最高分数，一般设置为100或者10分
  #fun.at 设置生存率的刻度
  #xfrac 设置数值轴与最左边标签的距离，可以调节下数值观察下图片变化情况
  f<-coxph(Surv(time,status)~ECBS.Group+TNM.Stage,
           x=T,y=T,data=mydata)
  sum.surv<-summary(f)
  c_index<-sum.surv$concordance
  return(list(c_index=c_index,mymodel=coxm))
}



# Yuan et.al. signature 6 genes
#[2017_CLINICAL.CANCER.RESEARCH_mydata_hcc_prognosis_6gene]289.full.pdf
get_score_M1=function(df){
  # M_gene=c("AHCYL2","LAMP2","SPRY1","SERPINA7","FGGY","YBX1P4")
  M_gene=c("ENSG00000158467","ENSG00000005893","ENSG00000164056","ENSG00000123561","ENSG00000172456","ENSG00000213188")
  # id_transfer(M_gene,"ENSEMBL","SYMBOL")
  MISSGENE=setdiff(M_gene,colnames(df))
  if(length(MISSGENE)==0){
    model_gene_df=df[,c("TIME_TO_EVENT","EVENT",M_gene)]
    model_gene_df %>% 
      rownames_to_column(var="ind") %>% 
      mutate(risk=0.51*ENSG00000158467+0.54*ENSG00000005893+0.36*ENSG00000164056+0.33*ENSG00000123561+0.33*ENSG00000172456+0.18*ENSG00000213188+0.001) %>% 
      dplyr::select(c("TIME_TO_EVENT","EVENT","risk","ind")) %>% 
      column_to_rownames(var="ind")
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}



# r Kim et.al. signature 65 genes
#[2012_HEPATOLOGY_65gene_signature].pdf

get_score_M2=function(df,kim_model,normalize=F){
  MISSGENE=setdiff(kim_model$gene,colnames(df))
  
  if(length(MISSGENE)==0){
    model_gene_df=df[,c("TIME_TO_EVENT","EVENT",kim_model$gene %>% as.character())] %>% na.omit()
    model_gene_df[,"risk"]=NA 
    for (i in 1:nrow(model_gene_df)) { 
      sub_risk_info=model_gene_df[i,kim_model$gene %>% as.character()] %>% 
        t() %>% as.data.frame() %>% set_colnames("ind") %>% 
        merge(.,kim_model,by.x=0,by.y="gene") %>% mutate(sub_risk=ind*coef) 
      if(normalize==F){ 
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk) 
      }else{ 
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)/sum(sub_risk_info$ind) 
      }
      # print(sum(sub_risk_info$sub_risk))
    }
    model_gene_df %>% dplyr::select(c("TIME_TO_EVENT","EVENT","risk"))
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}

# jiang et.al. 5 gene signature
# [2019_AGING_hcc_5gene_Glyco_signature]Glycolysis gene expression profilings screen for prognostic risk signature of hepatocellular carcinoma

get_score_M3=function(df,kim_model,normalize=F){
  MISSGENE=setdiff(kim_model$gene,colnames(df))
  if(length(MISSGENE)==0){
    model_gene_df=df[,c("TIME_TO_EVENT","EVENT",kim_model$gene %>% as.character(.))] %>% na.omit()
    # model_gene_df %>% mutate(risk=0.51*AHCYL2+0.54*LAMP2+0.36*SPRY1+0.33*SERPINA7+0.33*FGGY+0.18*YBX1P4+0.001) %>% 
    #   select(c("TIME_TO_EVENT","EVENT","risk"))
    model_gene_df[,"risk"]=NA
    for (i in 1:nrow(model_gene_df)) {
      sub_risk_info=model_gene_df[i,kim_model$gene %>% as.character()] %>% t() %>% as.data.frame() %>% set_colnames("ind") %>% 
        merge(.,kim_model,by.x=0,by.y="gene") %>% mutate(sub_risk=ind*coef)
      if(normalize==F){
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)
      }else{
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)/sum(sub_risk_info$ind)
      }
      # print(sum(sub_risk_info$sub_risk))
    }
    model_gene_df %>% dplyr::select(c("TIME_TO_EVENT","EVENT","risk"))
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}



# kong 3 gene signature
get_score_M4=function(df,kong_model,normalize=F){
  MISSGENE=setdiff(kong_model$gene,colnames(df))
  
  if(length(MISSGENE)==0){
    model_gene_df=df[,c("TIME_TO_EVENT","EVENT",kong_model$gene %>% as.character())] %>% na.omit()
    # model_gene_df %>% mutate(risk=0.51*AHCYL2+0.54*LAMP2+0.36*SPRY1+0.33*SERPINA7+0.33*FGGY+0.18*YBX1P4+0.001) %>% 
    #   select(c("TIME_TO_EVENT","EVENT","risk"))
    model_gene_df[,"risk"]=NA
    for (i in 1:nrow(model_gene_df)) {
      sub_risk_info=model_gene_df[i,kong_model$gene %>% as.character()] %>% t() %>% as.data.frame() %>% set_colnames("ind") %>% 
        merge(.,kong_model,by.x=0,by.y="gene") %>% mutate(sub_risk=ind*coef)
      if(normalize==F){
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)
      }else{
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)/sum(sub_risk_info$ind)
      }
      # print(sum(sub_risk_info$sub_risk))
    }
    model_gene_df %>% dplyr::select(c("TIME_TO_EVENT","EVENT","risk"))
    
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}


# liu 7 gene signature
get_score_M5=function(df,liu_model,normalize=F){
  MISSGENE=setdiff(liu_model$gene,colnames(df))
  
  if(length(MISSGENE)==0){
    model_gene_df=df[,c("TIME_TO_EVENT","EVENT",liu_model$gene%>% as.character())] %>% na.omit()
    model_gene_df[,"risk"]=NA
    for (i in 1:nrow(model_gene_df)) {
      sub_risk_info=model_gene_df[i,liu_model$gene%>% as.character()] %>% t() %>% as.data.frame() %>% set_colnames("ind") %>% 
        merge(.,liu_model,by.x=0,by.y="gene") %>% mutate(sub_risk=ind*coef)
      if(normalize==F){
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)
      }else{
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)/sum(sub_risk_info$ind)
      }
      # print(sum(sub_risk_info$sub_risk))
    }
    model_gene_df %>% dplyr::select(c("TIME_TO_EVENT","EVENT","risk"))
    
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}

# huanglu gene signature
get_score_M6=function(df, huanglu_model, normalize=F){
  MISSGENE=setdiff(huanglu_model$gene,colnames(df))
  
  if(length(MISSGENE)==0){
    model_gene_df=df[,c("TIME_TO_EVENT","EVENT",huanglu_model$gene %>% as.character())] %>% na.omit()
    # model_gene_df %>% mutate(risk=0.51*AHCYL2+0.54*LAMP2+0.36*SPRY1+0.33*SERPINA7+0.33*FGGY+0.18*YBX1P4+0.001) %>% 
    #   select(c("TIME_TO_EVENT","EVENT","risk"))
    model_gene_df[,"risk"]=NA
    for (i in 1:nrow(model_gene_df)) {
      sub_risk_info=model_gene_df[i,huanglu_model$gene %>% as.character()] %>% t() %>% as.data.frame() %>% set_colnames("ind") %>% 
        merge(.,huanglu_model,by.x=0,by.y="gene") %>% mutate(sub_risk=ind*coef)
      if(normalize==F){
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)
      }else{
        model_gene_df[i,"risk"]=sum(sub_risk_info$sub_risk)/sum(sub_risk_info$ind)
      }
      # print(sum(sub_risk_info$sub_risk))
    }
    model_gene_df %>% dplyr::select(c("TIME_TO_EVENT","EVENT","risk"))
    
  }else{
    cat(paste0("genes lost: ",paste0(MISSGENE,collapse = ",")))
  }
}



plotRiskScore <- function(y,ki,cancerName){
  cut = which(y$group == "risk low") %>% length() 
  yData = y
  yData$ind = rownames(y)
  yData$status=factor(yData$status,levels = c(0,1))
  yData$group=factor(yData$group,levels = c("risk low","risk high"))
  
  yData1 = yData[order(yData$riskScore),]
  yData1$ind = 1:nrow(y)
  p3 = ggplot(yData1,aes(x=ind,y=riskScore,color=group)) + geom_point(shape=19,size=1.4) + theme_bw() + scale_colour_lancet() +
    geom_vline(aes(xintercept=cut), colour="#000000", linetype="dashed") +
    labs(x="Patient",y="Risk Score") + theme(axis.text.x = element_blank(),panel.grid =element_blank(),axis.ticks = element_blank())

  ggsave(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/",cancerName,"/3.4.1.e.",ki,"riskScore.pdf",sep = ''), p3)
  
  # f)生存事件和生存状态随病人预后模型评分的分布也有一个图（纵轴是“time to event（用月表示的），横轴是病人的risk score，蓝色是活着，红色是死了）。
  p2=ggplot(yData,aes(x=riskScore,y=time/30,color=group)) + geom_point(shape=19,size=1.4) + theme_bw() + scale_colour_lancet() +
    labs(x="Risk Score",y="Time to Event") + theme(axis.text.x = element_blank(),panel.grid =element_blank(),axis.ticks = element_blank())

  ggsave(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/",cancerName,"/3.4.1.f.",ki,"TimeToEvent.pdf",sep = ''), p2)
}

plotRiskScore_Comp <- function(y,ki,cancerName){
  cut = which(y$group == "risk low") %>% length() 
  yData = y
  yData$ind = rownames(y)
  yData$status=factor(yData$status,levels = c(0,1))
  yData$group=factor(yData$group,levels = c("risk low","risk high"))
  
  yData1 = yData[order(yData$riskScore),]
  yData1$ind = 1:nrow(y)
  p3 = ggplot(yData1,aes(x=ind,y=riskScore,color=group)) + geom_point(shape=19,size=1.4) + theme_bw() + scale_colour_lancet() +
    geom_vline(aes(xintercept=cut), colour="#000000", linetype="dashed") +
    labs(x="Patient",y="Risk Score") + theme(axis.text.x = element_blank(),panel.grid =element_blank(),axis.ticks = element_blank())
  p3
  ggsave(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/",cancerName,"/3.4.1.e.",ki,"riskScore.pdf",sep = ''), p3)
  
  # f)生存事件和生存状态随病人预后模型评分的分布也有一个图（纵轴是“time to event（用月表示的），横轴是病人的risk score，蓝色是活着，红色是死了）。
  p2=ggplot(yData,aes(x=riskScore,y=time,color=group)) + geom_point(shape=19,size=1.4) + theme_bw() + scale_colour_lancet() +
    labs(x="Risk Score",y="Time to Event") + theme(axis.text.x = element_blank(),panel.grid =element_blank(),axis.ticks = element_blank())
  p2
  ggsave(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/",cancerName,"/3.4.1.f.",ki,"TimeToEvent.pdf",sep = ''), p2)
}


run_time_ROC=function(mytime,df,mycolor,line_width = 2){
  colnames(df)=c("TIME_TO_EVENT","EVENT","risk","group")
  ROC_result<-timeROC(T=df$TIME_TO_EVENT, delta=df$EVENT, marker=df$risk, cause=1,weighting="marginal", times=c(mytime),iid = TRUE)
  plot(ROC_result,time=mytime,col=mycolor,lw=line_width)
  ROC_sum=list(CI=(confint(ROC_result)$CI_AUC %>% .[1,1:2] %>% unname()),AUC=unname(ROC_result$AUC[2]))
  AUC_info=paste0(format(ROC_sum$AUC,digits = 2)," ", format(ROC_sum$CI[1]/100,digits = 2),"-", format(ROC_sum$CI[2]/100,digits = 2)," ")
  AUC_info
}

each_yr_roc=function(roc_data,colorName,ifsmooth=F){
  roc_data2=roc_data %>% set_colnames(c("TIME_TO_EVENT","EVENT","risk","group"))
  AUC_24=run_time_ROC(24,roc_data2,colorName[1]); par(new=TRUE)
  AUC_36=run_time_ROC(36,roc_data2,colorName[2]); par(new=TRUE)
  legend("bottomright", legend=c(paste0("AUC 24:\t",AUC_24), paste0("AUC 36:\t",AUC_36)),
         cex=0.8, col=colorName, lwd=2, lty=1, box.lwd=1, inset=.02)
}

#不太明白为什么使用survfit的时候一直说找不到输入的数据，可能是有冲突,只好把survfit这一步拿出来
surv_By_RiskScore <- function(survivalResult){
  p <- ggsurvplot(survivalResult, pval = TRUE, conf.int = FALSE, risk.table = TRUE, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "none", # Specify median survival
                  ggtheme = theme_bw(base_size = 16), # Change ggplot2 theme
                  palette = c("#015699", "#FAC00F"), xlim = c(0,2000),
                  legend.labs = c("Risk high", "Risk Low")
  )
  pdf("3.4.2.B.survival.pdf")
  print(p)
  dev.off()
}

surv_By_RiskScore_Comp <- function(survivalResult,ki,cancerName){
  p <- ggsurvplot(survivalResult, pval = T, conf.int = FALSE, risk.table = TRUE, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "none", # Specify median survival
                  ggtheme = theme_bw(base_size = 16), # Change ggplot2 theme
                  palette = c("#015699", "#FAC00F"),
                  legend.labs = c("risk high", "risk low")
  )
  pdf(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/",cancerName,"/3.4.2.B.",ki,".survival.pdf",sep = ''), onefile=F)
  print(p)
  dev.off()
}

surv_By_RiskScore_TCGA <- function(survivalResult,cancerName){
  p <- ggsurvplot(survivalResult, pval = TRUE, conf.int = FALSE, risk.table = TRUE, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "none", # Specify median survival
                  ggtheme = theme_bw(base_size = 16), # Change ggplot2 theme
                  palette = c("#015699", "#FAC00F"),
                  legend.labs = c("Risk high", "Risk Low"),
                  title = cancerName
  )
  pdf(paste("/home/bioinfo/TCGA/pancancer/4.model_data/results/TCGA/",cancerName,".survival.pdf",sep = ''), onefile=F)
  print(p)
  dev.off()
  return(p)
}





immune_Score_Heatmap = function(yData,iData,fileName,cancerName){
  imData = iData[c(match(substr(rownames(yData),1,12)[yData$group=="risk high"],iData$patient_id),
                   match(substr(rownames(yData),1,12)[yData$group=="risk low"],iData$patient_id)),] 
  imDataPosit <- imData$patient_id
  clusterGroup = data.frame( Group = match(imDataPosit,substr(rownames(yData),1,12)) %>% yData$group[.] %>% as.factor(.) )
  rownames(clusterGroup) <- imDataPosit %>% substr(.,1,12)
  
  imData = column_to_rownames(immune_Data,"patient_id") %>% t(.)
  pheatmap(imData,show_colnames = F,filename = fileName, annotation_col = clusterGroup,cluster_cols = F)

}

compOtherModel = function(geneExprData, my_surv_df, cancerName, filename){
  
  library(prodlim)
  library(survival)
  library(pec)
  # 载入其他模型基因数据
  expr_surv_df_3sets_rds="/home/bioinfo/TCGA/pancancer/4.model_data/PEC Curve/expr_surv_df_3sets.rds"
  jiang_model_file ="/home/bioinfo/TCGA/pancancer/4.model_data/PEC Curve/jiang_6_sig.csv"
  kim_model_file   ="/home/bioinfo/TCGA/pancancer/4.model_data/PEC Curve/Kim_65_gene_sig.csv"
  kong_model_file  ="/home/bioinfo/TCGA/pancancer/4.model_data/PEC Curve/kong_3_sig.txt"
  liu_model_file   ="/home/bioinfo/TCGA/pancancer/4.model_data/PEC Curve/Liu_7_sig.txt"
  huanglu_model_file   ="/home/bioinfo/TCGA/pancancer/4.model_data/PEC Curve/huanglu_sig.csv"
  
  # 其他模型的计算函数再公共函数载入
  M1_rs=get_score_M1(geneExprData)  #使用lihc，根据M1算出risk
  kim_model=read.csv(kim_model_file)[,c(8,2)] %>% set_colnames(c("gene","coef"))
  M2_rs=get_score_M2(geneExprData,kim_model)
  jiang_model=read.csv(jiang_model_file)[,c(1,3)] %>% set_colnames(c("coef","gene"))
  M3_rs=get_score_M3(geneExprData,jiang_model)
  kong_model=read.table(kong_model_file,header = T)[,c(1,3)] %>%
    set_colnames(c("coef","gene")) %>% mutate(coef=as.numeric(coef))
  M4_rs=get_score_M4(geneExprData,kong_model)
  liu_model=read.table(liu_model_file,header = T)[,c(3,2)] %>% 
    set_colnames(c("gene","coef")) %>%
    mutate(coef=as.numeric(coef))
  M5_rs=get_score_M5(geneExprData,liu_model)
  huanglu_model=read.csv(huanglu_model_file,header = T)[,c(3,2)] %>% set_colnames(c("gene","coef"))
  M6_rs=get_score_M6(geneExprData,huanglu_model)
  
  model7_list=list(my=my_surv_df,M1=M1_rs,M2=M2_rs,M3=M3_rs,M4=M4_rs,M5=M5_rs,M6=M6_rs)
  
  
  my=model7_list$my %>% plyr::rename(c("risk"="my_risk"))
  m2=model7_list$M2 %>% plyr::rename(c("risk"="m2_risk")) %>% .[,3,drop=F]
  m3=model7_list$M3 %>% plyr::rename(c("risk"="m3_risk")) %>% .[,3,drop=F]
  m4=model7_list$M4 %>% plyr::rename(c("risk"="m4_risk")) %>% .[,3,drop=F]
  m6=model7_list$M6 %>% plyr::rename(c("risk"="m6_risk")) %>% .[,3,drop=F]

  if (!(is.null(M1_rs) | is.null(M5_rs)) ) {
    m1=model7_list$M1 %>% plyr::rename(c("risk"="m1_risk")) %>% .[,3,drop=F]
    m5=model7_list$M5 %>% plyr::rename(c("risk"="m5_risk")) %>% .[,3,drop=F]
    merged_risk=merge(my,m1,by=0) %>% 
      merge(.,m2,by.x=1,by.y=0) %>% 
      merge(.,m3,by.x=1,by.y=0) %>% 
      merge(.,m4,by.x=1,by.y=0) %>% 
      merge(.,m5,by.x=1,by.y=0) %>% 
      merge(.,m6,by.x=1,by.y=0)
    set.seed(130971)
    Models <- list("ECBS"=coxph(Surv(TIME_TO_EVENT,EVENT)~my_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Yuan et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m1_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Kim et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m2_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Jiang et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m3_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Kong et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m4_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Liu et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m5_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Huang et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m6_risk,data=merged_risk,x=TRUE,y=TRUE)
    )
  }else if(is.null(M5_rs)){
    merged_risk=merge(my,m2,by=0) %>% 
      merge(.,m3,by.x=1,by.y=0) %>% 
      merge(.,m4,by.x=1,by.y=0) %>% 
      merge(.,m6,by.x=1,by.y=0)
    set.seed(130971)
    Models <- list("ECBS"=coxph(Surv(TIME_TO_EVENT,EVENT)~my_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Kim et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m2_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Jiang et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m3_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Kong et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m4_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Huang et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m6_risk,data=merged_risk,x=TRUE,y=TRUE)
    )
  }else{
    m5=model7_list$M5 %>% plyr::rename(c("risk"="m5_risk")) %>% .[,3,drop=F]
    merged_risk=merge(my,m2,by=0) %>% 
      merge(.,m3,by.x=1,by.y=0) %>% 
      merge(.,m4,by.x=1,by.y=0) %>% 
      merge(.,m5,by.x=1,by.y=0) %>% 
      merge(.,m6,by.x=1,by.y=0)
    set.seed(130971)
    Models <- list("ECBS"=coxph(Surv(TIME_TO_EVENT,EVENT)~my_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Kim et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m2_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Jiang et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m3_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Kong et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m4_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Liu et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m5_risk,data=merged_risk,x=TRUE,y=TRUE),
                   "Huang et al"=coxph(Surv(TIME_TO_EVENT,EVENT)~m6_risk,data=merged_risk,x=TRUE,y=TRUE)
    )
  }
  
  # compute the apparent prediction error
  CV = T
  if(CV==F){
    PredError <- pec(object=Models,
                     formula=Surv(TIME_TO_EVENT,EVENT)~1,
                     data=merged_risk,
                     exact=F,
                     cens.model="marginal",
                     splitcMethod="none",
                     testTimes=c(12,24,36),
                     B=0,
                     # splitMethod="Boot632plus",
                     # B=100,
                     verbose=TRUE)
  }else{
    PredError <- pec(object=Models,
                     formula=Surv(TIME_TO_EVENT,EVENT)~1,
                     data=merged_risk,
                     exact=F,
                     cens.model="marginal",
                     # splitcMethod="none",
                     # B=0,
                     splitMethod="Boot632plus",
                     B=100,
                     testTimes=c(12,24,36),
                     verbose=TRUE)
  }
  
  print(PredError,seq(12,36,12))
  summary(PredError,times=seq(12,36,12))
  pdf(filename)
  plot(PredError, 
       col=c("#E41A1C","#4DAF4A","#1E90FF","#FF8C00","#984EA3","#EB6046","#AFEEEE","#A9A9A9"), 
       xlim=c(0,60), cex.lab=1.4, cex.axis = 1.2, cex = 1,)
  if (!(is.null(M1_rs) | is.null(M5_rs))) {
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
           lwd = 2, cex = 1.2)   
  }else if(is.null(M5_rs)){
    legend(x = "bottomright",          # Position
           legend = c("Reference",
                      "ECBS", 
                      "Kim et al", 
                      "Jiang et al",
                      "Kong et al",
                      "Huang et al"),  # Legend texts
           lty = c(1,1,1,1,1,1,1),           # Line types
           col = c("#E41A1C","#4DAF4A","#1E90FF","#FF8C00","#984EA3","#EB6046","#AFEEEE","#A9A9A9"),           # Line colors
           lwd = 2, cex = 1.2)  
  }else{
    legend(x = "bottomright",          # Position
           legend = c("Reference",
                      "ECBS", 
                      "Kim et al", 
                      "Jiang et al",
                      "Kong et al",
                      "Liu et al",
                      "Huang et al"),  # Legend texts
           lty = c(1,1,1,1,1,1,1),           # Line types
           col = c("#E41A1C","#4DAF4A","#1E90FF","#FF8C00","#984EA3","#EB6046","#AFEEEE","#A9A9A9"),           # Line colors
           lwd = 2, cex = 1.2)  
  }
  return(PredError)
}


run_time_ROC_mul=function(mytime,df,mycolor,line_width = 2){
  colnames(df)=c("TIME_TO_EVENT","EVENT","risk","group","TNM","nomp Score")
  ROC_result_risk<-timeROC(T=df$TIME_TO_EVENT, 
                           delta=df$EVENT, 
                           marker=df$group, 
                           cause=1, 
                           weighting="marginal", 
                           times=c(mytime),
                           iid = TRUE)
  ROC_sum=list(CI=(confint(ROC_result_risk)$CI_AUC %>% .[1,1:2] %>% unname()),AUC=unname(ROC_result_risk$AUC[2]))
  AUC_info_risk =paste0(format(ROC_sum$AUC,digits = 2)," ", format(ROC_sum$CI[1]/100,digits = 2),"-", format(ROC_sum$CI[2]/100,digits = 2)," ")
  AUC_info_risk
  ROC_result_TNM<-timeROC(T=df$TIME_TO_EVENT, 
                          delta=df$EVENT, 
                          marker=df$TNM, 
                          cause=1,
                          weighting="marginal", 
                          times=c(mytime),
                          iid = TRUE)
  ROC_sum=list(CI=(confint(ROC_result_TNM)$CI_AUC %>% .[1,1:2] %>% unname()),AUC=unname(ROC_result_TNM$AUC[2]))
  AUC_info_TMN =paste0(format(ROC_sum$AUC,digits = 2)," ", format(ROC_sum$CI[1]/100,digits = 2),"-", format(ROC_sum$CI[2]/100,digits = 2)," ")
  AUC_info_TMN
  ROC_result_nomo<-timeROC(T=df$TIME_TO_EVENT, 
                           delta=df$EVENT, 
                           marker=df$`nomp Score`, 
                           cause=1,
                           weighting="marginal", 
                           times=c(mytime),
                           iid = TRUE)
  ROC_sum=list(CI=(confint(ROC_result_nomo)$CI_AUC %>% .[1,1:2] %>% unname()),AUC=unname(ROC_result_nomo$AUC[2]))
  AUC_info_nomo =paste0(format(ROC_sum$AUC,digits = 2)," ", format(ROC_sum$CI[1]/100,digits = 2),"-", format(ROC_sum$CI[2]/100,digits = 2)," ")
  AUC_info_nomo
  plot(ROC_result_risk,time=mytime,col=mycolor[1],lw=line_width, title="", ) ; par(new=TRUE)
  plot(ROC_result_nomo,time=mytime,col=mycolor[2],lw=line_width, title="", ) ; par(new=TRUE)
  plot(ROC_result_TNM,time=mytime,col=mycolor[3],lw=line_width, title="", )
  legend("bottomright", legend=c(paste0("\t \t \t \t ", "AUC  95% CI"),
                                 paste0("ECBS Group:\t",AUC_info_risk),
                                 paste0("nomogram: \t \t \t",AUC_info_nomo),
                                 paste0("TNM stage: \t \t",AUC_info_TMN)),
         cex=1.2, col=c("#000000",colorName), lwd=2, lty=c(0,1,1,1), box.lwd=1, inset=.02)
}

each_yr_roc_mul=function(roc_data,colorName,cancerName,filename, ifsmooth=F){
  roc_data2=roc_data %>% set_colnames(c("TIME_TO_EVENT","EVENT","risk","group","TNM","nomp Score"))
  
  pdf(paste(filename,"/24.ROC.pdf",sep = ""),paper = "a4") 
  run_time_ROC_mul(24,roc_data2,colorName)
  title(main = paste(cancerName," at 24 months",sep = ""))
  dev.off()
  
  pdf(paste(filename,"/36.ROC.pdf",sep = "") ,paper = "a4") 
  run_time_ROC_mul(36,roc_data2,colorName)
  title(paste(cancerName," at 36 months",sep = ""))
  dev.off()
}
