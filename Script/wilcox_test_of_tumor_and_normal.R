###### wilcox test of editing level between tumor and normal samples (paired Wilcoxon test, false discovery rate [FDR] < 0.05, and mean editing level difference among comparison groups [Diff] <= 5%)

#整体步骤：
#1. 根据EF——info数据中的列名，区分出tumor和normal，两者的区别在sampleID的第14、15位，小于等于9 的是tumor，大于9的则是normal。
#2. 进行tumor和normal之间的wilcox test分析，这里两组间所有的RESs（RNA editing sites）是一致的，只是0和非0 的区别，所以每个sites都会有一个pvalue和FDR
library(dplyr)
library(data.table)
library(tidyr)

#1 这个是之前分组的方法，供参考。
all_res_files <- list.files("/home/bioinfo/TCGA/SPRINT/hyper_editing/STAD/hyper_res/")
dir.create("/home/bioinfo/TCGA/SPRINT/hyper_editing/STAD/tumor_hyper_res") 
dir.create("/home/bioinfo/TCGA/SPRINT/hyper_editing/STAD/normal_hyper_res/")
tumor_folder <- "/home/bioinfo/TCGA/SPRINT/hyper_editing/STAD/tumor_hyper_res/"
normal_folder <- "/home/bioinfo/TCGA/SPRINT/hyper_editing/STAD/normal_hyper_res/"

setwd("/home/bioinfo/TCGA/SPRINT/hyper_editing/STAD/hyper_res/")
for(f in all_res_files) {
  name_list <- f %>% substring(14,15) %>% as.numeric()
  if(name_list>=1&name_list<=9) {
    file.copy(f,tumor_folder)
  }else {
    file.copy(f,normal_folder)
  }
}

#2 计算两组间wilcox test 的。FDR
# *_select_res就是根据上一步找到的tumor和normal对应的列。
pvalue_tumor_normal = data.table()
# for (i in dim(tumor_select_res[i, ])){
#   a = wilcox.test(as.numeric(tumor_select_res[i, ]), as.numeric(normal_select_res[i, ]),  paired = F)
#   pvalue_tumor_normal_tpm = data.table(a = a$p.value, b = a$statistic)
#   pvalue_tumor_normal = rbind(pvalue_tumor_normal_tpm, pvalue_tumor_normal)
# }
for (i in (1:7257)){
  a = wilcox.test(as.numeric(tumor_select_res[i, ]), as.numeric(normal_select_res[i, ]),  paired = F)
  pvalue_tumor_normal_tpm = data.table(a = a$p.value, b = a$statistic)
  pvalue_tumor_normal = rbind(pvalue_tumor_normal_tpm, pvalue_tumor_normal)
}

setwd("/dsk2/who/zhuhm132/npc/all_res/ori_all_res/EditingRate")
write.table(pvalue_tumor_normal, file = "pvalue_tumor_normal_3.txt", col.names = T, row.names = T, quote = TRUE, append = TRUE)
pvalue_all <- cbind(inter_res[ ,1:3 ], pvalue_tumor_normal)
write.table(pvalue_all, file = "pvalue_all_3.txt", col.names = T, row.names = F, quote = TRUE, append = TRUE)

sign_pvalue <- pvalue_all[pvalue_all$a<= 0.05, ]
write.table(sign_pvalue, file = "sign_pvalue_3.txt", col.names = T, row.names = F, quote = TRUE, append = TRUE)
high_sign_pvalue <- pvalue_all[pvalue_all$a<= 0.01, ]
write.table(high_sign_pvalue, file = "high_sign_pvalue_3.txt", col.names = T, row.names = F, quote = TRUE, append = TRUE)
