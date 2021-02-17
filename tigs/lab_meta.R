library(tidyverse)
library(pheatmap)

ti1<-read.table("TI_bang/genus/TI1_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti2<-read.table("TI_bang/genus/TI2_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti3<-read.table("TI_bang/genus/TI3_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti4<-read.table("TI_bang/genus/TI4_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
#ti5<-read.table("TI_bang/genus/TI5_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti6<-read.table("TI_bang/genus/TI6_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti7<-read.table("TI_bang/genus/TI7_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti8<-read.table("TI_bang/genus/TI8_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti9<-read.table("TI_bang/genus/TI9_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti10<-read.table("TI_bang/genus/TI10_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti11<-read.table("TI_bang/genus/TI11_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti12<-read.table("TI_bang/genus/TI12_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti13<-read.table("TI_bang/genus/TI13_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti14<-read.table("TI_bang/genus/TI14_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti15<-read.table("TI_bang/genus/TI14_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
ti16<-read.table("TI_bang/genus/TI15_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")

reduce_dataframe<-function(df){
  return(df[,c(3,5)])
}
  
ti1<-reduce_dataframe(ti1)
ti2<-reduce_dataframe(ti2)
ti3<-reduce_dataframe(ti3)
ti4<-reduce_dataframe(ti4)
ti6<-reduce_dataframe(ti6)
ti7<-reduce_dataframe(ti7)
ti8<-reduce_dataframe(ti8)
ti9<-reduce_dataframe(ti9)
ti10<-reduce_dataframe(ti10)
ti11<-reduce_dataframe(ti11)
ti12<-reduce_dataframe(ti12)
ti13<-reduce_dataframe(ti13)
ti14<-reduce_dataframe(ti14)
ti15<-reduce_dataframe(ti15)
ti16<-reduce_dataframe(ti16)

ti <- list(ti1,ti2,ti3,ti4,ti6,ti7,ti8,ti9,ti10,ti11,ti12,ti13,ti14,ti15,ti16) %>% reduce(inner_join, by = "taxon_name")
ti <- data.frame(ti)
ti <- ti[c(2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
colnames(ti)<-c('taxon_name','T1','T2','T3','T4','T6','T7','T8',
                'T9','T10','T11','T12','T13','T14','T15','T16')

ti$rowsum <- rowSums(ti[,2:16])

ti_top200 <- ti[order(-ti[1:200,16]),]
rownames(ti_top200) <- ti_top200$taxon_name
pheatmap(log2(ti_top200[1:40,2:16]), #display_numbers = order_by_A[1:50,2:6],
         scale = "row", #You can try to change scale if you want
         cluster_rows=T,
         cluster_cols=F,
         cutree_rows = 5,
         number_color = "grey30", #cellwidth = 90, cellheight = 10,
         margins=c(3,25),fontsize = 10,show_colnames = TRUE,show_rownames = TRUE,main = "Top 200")

######### ######### ######### ######### ######### ######### ######### 

tii1<-read.table("TII_chen/genus/TII1_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii2<-read.table("TII_chen/genus/TII2_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii3<-read.table("TII_chen/genus/TII3_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii4<-read.table("TII_chen/genus/TII4_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii6<-read.table("TII_chen/genus/TII6_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii7<-read.table("TII_chen/genus/TII7_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii8<-read.table("TII_chen/genus/TII8_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii9<-read.table("TII_chen/genus/TII9_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii10<-read.table("TII_chen/genus/TII10_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii11<-read.table("TII_chen/genus/TII11_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii12<-read.table("TII_chen/genus/TII12_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii13<-read.table("TII_chen/genus/TII13_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii14<-read.table("TII_chen/genus/TII14_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii15<-read.table("TII_chen/genus/TII14_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tii16<-read.table("TII_chen/genus/TII15_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")

tii1<-reduce_dataframe(tii1)
tii2<-reduce_dataframe(tii2)
tii3<-reduce_dataframe(tii3)
tii4<-reduce_dataframe(tii4)
tii6<-reduce_dataframe(tii6)
tii7<-reduce_dataframe(tii7)
tii8<-reduce_dataframe(tii8)
tii9<-reduce_dataframe(tii9)
tii10<-reduce_dataframe(tii10)
tii11<-reduce_dataframe(tii11)
tii12<-reduce_dataframe(tii12)
tii13<-reduce_dataframe(tii13)
tii14<-reduce_dataframe(tii14)
tii15<-reduce_dataframe(tii15)
tii16<-reduce_dataframe(tii16)

tii <- list(tii1,tii2,tii3,tii4,tii6,tii7,tii8,tii9,tii10,tii11,tii12,tii13,tii14,tii15,tii16) %>% reduce(inner_join, by = "taxon_name")
tii <- data.frame(tii)
tii <- tii[c(2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
colnames(tii)<-c('taxon_name','T1','T2','T3','T4','T6','T7','T8',
                 'T9','T10','T11','T12','T13','T14','T15','T16')

tii$rowsum <- rowSums(tii[,2:16])
head(tii_top200)
tii_top200 <- tii[order(-tii[1:200,16]),]
rownames(tii_top200) <- tii_top200$taxon_name
pheatmap(log2(tii_top200[1:40,2:16]), #display_numbers = order_by_A[1:50,2:6],
         scale = "row", #You can try to change scale if you want
         cluster_rows=T,
         cluster_cols=F,
         cutree_rows = 5,
         number_color = "grey30", #cellwidth = 90, cellheight = 10,
         margins=c(3,25),fontsize = 10,show_colnames = TRUE,show_rownames = TRUE,main = "Top 200")

######### ######### ######### ######### ######### ######### ######### ######### 
t3_filter<-function(df){
   return((df[2:200,2:3]))
}

tiii1<-read.csv("TIII_del/genus/TIII1_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii2<-read.csv("TIII_del/genus/TIII2_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii3<-read.csv("TIII_del/genus/TIII3_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii4<-read.csv("TIII_del/genus/TIII4_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii6<-read.csv("TIII_del/genus/TIII6_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii7<-read.csv("TIII_del/genus/TIII7_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii8<-read.csv("TIII_del/genus/TIII8_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii9<-read.csv("TIII_del/genus/TIII9_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii10<-read.csv("TIII_del/genus/TIII10_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii11<-read.csv("TIII_del/genus/TIII11_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii12<-read.csv("TIII_del/genus/TIII12_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii13<-read.csv("TIII_del/genus/TIII13_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii14<-read.csv("TIII_del/genus/TIII14_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii15<-read.csv("TIII_del/genus/TIII14_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiii16<-read.csv("TIII_del/genus/TIII15_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")

tiii1<-t3_filter(tiii1)
tiii2<-t3_filter(tiii2)
tiii3<-t3_filter(tiii3)
tiii4<-t3_filter(tiii4)
tiii6<-t3_filter(tiii6)
tiii7<-t3_filter(tiii7)
tiii8<-t3_filter(tiii8)
tiii9<-t3_filter(tiii9)
tiii10<-t3_filter(tiii10)
tiii11<-t3_filter(tiii11)
tiii12<-t3_filter(tiii12)
tiii13<-t3_filter(tiii13)
tiii14<-t3_filter(tiii14)
tiii15<-t3_filter(tiii15)
tiii16<-t3_filter(tiii16)

tiii <- list(tiii1,tiii2,tiii3,tiii4,tiii6,tiii7,tiii8,tiii9,tiii10,tiii11,tiii12,tiii13,tiii14,tiii15,tiii16) %>% reduce(inner_join, by = "genus")
tiii <- data.frame(tiii)
tiii <- tiii[c(2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
colnames(tiii)<-c('taxon_name','T1','T2','T3','T4','T6','T7','T8',
                  'T9','T10','T11','T12','T13','T14','T15','T16')

tiii$rowsum <- rowSums(tiii[,2:16])

tiii_top200 <- tiii[order(-tiii[1:156,16]),]
rownames(tiii_top200) <- tiii_top200$taxon_name
pheatmap(log2(tiii_top200[1:40,2:16]), #display_numbers = order_by_A[1:50,2:6],
         scale = "row", #You can try to change scale if you want
         cluster_rows=T,
         cluster_cols=F,
         cutree_rows = 5,
         number_color = "grey30", #cellwidth = 90, cellheight = 10,
         margins=c(3,25),fontsize = 10,show_colnames = TRUE,show_rownames = TRUE,main = "Top 200")

######### ######### ######### ######### ######### ######### ######### #########

tiv1<-read.table("TIV_mgl/genus/TIV1_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv2<-read.table("TIV_mgl/genus/TIV2_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv3<-read.table("TIV_mgl/genus/TIV3_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv4<-read.table("TIV_mgl/genus/TIV4_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv6<-read.table("TIV_mgl/genus/TIV6_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv7<-read.table("TIV_mgl/genus/TIV7_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv8<-read.table("TIV_mgl/genus/TIV8_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv9<-read.table("TIV_mgl/genus/TIV9_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv10<-read.table("TIV_mgl/genus/TIV10_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv11<-read.table("TIV_mgl/genus/TIV11_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv12<-read.table("TIV_mgl/genus/TIV12_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv13<-read.table("TIV_mgl/genus/TIV13_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv14<-read.table("TIV_mgl/genus/TIV14_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv15<-read.table("TIV_mgl/genus/TIV14_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")
tiv16<-read.table("TIV_mgl/genus/TIV15_R1.fastq.gz.kaiju.out.kaiju_genus.tsv",header = T,sep="\t")

tiv1<-reduce_dataframe(tiv1)
tiv2<-reduce_dataframe(tiv2)
tiv3<-reduce_dataframe(tiv3)
tiv4<-reduce_dataframe(tiv4)
tiv6<-reduce_dataframe(tiv6)
tiv7<-reduce_dataframe(tiv7)
tiv8<-reduce_dataframe(tiv8)
tiv9<-reduce_dataframe(tiv9)
tiv10<-reduce_dataframe(tiv10)
tiv11<-reduce_dataframe(tiv11)
tiv12<-reduce_dataframe(tiv12)
tiv13<-reduce_dataframe(tiv13)
tiv14<-reduce_dataframe(tiv14)
tiv15<-reduce_dataframe(tiv15)
tiv16<-reduce_dataframe(tiv16)

tiv <- list(tiv1,tiv2,tiv3,tiv4,tiv6,tiv7,tiv8,tiv9,tiv10,tiv11,tiv12,tiv13,tiv14,tiv15,tiv16) %>% reduce(inner_join, by = "taxon_name")
tiv <- data.frame(tiv)
tiv <- tiv[c(2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
colnames(tiv)<-c('taxon_name','T1','T2','T3','T4','T6','T7','T8',
                 'T9','T10','T11','T12','T13','T14','T15','T16')

tiv$rowsum <- rowSums(tiv[,2:16])

tiv_top200 <- tiv[order(-tiv[1:200,16]),]
rownames(tiv_top200) <- tiv_top200$taxon_name
pheatmap(log2(tiv_top200[1:40,2:16]), #display_numbers = order_by_A[1:50,2:6],
         scale = "row", #You can try to change scale if you want
         cluster_rows=T,
         cluster_cols=F,
         cutree_rows = 5,
         number_color = "grey30", #cellwidth = 90, cellheight = 10,
         margins=c(3,25),fontsize = 10,show_colnames = TRUE,show_rownames = TRUE,main = "Top 200")


head(ti)
head(tii)
head(tiii)
head(tiv)

all_ti_join <- list(ti,tii,tiii,tiv) %>% reduce(inner_join, by = "taxon_name")
all_ti_join <- data.frame(all_ti_join)
head(all_ti_join)

rownames(all_ti_join) <- all_ti_join$taxon_name


pheatmap(log2(all_ti_join[5:50,2:16]), #display_numbers = order_by_A[1:50,2:6],
         scale = "row", #You can try to change scale if you want
         cluster_rows=T,
         cluster_cols=F,
         cutree_rows = 5,
         number_color = "grey30", #cellwidth = 90, cellheight = 10,
         margins=c(3,25),fontsize = 10,show_colnames = TRUE,show_rownames = TRUE,main = "Top 200")