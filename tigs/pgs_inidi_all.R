# WILD Bangalore Poppulations
library(dplyr)
library(tidyr)
library(DESeq2)
library(genefilter)
library(pheatmap)
library(ggplot2)


b1 <- read.csv("B/B1_kaiju.PGS.table", header = T, sep="\t")
b2 <- read.csv("B/B2_kaiju.PGS.table", header = T, sep="\t")
b3 <- read.csv("B/B3_kaiju.PGS.table", header = T, sep="\t")
b4 <- read.csv("B/B4_kaiju.PGS.table", header = T, sep="\t")
b5 <- read.csv("B/B5_kaiju.PGS.table", header = T, sep="\t")
b6 <- read.csv("B/B6_kaiju.PGS.table", header = T, sep="\t")
b7 <- read.csv("B/B7_kaiju.PGS.table", header = T, sep="\t")
b8 <- read.csv("B/B8_kaiju.PGS.table", header = T, sep="\t")
b9 <- read.csv("B/B9_kaiju.PGS.table", header = T, sep="\t")
b10 <- read.csv("B/B10_kaiju.PGS.table", header = T, sep="\t")
head(b1[,c(3,5)])

b1 <- separate(data = b1[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
b2 <- separate(data = b2[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
b3 <- separate(data = b3[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
b4 <- separate(data = b4[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
b5 <- separate(data = b5[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
b6 <- separate(data = b6[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
b7 <- separate(data = b7[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
b8 <- separate(data = b8[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
b9 <- separate(data = b9[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
b10 <- separate(data = b10[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")

join_species <- inner_join(b1,b2,by="Species") %>% inner_join(b3,by="Species")%>% inner_join(b4,by="Species") %>% 
  inner_join(b5,by="Species") %>% inner_join(b6,by="Species") %>% inner_join(b7,by="Species") %>% inner_join(b8,by="Species") %>% 
  inner_join(b9,by="Species") %>% inner_join(b10,by="Species")

join_species <- join_species[,c(2,3,4,1,5,8,11,14,17,20,23)]

# Wild Mangalore Poppulations

m1 <- read.csv("M/M1_kaiju.PGS.table", header = T, sep="\t")
m2 <- read.csv("M/M2_kaiju.PGS.table", header = T, sep="\t")
m3 <- read.csv("M/M3_kaiju.PGS.table", header = T, sep="\t")
m4 <- read.csv("M/M4_kaiju.PGS.table", header = T, sep="\t")
m5 <- read.csv("M/M5_kaiju.PGS.table", header = T, sep="\t")
m6 <- read.csv("M/M6_kaiju.PGS.table", header = T, sep="\t")
m7 <- read.csv("M/M7_kaiju.PGS.table", header = T, sep="\t")
m8 <- read.csv("M/M8_kaiju.PGS.table", header = T, sep="\t")
m9 <- read.csv("M/M9_kaiju.PGS.table", header = T, sep="\t")
m10 <- read.csv("M/M10_kaiju.PGS.table", header = T, sep="\t")
head(m9[,])

m1 <- separate(data = m1[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
m2 <- separate(data = m2[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
m3 <- separate(data = m3[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
m4 <- separate(data = m4[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
m5 <- separate(data = m5[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
m6 <- separate(data = m6[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
m7 <- separate(data = m7[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
m8 <- separate(data = m8[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
m9 <- separate(data = m9[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
m10 <- separate(data = m10[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")

m1 <- m1[m1$reads>100,]
m2 <- m2[m2$reads>100,]
m3 <- m3[m3$reads>100,]
m4 <- m4[m4$reads>100,]
m5 <- m5[m5$reads>100,]
m6 <- m6[m6$reads>100,]
m7 <- m7[m7$reads>100,]
m8 <- m8[m8$reads>100,]
m9 <- m9[m9$reads>100,]
m10 <- m10[m10$reads>100,]

m1 <- m1[complete.cases(m1[,3:4]),]
m2 <- m2[complete.cases(m2[,3:4]),]
m3 <- m3[complete.cases(m3[,3:4]),]
m4 <- m4[complete.cases(m4[,3:4]),]
m5 <- m5[complete.cases(m5[,3:4]),]
m6 <- m6[complete.cases(m6[,3:4]),]
m7 <- m7[complete.cases(m7[,3:4]),]
m8 <- m8[complete.cases(m8[,3:4]),]
m9 <- m9[complete.cases(m9[,3:4]),]
m10 <- m10[complete.cases(m10[,3:4]),]

join_species <- inner_join(m1,m2,by="Species") %>% inner_join(m3,by="Species")%>% inner_join(m4,by="Species") %>% 
  inner_join(m5,by="Species") %>% inner_join(m6,by="Species") %>% inner_join(m7,by="Species") %>% inner_join(m8,by="Species") %>% 
  inner_join(m9,by="Species") %>% inner_join(m10,by="Species")

wild_m_join_species <- join_species[,c(2,3,4,1,5,8,11,14,17,20,23,26,29)]
colnames(wild_m_join_species)<-c("phylum","genus","species","mi","mii","miii","miv","mv","mvi","mvii","mviii","mixm","mx")

# Lab Mangloare     

ml1 <- read.csv("TIV/TIV1_kaiju.PGS.table", header = T, sep="\t")
ml2 <- read.csv("TIV/TIV2_kaiju.PGS.table", header = T, sep="\t")
ml3 <- read.csv("TIV/TIV3_kaiju.PGS.table", header = T, sep="\t")
ml4 <- read.csv("TIV/TIV4_kaiju.PGS.table", header = T, sep="\t")
ml5 <- read.csv("TIV/TIV5_kaiju.PGS.table", header = T, sep="\t")
ml6 <- read.csv("TIV/TIV6_kaiju.PGS.table", header = T, sep="\t")
ml7 <- read.csv("TIV/TIV7_kaiju.PGS.table", header = T, sep="\t")
ml8 <- read.csv("TIV/TIV8_kaiju.PGS.table", header = T, sep="\t")
ml9 <- read.csv("TIV/TIV9_kaiju.PGS.table", header = T, sep="\t")
ml10 <- read.csv("TIV/TIV10_kaiju.PGS.table", header = T, sep="\t")
ml11 <- read.csv("TIV/TIV11_kaiju.PGS.table", header = T, sep="\t")
ml12 <- read.csv("TIV/TIV12_kaiju.PGS.table", header = T, sep="\t")
ml13 <- read.csv("TIV/TIV13_kaiju.PGS.table", header = T, sep="\t")
ml14 <- read.csv("TIV/TIV14_kaiju.PGS.table", header = T, sep="\t")
ml15 <- read.csv("TIV/TIV15_kaiju.PGS.table", header = T, sep="\t")
ml16 <- read.csv("TIV/TIV16_kaiju.PGS.table", header = T, sep="\t")
head(ml9[,])


ml1 <- separate(data = ml1[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml2 <- separate(data = ml2[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml3 <- separate(data = ml3[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml4 <- separate(data = ml4[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml5 <- separate(data = ml5[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml6 <- separate(data = ml6[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml7 <- separate(data = ml7[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml8 <- separate(data = ml8[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml9 <- separate(data = ml9[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml10 <- separate(data = ml10[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml11 <- separate(data = ml11[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml12 <- separate(data = ml12[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml13 <- separate(data = ml13[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml14 <- separate(data = ml14[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml15 <- separate(data = ml15[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
ml16 <- separate(data = ml16[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")

ml1 <- ml1[ml1$reads>100,]
ml2 <- ml2[ml2$reads>100,]
ml3 <- ml3[ml3$reads>100,]
ml4 <- ml4[ml4$reads>100,]
ml5 <- ml5[ml5$reads>100,]
ml6 <- ml6[ml6$reads>100,]
ml7 <- ml7[ml7$reads>100,]
ml8 <- ml8[ml8$reads>100,]
ml9 <- ml9[ml9$reads>100,]
ml10 <- ml10[ml10$reads>100,]
ml11 <- ml11[ml11$reads>100,]
ml12 <- ml12[ml12$reads>100,]
ml13 <- ml13[ml13$reads>100,]
ml14 <- ml14[ml14$reads>100,]
ml15 <- ml15[ml15$reads>100,]
ml16 <- ml16[ml16$reads>100,]

ml1 <- ml1[complete.cases(ml1[,3:4]),]
ml2 <- ml2[complete.cases(ml2[,3:4]),]
ml3 <- ml3[complete.cases(ml3[,3:4]),]
ml4 <- ml4[complete.cases(ml4[,3:4]),]
ml5 <- ml5[complete.cases(ml5[,3:4]),]
ml6 <- ml6[complete.cases(ml6[,3:4]),]
ml7 <- ml7[complete.cases(ml7[,3:4]),]
ml8 <- ml8[complete.cases(ml8[,3:4]),]
ml9 <- ml9[complete.cases(ml9[,3:4]),]
ml10 <- ml10[complete.cases(ml10[,3:4]),]
ml11 <- ml11[complete.cases(ml11[,3:4]),]
ml12 <- ml12[complete.cases(ml12[,3:4]),]
ml13 <- ml13[complete.cases(ml13[,3:4]),]
ml14 <- ml14[complete.cases(ml14[,3:4]),]
ml15 <- ml15[complete.cases(ml15[,3:4]),]
ml16 <- ml16[complete.cases(ml16[,3:4]),]

join_species <- inner_join(ml1,ml2,by="Species") %>% inner_join(ml3,by="Species") %>% inner_join(ml4,by="Species") %>%  inner_join(ml5,by="Species") %>% 
  inner_join(ml6,by="Species") %>% inner_join(ml7,by="Species") %>% inner_join(ml8,by="Species") %>% 
  inner_join(ml9,by="Species") %>% inner_join(ml10,by="Species") %>% inner_join(ml11,by="Species") %>%  inner_join(ml12,by="Species") %>%
inner_join(ml13,by="Species") %>%  inner_join(ml14,by="Species") %>%  inner_join(ml15,by="Species") %>%  inner_join(ml16,by="Species") 

lab_m_join_species <- join_species[,c(2,3,4,1,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47)]
colnames(lab_m_join_species)<-c("phylum","genus","species","mli","mlii","mliii","mliv","mlv","mlvi","mlvii","mlviii","mlixm","mlx",
                                "mlxi","mlxii","mlxiii","mlxiv","mlxv","mlxvi")


wild_lab_full_join <- full_join(wild_m_join_species, lab_m_join_species, by="species")
wild_lab_anti_join <- anti_join(wild_m_join_species,lab_m_join_species,by="species")
wild_lab_inner_join <- inner_join(wild_m_join_species,lab_m_join_species,by="species")

for_deseq <- wild_lab_inner_join[,-c(1,2,14,15)]
print(colnames(for_deseq))
mycounts <- wild_lab_inner_join


metadata <-  read.csv("coldata.csv",header = T, sep = "\t", stringsAsFactors = F)


dds <- DESeqDataSetFromMatrix(countData=for_deseq, 
                              colData=metadata, 
                              design=~population, 
                              tidy=TRUE)
dds
dds <- DESeq(dds)
res <- results(dds)
res
plotMA(res, ylim=c(-2,2))
res %>% ggplot(aes(baseMean, log2FoldChange, col=sig)) + geom_point() + scale_x_log10() + ggtitle("MA plot")
res <- as.data.frame(res)
res %>% filter(padj<0.05)
View(res)

rld <- rlogTransformation(dds)
plotPCA(rld, intgroup="population")
NMF::aheatmap(assay(rld)[arrange(res, padj, pvalue)$row[1:50], ], labRow=arrange(res, padj, pvalue)$symbol[1:50], scale="row", distfun="pearson", annCol=dplyr::select(metadata, population, celltype), col=c("green","black","black","red"))

res %>% ggplot(aes(log2FoldChange, -1*log10(pvalue), col=padj)) + geom_point() + ggtitle("Volcano plot")


topVarGenes <- head(order(-rowVars(assay(rld))),20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("id","population")])
pheatmap(mat, annotation_col=df)

write.table(res,"dds_result",quote = F)
inner_join()


write.table(res_join,"dds_join_kaiju",quote = F,sep="\t")
pheatmap(log2(res_join[10:100,2:25]))


colnames(res_join[,2:27])
res$species <- rownames(res)
res_join <- inner_join(for_deseq,res,by="species",sep="\t")
pheatmap(res_join[100:170,2:27])
  