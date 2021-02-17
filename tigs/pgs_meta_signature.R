
library(dplyr)
library(tidyr)
library(ggsci)
library(pheatmap)
library(ggplot2)


#########################################################
#
#     WILD Bangalore Poppulations
#
#########################################################
TI1 <- read.csv("TI/TI1_kaiju.PGS.table", header = T, sep="\t")
TI2 <- read.csv("TI/TI2_kaiju.PGS.table", header = T, sep="\t")
TI3 <- read.csv("TI/TI3_kaiju.PGS.table", header = T, sep="\t")
TI4 <- read.csv("TI/TI4_kaiju.PGS.table", header = T, sep="\t")
TI6 <- read.csv("TI/TI6_kaiju.PGS.table", header = T, sep="\t")
TI7 <- read.csv("TI/TI7_kaiju.PGS.table", header = T, sep="\t")
TI8 <- read.csv("TI/TI8_kaiju.PGS.table", header = T, sep="\t")
TI9 <- read.csv("TI/TI9_kaiju.PGS.table", header = T, sep="\t")
TI10 <- read.csv("TI/TI10_kaiju.PGS.table", header = T, sep="\t")
TI11 <- read.csv("TI/TI11_kaiju.PGS.table", header = T, sep="\t")
TI12 <- read.csv("TI/TI12_kaiju.PGS.table", header = T, sep="\t")
TI13 <- read.csv("TI/TI13_kaiju.PGS.table", header = T, sep="\t")
TI14 <- read.csv("TI/TI14_kaiju.PGS.table", header = T, sep="\t")
TI15 <- read.csv("TI/TI15_kaiju.PGS.table", header = T, sep="\t")
TI16 <- read.csv("TI/TI16_kaiju.PGS.table", header = T, sep="\t")

TI1 <- separate(data = TI1[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI2 <- separate(data = TI2[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI3 <- separate(data = TI3[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI4 <- separate(data = TI4[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI5 <- separate(data = TI5[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI6 <- separate(data = TI6[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI7 <- separate(data = TI7[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI8 <- separate(data = TI8[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI9 <- separate(data = TI9[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI10 <- separate(data = TI10[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI11 <- separate(data = TI11[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI12 <- separate(data = TI12[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI13 <- separate(data = TI13[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI14 <- separate(data = TI14[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI15 <- separate(data = TI15[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TI16 <- separate(data = TI16[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")

TI1 <- TI1[TI1$reads>100,]
TI2 <- TI2[TI2$reads>100,]
TI3 <- TI3[TI3$reads>100,]
TI4 <- TI4[TI4$reads>100,]
TI6 <- TI6[TI6$reads>100,]
TI7 <- TI7[TI7$reads>100,]
TI8 <- TI8[TI8$reads>100,]
TI9 <- TI9[TI9$reads>100,]
TI10 <- TI10[TI10$reads>100,]
TI11 <- TI11[TI11$reads>100,]
TI12 <- TI12[TI12$reads>100,]
TI13 <- TI13[TI13$reads>100,]
TI14 <- TI14[TI14$reads>100,]
TI15 <- TI15[TI15$reads>100,]
TI16 <- TI16[TI16$reads>100,]

TI1 <- TI1[complete.cases(TI1[,3:4]),]
TI2 <- TI2[complete.cases(TI2[,3:4]),]
TI3 <- TI3[complete.cases(TI3[,3:4]),]
TI4 <- TI4[complete.cases(TI4[,3:4]),]
TI6 <- TI6[complete.cases(TI6[,3:4]),]
TI7 <- TI7[complete.cases(TI7[,3:4]),]
TI8 <- TI8[complete.cases(TI8[,3:4]),]
TI9 <- TI9[complete.cases(TI9[,3:4]),]
TI10 <- TI10[complete.cases(TI10[,3:4]),]
TI11 <- TI11[complete.cases(TI11[,3:4]),]
TI12 <- TI12[complete.cases(TI12[,3:4]),]
TI13 <- TI13[complete.cases(TI13[,3:4]),]
TI14 <- TI14[complete.cases(TI14[,3:4]),]
TI15 <- TI15[complete.cases(TI15[,3:4]),]
TI16 <- TI16[complete.cases(TI16[,3:4]),]

TI1_species <- inner_join(TI1,TI2,by="Species") %>% inner_join(TI3, by="Species")%>% inner_join(TI4,by="Species") %>% 
  inner_join(TI6,by="Species") %>% inner_join(TI7,by="Species") %>% inner_join(TI8,by="Species") %>% 
  inner_join(TI9,by="Species") %>% inner_join(TI10,by="Species")  %>% inner_join(TI11,by="Species")  %>% inner_join(TI12,by="Species")  %>%
  inner_join(TI13,by="Species")  %>% inner_join(TI14,by="Species") %>% inner_join(TI15,by="Species") %>% inner_join(TI16,by="Species")

TI1_species <- TI1_species[,c(2,3,4,1,5,8,11,14,17,20,23,26,29,32,35,38,41,44)]
colnames(TI1_species)<-c("phylum","genus","species","TI1","TI2","TI3","TI4","TI6","TI7","TI8","TI9","TI10","TI11","TI12","TI13","TI14","TI15","TI16")

#########################################################
#
#     Lab Chennai Poppulations
#
#########################################################

TII1 <- read.csv("TII/TII1_kaiju.PGS.table", header = T, sep="\t")
TII2 <- read.csv("TII/TII2_kaiju.PGS.table", header = T, sep="\t")
TII3 <- read.csv("TII/TII3_kaiju.PGS.table", header = T, sep="\t")
TII4 <- read.csv("TII/TII4_kaiju.PGS.table", header = T, sep="\t")
TII6 <- read.csv("TII/TII6_kaiju.PGS.table", header = T, sep="\t")
TII7 <- read.csv("TII/TII7_kaiju.PGS.table", header = T, sep="\t")
TII8 <- read.csv("TII/TII8_kaiju.PGS.table", header = T, sep="\t")
TII9 <- read.csv("TII/TII9_kaiju.PGS.table", header = T, sep="\t")
TII10 <- read.csv("TII/TII10_kaiju.PGS.table", header = T, sep="\t")
TII11 <- read.csv("TII/TII11_kaiju.PGS.table", header = T, sep="\t")
TII12 <- read.csv("TII/TII12_kaiju.PGS.table", header = T, sep="\t")
TII13 <- read.csv("TII/TII13_kaiju.PGS.table", header = T, sep="\t")
TII14 <- read.csv("TII/TII14_kaiju.PGS.table", header = T, sep="\t")
TII15 <- read.csv("TII/TII15_kaiju.PGS.table", header = T, sep="\t")
TII16 <- read.csv("TII/TII16_kaiju.PGS.table", header = T, sep="\t")

TII1 <- separate(data = TII1[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII2 <- separate(data = TII2[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII3 <- separate(data = TII3[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII4 <- separate(data = TII4[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII5 <- separate(data = TII5[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII6 <- separate(data = TII6[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII7 <- separate(data = TII7[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII8 <- separate(data = TII8[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII9 <- separate(data = TII9[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII10 <- separate(data = TII10[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII11 <- separate(data = TII11[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII12 <- separate(data = TII12[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII13 <- separate(data = TII13[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII14 <- separate(data = TII14[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII15 <- separate(data = TII15[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TII16 <- separate(data = TII16[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")

TII1 <- TII1[TII1$reads>100,]
TII2 <- TII2[TII2$reads>100,]
TII3 <- TII3[TII3$reads>100,]
TII4 <- TII4[TII4$reads>100,]
TII6 <- TII6[TII6$reads>100,]
TII7 <- TII7[TII7$reads>100,]
TII8 <- TII8[TII8$reads>100,]
TII9 <- TII9[TII9$reads>100,]
TII10 <- TII10[TII10$reads>100,]
TII11 <- TII11[TII11$reads>100,]
TII12 <- TII12[TII12$reads>100,]
TII13 <- TII13[TII13$reads>100,]
TII14 <- TII14[TII14$reads>100,]
TII15 <- TII15[TII15$reads>100,]
TII16 <- TII16[TII16$reads>100,]

TII1 <- TII1[complete.cases(TII1[,3:4]),]
TII2 <- TII2[complete.cases(TII2[,3:4]),]
TII3 <- TII3[complete.cases(TII3[,3:4]),]
TII4 <- TII4[complete.cases(TII4[,3:4]),]
TII6 <- TII6[complete.cases(TII6[,3:4]),]
TII7 <- TII7[complete.cases(TII7[,3:4]),]
TII8 <- TII8[complete.cases(TII8[,3:4]),]
TII9 <- TII9[complete.cases(TII9[,3:4]),]
TII10 <- TII10[complete.cases(TII10[,3:4]),]
TII11 <- TII11[complete.cases(TII11[,3:4]),]
TII12 <- TII12[complete.cases(TII12[,3:4]),]
TII13 <- TII13[complete.cases(TII13[,3:4]),]
TII14 <- TII14[complete.cases(TII14[,3:4]),]
TII15 <- TII15[complete.cases(TII15[,3:4]),]
TII16 <- TII16[complete.cases(TII16[,3:4]),]

TII1_species <- inner_join(TII1,TII2,by="Species") %>% inner_join(TII3, by="Species")%>% inner_join(TII4,by="Species") %>% 
  inner_join(TII6,by="Species") %>% inner_join(TII7,by="Species") %>% inner_join(TII8,by="Species") %>% 
  inner_join(TII9,by="Species") %>% inner_join(TII10,by="Species")  %>% inner_join(TII11,by="Species")  %>% inner_join(TII12,by="Species")  %>%
  inner_join(TII13,by="Species")  %>% inner_join(TII14,by="Species") %>% inner_join(TII15,by="Species") %>% inner_join(TII16,by="Species")

TII1_species <- TII1_species[,c(2,3,4,1,5,8,11,14,17,20,23,26,29,32,35,38,41,44)]
colnames(TII1_species)<-c("phylum","genus","species","TII1","TII2","TII3","TII4","TII6","TII7","TII8","TII9","TII10","TII11","TII12","TII13","TII14","TII15","TII16")

#########################################################
#
#     Lab Delhi Poppulations
#
#########################################################
TIII1 <- read.csv("TIII/TIII1_kaiju.PGS.table", header = T, sep="\t")
TIII2 <- read.csv("TIII/TIII2_kaiju.PGS.table", header = T, sep="\t")
TIII3 <- read.csv("TIII/TIII3_kaiju.PGS.table", header = T, sep="\t")
TIII4 <- read.csv("TIII/TIII4_kaiju.PGS.table", header = T, sep="\t")
TIII5 <- read.csv("TIII/TIII5_kaiju.PGS.table", header = T, sep="\t")
TIII6 <- read.csv("TIII/TIII6_kaiju.PGS.table", header = T, sep="\t")
TIII7 <- read.csv("TIII/TIII7_kaiju.PGS.table", header = T, sep="\t")
TIII8 <- read.csv("TIII/TIII8_kaiju.PGS.table", header = T, sep="\t")
TIII9 <- read.csv("TIII/TIII9_kaiju.PGS.table", header = T, sep="\t")
TIII10 <- read.csv("TIII/TIII10_kaiju.PGS.table", header = T, sep="\t")
TIII11 <- read.csv("TIII/TIII11_kaiju.PGS.table", header = T, sep="\t")
TIII12 <- read.csv("TIII/TIII12_kaiju.PGS.table", header = T, sep="\t")
TIII13 <- read.csv("TIII/TIII13_kaiju.PGS.table", header = T, sep="\t")
TIII14 <- read.csv("TIII/TIII14_kaiju.PGS.table", header = T, sep="\t")
TIII15 <- read.csv("TIII/TIII15_kaiju.PGS.table", header = T, sep="\t")
TIII16 <- read.csv("TIII/TIII16_kaiju.PGS.table", header = T, sep="\t")

TIII1 <- separate(data = TIII1[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII2 <- separate(data = TIII2[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII3 <- separate(data = TIII3[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII4 <- separate(data = TIII4[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII5 <- separate(data = TIII5[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII6 <- separate(data = TIII6[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII7 <- separate(data = TIII7[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII8 <- separate(data = TIII8[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII9 <- separate(data = TIII9[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII10 <- separate(data = TIII10[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII11 <- separate(data = TIII11[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII12 <- separate(data = TIII12[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII13 <- separate(data = TIII13[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII14 <- separate(data = TIII14[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII15 <- separate(data = TIII15[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIII16 <- separate(data = TIII16[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")

TIII1 <- TIII1[TIII1$reads>100,]
TIII2 <- TIII2[TIII2$reads>100,]
TIII3 <- TIII3[TIII3$reads>100,]
TIII4 <- TIII4[TIII4$reads>100,]
TIII5 <- TIII5[TIII5$reads>100,]
TIII6 <- TIII6[TIII6$reads>100,]
TIII7 <- TIII7[TIII7$reads>100,]
TIII8 <- TIII8[TIII8$reads>100,]
TIII9 <- TIII9[TIII9$reads>100,]
TIII10 <- TIII10[TIII10$reads>100,]
TIII11 <- TIII11[TIII11$reads>100,]
TIII12 <- TIII12[TIII12$reads>100,]
TIII13 <- TIII13[TIII13$reads>100,]
TIII14 <- TIII14[TIII14$reads>100,]
TIII15 <- TIII15[TIII15$reads>100,]
TIII16 <- TIII16[TIII16$reads>100,]

TIII1 <- TIII1[complete.cases(TIII1[,3:4]),]
TIII2 <- TIII2[complete.cases(TIII2[,3:4]),]
TIII3 <- TIII3[complete.cases(TIII3[,3:4]),]
TIII4 <- TIII4[complete.cases(TIII4[,3:4]),]
TIII5 <- TIII5[complete.cases(TIII5[,3:4]),]
TIII6 <- TIII6[complete.cases(TIII6[,3:4]),]
TIII7 <- TIII7[complete.cases(TIII7[,3:4]),]
TIII8 <- TIII8[complete.cases(TIII8[,3:4]),]
TIII9 <- TIII9[complete.cases(TIII9[,3:4]),]
TIII10 <- TIII10[complete.cases(TIII10[,3:4]),]
TIII11 <- TIII11[complete.cases(TIII11[,3:4]),]
TIII12 <- TIII12[complete.cases(TIII12[,3:4]),]
TIII13 <- TIII13[complete.cases(TIII13[,3:4]),]
TIII14 <- TIII14[complete.cases(TIII14[,3:4]),]
TIII15 <- TIII15[complete.cases(TIII15[,3:4]),]
TIII16 <- TIII16[complete.cases(TIII16[,3:4]),]

TIII1_species <- inner_join(TIII1,TIII2,by="Species") %>% inner_join(TIII3, by="Species")%>% inner_join(TIII4,by="Species") %>%  inner_join(TIII5,by="Species") %>% 
  inner_join(TIII6,by="Species") %>% inner_join(TIII7,by="Species") %>% inner_join(TIII8,by="Species") %>% 
  inner_join(TIII9,by="Species") %>% inner_join(TIII10,by="Species")  %>% inner_join(TIII11,by="Species")  %>% inner_join(TIII12,by="Species")  %>%
  inner_join(TIII13,by="Species")  %>% inner_join(TIII14,by="Species") %>% inner_join(TIII15,by="Species") %>% inner_join(TIII16,by="Species")

TIII1_species <- TIII1_species[,c(2,3,4,1,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47)]
colnames(TIII1_species)<-c("phylum","genus","species","TIII1","TIII2","TIII3","TIII4","TIII5","TIII6","TIII7","TIII8","TIII9","TIII10","TIII11","TIII12","TIII13","TIII14","TIII15","TIII16")


#########################################################
#
#     Lab Mangalore Poppulations
#
#########################################################

TIV1 <- read.csv("TIV/TIV1_kaiju.PGS.table", header = T, sep="\t")
TIV2 <- read.csv("TIV/TIV2_kaiju.PGS.table", header = T, sep="\t")
TIV3 <- read.csv("TIV/TIV3_kaiju.PGS.table", header = T, sep="\t")
TIV4 <- read.csv("TIV/TIV4_kaiju.PGS.table", header = T, sep="\t")
TIV5 <- read.csv("TIV/TIV5_kaiju.PGS.table", header = T, sep="\t")
TIV6 <- read.csv("TIV/TIV6_kaiju.PGS.table", header = T, sep="\t")
TIV7 <- read.csv("TIV/TIV7_kaiju.PGS.table", header = T, sep="\t")
TIV8 <- read.csv("TIV/TIV8_kaiju.PGS.table", header = T, sep="\t")
TIV9 <- read.csv("TIV/TIV9_kaiju.PGS.table", header = T, sep="\t")
TIV10 <- read.csv("TIV/TIV10_kaiju.PGS.table", header = T, sep="\t")
TIV11 <- read.csv("TIV/TIV11_kaiju.PGS.table", header = T, sep="\t")
TIV12 <- read.csv("TIV/TIV12_kaiju.PGS.table", header = T, sep="\t")
TIV13 <- read.csv("TIV/TIV13_kaiju.PGS.table", header = T, sep="\t")
TIV14 <- read.csv("TIV/TIV14_kaiju.PGS.table", header = T, sep="\t")
TIV15 <- read.csv("TIV/TIV15_kaiju.PGS.table", header = T, sep="\t")
TIV16 <- read.csv("TIV/TIV16_kaiju.PGS.table", header = T, sep="\t")

TIV1 <- separate(data = TIV1[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV2 <- separate(data = TIV2[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV3 <- separate(data = TIV3[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV4 <- separate(data = TIV4[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV5 <- separate(data = TIV5[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV6 <- separate(data = TIV6[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV7 <- separate(data = TIV7[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV8 <- separate(data = TIV8[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV9 <- separate(data = TIV9[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV10 <- separate(data = TIV10[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV11 <- separate(data = TIV11[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV12 <- separate(data = TIV12[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV13 <- separate(data = TIV13[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV14 <- separate(data = TIV14[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV15 <- separate(data = TIV15[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")
TIV16 <- separate(data = TIV16[,c(3,5)], col = taxon_name, into = c("Phylum","Genus","Species"), sep = ";")

TIV1 <- TIV1[TIV1$reads>100,]
TIV2 <- TIV2[TIV2$reads>100,]
TIV3 <- TIV3[TIV3$reads>100,]
TIV4 <- TIV4[TIV4$reads>100,]
TIV5 <- TIV5[TIV5$reads>100,]
TIV6 <- TIV6[TIV6$reads>100,]
TIV7 <- TIV7[TIV7$reads>100,]
TIV8 <- TIV8[TIV8$reads>100,]
TIV9 <- TIV9[TIV9$reads>100,]
TIV10 <- TIV10[TIV10$reads>100,]
TIV11 <- TIV11[TIV11$reads>100,]
TIV12 <- TIV12[TIV12$reads>100,]
TIV13 <- TIV13[TIV13$reads>100,]
TIV14 <- TIV14[TIV14$reads>100,]
TIV15 <- TIV15[TIV15$reads>100,]
TIV16 <- TIV16[TIV16$reads>100,]

TIV1 <- TIV1[complete.cases(TIV1[,3:4]),]
TIV2 <- TIV2[complete.cases(TIV2[,3:4]),]
TIV3 <- TIV3[complete.cases(TIV3[,3:4]),]
TIV4 <- TIV4[complete.cases(TIV4[,3:4]),]
TIV5 <- TIV5[complete.cases(TIV5[,3:4]),]
TIV6 <- TIV6[complete.cases(TIV6[,3:4]),]
TIV7 <- TIV7[complete.cases(TIV7[,3:4]),]
TIV8 <- TIV8[complete.cases(TIV8[,3:4]),]
TIV9 <- TIV9[complete.cases(TIV9[,3:4]),]
TIV10 <- TIV10[complete.cases(TIV10[,3:4]),]
TIV11 <- TIV11[complete.cases(TIV11[,3:4]),]
TIV12 <- TIV12[complete.cases(TIV12[,3:4]),]
TIV13 <- TIV13[complete.cases(TIV13[,3:4]),]
TIV14 <- TIV14[complete.cases(TIV14[,3:4]),]
TIV15 <- TIV15[complete.cases(TIV15[,3:4]),]
TIV16 <- TIV16[complete.cases(TIV16[,3:4]),]

TIV1_species <- inner_join(TIV1,TIV2,by="Species") %>% inner_join(TIV3, by="Species")%>% inner_join(TIV4,by="Species") %>% inner_join(TIV5,by="Species") %>% 
  inner_join(TIV6,by="Species") %>% inner_join(TIV7,by="Species") %>% inner_join(TIV8,by="Species") %>% 
  inner_join(TIV9,by="Species") %>% inner_join(TIV10,by="Species")  %>% inner_join(TIV11,by="Species")  %>% inner_join(TIV12,by="Species")  %>%
  inner_join(TIV13,by="Species")  %>% inner_join(TIV14,by="Species") %>% inner_join(TIV15,by="Species") %>% inner_join(TIV16,by="Species")

TIV1_species <- TIV1_species[,c(2,3,4,1,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47)]
colnames(TIV1_species)<-c("phylum","genus","species","TIV1","TIV2","TIV3","TIV4","TIV5","TIV6","TIV7","TIV8","TIV9","TIV10","TIV11","TIV12","TIV13","TIV14","TIV15","TIV16")


#########################################################
#
#     Writing in table all the population Poppulations
#
#########################################################

write.table(TI1_species, "TI1_Species", quote = F, row.names = F, col.names = T, sep="\t")
write.table(TII1_species, "TII1_Species", quote = F, row.names = F, col.names = T, sep="\t")
write.table(TIII1_species, "TIII1_Species", quote = F, row.names = F, col.names = T, sep="\t")
write.table(TIV1_species, "TIV1_Species", quote = F, row.names = F, col.names = T, sep="\t")


all_population_join <- full_join(TI1_species,TII1_species,by="species") %>% full_join(TIII1_species,by="species") %>%  full_join(TIV1_species,by="species")
all_population_join <- all_population_join[,-c(19,20,36,37,54,55)]
all_population_join[is.na(all_population_join),] <- 0
all_population_join$T1Rowsum <- rowSums(all_population_join[,9:18])
all_population_join$T2Rowsum <- rowSums(all_population_join[,19:33])
all_population_join$T3Rowsum <- rowSums(all_population_join[,34:49])
all_population_join$T4Rowsum <- rowSums(all_population_join[,50:65])

for_fiter <- all_population_join[,c(1,2,3,66:69)]

zero_row <- for_fiter[rowSums(for_fiter[4:7] == 0) > 0, ]
zero_row <- zero_row[,3:7]
rownames(zero_row) <- zero_row$species

zero_row[zero_row == 0] <- 1
nrow(zero_row)
pheatmap(log2(zero_row[1:292,2:5]),
         show_rownames = F,
         cluster_rows = T,
         cluster_cols = F)
