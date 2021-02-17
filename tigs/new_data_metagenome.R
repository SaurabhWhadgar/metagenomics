setwd("D:/IBAB_Data/tigs/new_data/pooled")
source("../genus/kaju_function.R")

library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(ggsci)


bang <- read.csv("species/Bangalore_All_kaiju.All.table",header = T,sep="\t")
mang <- read.csv("species/Mangalore_All_kaiju.All.table",header = T,sep="\t")
laksh <- read.csv("species/Lakshadweep_All_kaiju.All.table",header = T,sep="\t")
TI <- read.csv("species/TI_All_kaiju.All.table",header = T,sep="\t")
TII <- read.csv("species/TII_All_kaiju.All.table",header = T,sep="\t")
TIII <- read.csv("species/TIII_All_kaiju.All.table",header = T,sep="\t")
TIV <- read.csv("species/TIV_All_kaiju.All.table",header = T,sep="\t")

### Separate kaiju data ###
bang_sep_taxon <- kaiju_separate(bang,10)
mang_sep_taxon <- kaiju_separate(mang,10)
sapply(mang_sep_taxon,class)

laksh_sep_taxon <- kaiju_separate(laksh,10)
ti_sep_taxon <- kaiju_separate(TI,10)
tii_sep_taxon <- kaiju_separate(TII,10)
tiii_sep_taxon <- kaiju_separate(TIII,10)
tiv_sep_taxon <- kaiju_separate(TIV,10)

### 
bang_phylum <- head(kaiju_taxon_stat(bang_sep_taxon, phylum),10)
mang_phylum <- head(kaiju_taxon_stat(mang_sep_taxon, phylum),10)
laksh_phylum <- head(kaiju_taxon_stat(laksh_sep_taxon, phylum),10)
ti_phylum <- head(kaiju_taxon_stat(ti_sep_taxon, phylum),10)
tii_phylum <- head(kaiju_taxon_stat(tii_sep_taxon, phylum),10)
tiii_phlyum <- head(kaiju_taxon_stat(tiii_sep_taxon, phylum),10)
tiv_phylum <- head(kaiju_taxon_stat(tiv_sep_taxon, phylum),10)

phylum_data <- data.frame(kaiju_join("full",bang_phylum, mang_phylum, laksh_phylum, ti_phylum, tii_phylum, tiii_phlyum, tiv_phylum))

colnames(phylum_data) <- c("Phylum","Bang","Mang","Laksh","TI","TII","TIII","TIV")

data_percentage <- apply(phylum_data[,-1], 2, function(x){x*100/sum(x,na.rm=T)})

rownames(data_percentage) <- phylum_data[,1]

write.table(phylum_data,"Pooled_Phylum_joined.tsv",sep="\t",quote = F,row.names = F)

data = melt(data_percentage)

ggplot(data=data, aes(x = Var2, y = value, fill = Var1)) +  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Paired")

ggplot(data, aes(x = Var2)) +   geom_bar(aes(y = value/sum(value), fill = cut)) + 
  scale_fill_brewer(palette = "Set3") +   ylab("Percent") + ggtitle("Show precentages in bar chart")

data

bang_taxon <- kaiju_groupby(bang,"Bacteria",genus)
mang_taxon <- kaiju_groupby(mang,"Bacteria",genus)
laksh_taxon <- kaiju_groupby(laksh,"Bacteria",genus)
TI_taxon <- kaiju_groupby(TI,"Bacteria",genus)
TII_taxon <- kaiju_groupby(TII,"Bacteria",genus)
TIII_taxon <- kaiju_groupby(TIII,"Bacteria",genus)
TIV_taxon <- kaiju_groupby(TIV,"Bacteria",genus)

join_taxon <- data.frame(inner_join(bang_taxon,mang_taxon,by="phylum") %>% inner_join(laksh_taxon,by="phylum")%>% 
  inner_join(TI_taxon,by="phylum") %>% inner_join(TII_taxon,by="phylum") %>% inner_join(TIII_taxon,by="phylum") %>% 
  inner_join(TIV_taxon,by="phylum"))

join_taxon <- data.frame(inner_join(bang_taxon,mang_taxon,by="genus") %>% inner_join(laksh_taxon,by="genus")%>% 
                           inner_join(TI_taxon,by="genus") %>% inner_join(TII_taxon,by="genus") %>% inner_join(TIII_taxon,by="genus") %>% 
                           inner_join(TIV_taxon,by="genus"))

colnames(join_taxon)<-c("taxon","Bang","Mang","Laskh","TI","TII","TIII","TIV")

#write.table(join_taxon,"Pooled_All_joined_phylum.csv",sep=",",quote = F,row.names = F)

data_abundance <- data.frame(join_taxon[1:20,-1])
data2 <- data.frame(apply(data_abundance, 2, function(x){x*100/sum(x,na.rm=T)}))
data2
data2$taxon <- join_taxon[1:20,1]
data2
colnames(data2)<-c("Bang","Mang","Laskh","TI","TII","TIII","TIV","taxon")
data2 = melt(data2,id.vars = "taxon")

data2

ggplot(data=data2, aes(x = variable, y = value, fill = taxon)) +  geom_bar(stat="identity",position=position_dodge())

ggplot(data=data2, aes(x = variable, y = value, fill = taxon)) +  geom_bar(stat = "identity") 


bang_diversity <- kaiju_diversity(bang_phlyum)
bind_phylum <- kaiju_join("left",bang_phlyum,mang_phlyum,laksh_phlyum,TI_phlyum,TII_phlyum,TIII_phlyum,TIV_phlyum)

#---------------------
#
# Diversity Plot
#
#----------------------

#===================== For Paper at genus IN phylum Level ================

bang_phylum <- kaiju_taxon_stat(bang_sep_taxon, phylum)
mang_phylum <- kaiju_taxon_stat(mang_sep_taxon, phylum)
laksh_phylum <- kaiju_taxon_stat(laksh_sep_taxon, phylum)
ti_phylum <- kaiju_taxon_stat(ti_sep_taxon, phylum)
tii_phylum <- kaiju_taxon_stat(tii_sep_taxon, phylum)
tiii_phlyum <- kaiju_taxon_stat(tiii_sep_taxon, phylum)
tiv_phylum <- kaiju_taxon_stat(tiv_sep_taxon, phylum)
phylum_data <- data.frame(kaiju_join("full",bang_phylum, mang_phylum, laksh_phylum, ti_phylum, tii_phylum, tiii_phlyum, tiv_phylum))

for_phylumGenusDiversityPlot <- phylum_data[1:10,]
for_phylumGenusDiversityPlot[11,] <- c("Other",as.numeric(colSums(phylum_data[11:nrow(phylum_data),2:ncol(phylum_data)])))

colnames(for_phylumGenusDiversityPlot) <- c("Phylum","Bang","Mang","Laksh","TI","TII","TIII","TIV")

data_percentage <- apply(as.data.frame(sapply(for_phylumGenusDiversityPlot[,-1], as.numeric)), 2, function(x){x*100/sum(x,na.rm=T)})

rownames(data_percentage) <- for_phylumGenusDiversityPlot[,1]

data = melt(data_percentage)

ggplot(data=data, aes(x = Var2, y = value, fill = Var1)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_d3(palette = "category20c") +  theme(legend.position="right") + labs(fill = "Phylum") + xlab("Population") + ylab("Diversity - Genus Level within Phyla (%)")


#===================== For Paper at species IN genus Level ================

bang_genus <- kaiju_taxon_stat(bang_sep_taxon, genus)
mang_genus <- kaiju_taxon_stat(mang_sep_taxon, genus)
laksh_genus <- kaiju_taxon_stat(laksh_sep_taxon, genus)
ti_genus <- kaiju_taxon_stat(ti_sep_taxon, genus)
tii_genus <- kaiju_taxon_stat(tii_sep_taxon, genus)
tiii_phlyum <- kaiju_taxon_stat(tiii_sep_taxon, genus)
tiv_genus <- kaiju_taxon_stat(tiv_sep_taxon, genus)



#---------------------
#
# Abundance Plot
#
#----------------------

bang_taxon <- kaiju_groupby(bang,"Bacteria",species,10)
mang_taxon <- kaiju_groupby(mang,"Bacteria",species,10)
laksh_taxon <- kaiju_groupby(laksh,"Bacteria",species,10)
TI_taxon <- kaiju_groupby(TI,"Bacteria",species,10)
TII_taxon <- kaiju_groupby(TII,"Bacteria",species,10)
TIII_taxon <- kaiju_groupby(TIII,"Bacteria",species,10)
TIV_taxon <- kaiju_groupby(TIV,"Bacteria",species,10)

join_taxon <- data.frame(inner_join(bang_taxon,mang_taxon,by="species") %>% inner_join(laksh_taxon,by="species")%>% inner_join(TI_taxon,by="species") %>% inner_join(TII_taxon,by="species") %>% inner_join(TIII_taxon,by="species") %>% inner_join(TIV_taxon,by="species"))

colnames(join_taxon)<-c("taxon","Bang","Mang","Laskh","TI","TII","TIII","TIV")

join_taxon <- join_taxon %>% drop_na()
#for_genusAbundacePlot<- join_taxon[-c(3),]

data_percentage <- apply(as.data.frame(sapply(join_taxon[1:20,-1], as.numeric)), 2, function(x){x*100/sum(x,na.rm=T)})

rownames(data_percentage) <- join_taxon[1:20,1]

data = melt(data_percentage)

ggplot(data=data, aes(x = Var2, y = value, fill = Var1)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_igv()+  theme(legend.position="right") + labs(fill = "Species") + xlab("Population") + ylab("Abundance - Spceis Level within Phyla (%)")









#taxon_data <- separate(data = mang, col = taxon_name, into = c("phylum","class","order","family","genus","species"), sep = ";")
#taxon_data_filter <- delete.na(taxon_data, 3)
# taxon_data_filter <- as.data.frame(apply(taxon_data_filter,2,function(x)gsub('\\s+', '',(gsub("'", '',x)))))
#taxon_data_filter['reads'] <- as.data.frame(sapply(taxon_data_filter['reads'],as.integer))
# sapply(taxon_data_filter, class)
#return(as.data.frame(taxon_data_filter))
#taxon_data_filter %>% group_by(phylum) 
#%>%  summarise(reads, sum) %>% arrange(desc(reads)
