setwd("D:/IBAB_Data/tigs/new_data/pooled")
source("../genus/kaju_function.R")

bang <- read.csv("../../new_refseq/species/Bangalore_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
mang <- read.csv("../../new_refseq/species/Mangalore_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
laksh <- read.csv("../../new_refseq/species/Lakshadweep_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
TI <- read.csv("../../new_refseq/species/TI_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
TII <- read.csv("../../new_refseq/species/TII_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
TIII <- read.csv("../../new_refseq/species/TIII_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
TIV <- read.csv("../../new_refseq/species/TIV_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)


bang_sep_taxon <- kaiju_separate(bang, 12)
mang_sep_taxon <- kaiju_separate(mang, 12)
laksh_sep_taxon <- kaiju_separate(laksh, 12)
ti_sep_taxon <- kaiju_separate(TI, 12)
tii_sep_taxon <- kaiju_separate(TII, 12)
tiii_sep_taxon <- kaiju_separate(TIII, 12)
tiv_sep_taxon <- kaiju_separate(TIV, 12)

bang_phylum <- kaiju_taxon_stat(bang_sep_taxon, phylum)
mang_phylum <- kaiju_taxon_stat(mang_sep_taxon, phylum)
laksh_phylum <- kaiju_taxon_stat(laksh_sep_taxon, phylum)
ti_phylum <- kaiju_taxon_stat(ti_sep_taxon, phylum)
tii_phylum <- kaiju_taxon_stat(tii_sep_taxon, phylum)
tiii_phlyum <- kaiju_taxon_stat(tiii_sep_taxon, phylum)
tiv_phylum <- kaiju_taxon_stat(tiv_sep_taxon, phylum)

write.table(bang_sep_taxon, "Bang_phylum_new_refseq.csv",quote = F, row.names = F,sep=",")
write.table(mang_sep_taxon, "Mang_phylum_new_refseq.csv",quote = F, row.names = F,sep=",")

phylum_data <- as.data.frame(full_join(bang_phylum, mang_phylum, by ='phylum') %>% full_join(laksh_phylum, by='phylum')  %>% 
  full_join(ti_phylum, by='phylum') %>% full_join(tii_phylum, by='phylum') %>% full_join(tiii_phlyum, by='phylum') %>% full_join(tiv_phylum, by='phylum'))
  

for_phylumGenusDiversityPlot <- phylum_data[1:10,]

for_phylumGenusDiversityPlot[11,] <- c("Other",as.numeric(colSums(phylum_data[11:nrow(phylum_data),2:ncol(phylum_data)], na.rm = TRUE)))


colnames(for_phylumGenusDiversityPlot) <- c("Phylum","Bang","Mang","Laksh","TI","TII","TIII","TIV")

data_percentage <- apply(as.data.frame(sapply(for_phylumGenusDiversityPlot[,-1], as.numeric)), 2, function(x){x*100/sum(x,na.rm=T)})

rownames(data_percentage) <- for_phylumGenusDiversityPlot[,1]

data = melt(data_percentage)

ggplot(data=data, aes(x = Var2, y = value, fill = Var1)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_d3(palette = "category20c") +  theme(legend.position="right") + labs(fill = "Phylum") + xlab("Population") + ylab("Diversity - Genus Level within Phyla (%)")


#---------------------
#
# Abundance Plot
#
#----------------------

bang_taxon <- kaiju_groupby(bang,"Bacteria",genus,12)
mang_taxon <- kaiju_groupby(mang,"Bacteria",genus,12)
laksh_taxon <- kaiju_groupby(laksh,"Bacteria",genus,12)
TI_taxon <- kaiju_groupby(TI,"Bacteria",genus,12)
TII_taxon <- kaiju_groupby(TII,"Bacteria",genus,12)
TIII_taxon <- kaiju_groupby(TIII,"Bacteria",genus,12)
TIV_taxon <- kaiju_groupby(TIV,"Bacteria",genus,12)

join_taxon <- data.frame(inner_join(bang_taxon,mang_taxon,by="genus") %>% inner_join(laksh_taxon,by="genus")%>% inner_join(TI_taxon,by="genus") %>% inner_join(TII_taxon,by="genus") %>% inner_join(TIII_taxon,by="genus") %>% inner_join(TIV_taxon,by="genus"))

colnames(join_taxon)<-c("taxon","Bang","Mang","Laskh","TI","TII","TIII","TIV")

join_taxon <- join_taxon %>% drop_na()
for_genusAbundacePlot<- join_taxon[-c(4),]

data_percentage <- apply(as.data.frame(sapply(for_genusAbundacePlot[1:20,-1], as.numeric)), 2, function(x){x*100/sum(x,na.rm=T)})

rownames(data_percentage) <- for_genusAbundacePlot[1:20,1]

data = melt(data_percentage)

ggplot(data=data, aes(x = Var2, y = value, fill = Var1)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_igv()+  theme(legend.position="right") + labs(fill = "Genus") + xlab("Population") + ylab("Abundance - Genus Level within Phyla (%)")

data

#---------------------
#
# Abundance Plot
#
#----------------------


bang_taxon <- kaiju_groupby(bang,"Bacteria",species,12)
mang_taxon <- kaiju_groupby(mang,"Bacteria",species,12)
laksh_taxon <- kaiju_groupby(laksh,"Bacteria",species,12)
TI_taxon <- kaiju_groupby(TI,"Bacteria",species,12)
TII_taxon <- kaiju_groupby(TII,"Bacteria",species,12)
TIII_taxon <- kaiju_groupby(TIII,"Bacteria",species,12)
TIV_taxon <- kaiju_groupby(TIV,"Bacteria",species,12)

join_taxon <- data.frame(inner_join(bang_taxon,mang_taxon,by="species") %>% inner_join(laksh_taxon,by="species")%>% inner_join(TI_taxon,by="species") %>% inner_join(TII_taxon,by="species") %>% inner_join(TIII_taxon,by="species") %>% inner_join(TIV_taxon,by="species"))

colnames(join_taxon)<-c("taxon","Bang","Mang","Laskh","TI","TII","TIII","TIV")

join_taxon <- join_taxon %>% drop_na()
for_genusAbundacePlot<- join_taxon[-c(3),]

data_percentage <- apply(as.data.frame(sapply(for_genusAbundacePlot[1:20,-1], as.numeric)), 2, function(x){x*100/sum(x,na.rm=T)})


data_percentage <- for_genusAbundacePlot[1:20,]
rownames(data_percentage) <- for_genusAbundacePlot[1:20,1]

data = melt(data_percentage,by='taxon')

ggplot(data=data, aes(x = variable, y = value, fill = taxon)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_igv()+  theme(legend.position="right") + labs(fill = "Genus") + xlab("Population") + ylab("Abundance - Species Level within Phyla (%)")


# nrow(bang_sep_taxon)
# nrow(bang_sep_taxon %>% select(phylum,genus) %>% unique())
# kaiju_taxa_count <- kaiju_dataframe %>% group_by(!!taxon_level) %>%  summarise(total=n()) %>% arrange(desc(total))
# return(kaiju_taxa_count)
# 



bang <- read.csv("../../new_refseq/species_filter/Bangalore_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
mang <- read.csv("../../new_refseq/species_filter/Mangalore_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
laksh <- read.csv("../../new_refseq/species_filter/Lakshadweep_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
TI <- read.csv("../../new_refseq/species_filter/TI_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
TII <- read.csv("../../new_refseq/species_filter/TII_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
TIII <- read.csv("../../new_refseq/species_filter/TIII_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)
TIV <- read.csv("../../new_refseq/species_filter/TIV_All_kaiju.All.table",header = T,sep="\t", stringsAsFactors = FALSE)

bang_sep_taxon <- kaiju_separate(bang, 12)
mang_sep_taxon <- kaiju_separate(mang, 12)
laksh_sep_taxon <- kaiju_separate(laksh, 12)
ti_sep_taxon <- kaiju_separate(TI, 12)
tii_sep_taxon <- kaiju_separate(TII, 12)
tiii_sep_taxon <- kaiju_separate(TIII, 12)
tiv_sep_taxon <- kaiju_separate(TIV, 12)

bang_phylum <- kaiju_taxon_stat(bang_sep_taxon, phylum)
mang_phylum <- kaiju_taxon_stat(mang_sep_taxon, phylum)
laksh_phylum <- kaiju_taxon_stat(laksh_sep_taxon, phylum)
ti_phylum <- kaiju_taxon_stat(ti_sep_taxon, phylum)
tii_phylum <- kaiju_taxon_stat(tii_sep_taxon, phylum)
tiii_phlyum <- kaiju_taxon_stat(tiii_sep_taxon, phylum)
tiv_phylum <- kaiju_taxon_stat(tiv_sep_taxon, phylum)

phylum_data <- as.data.frame(full_join(bang_phylum, mang_phylum, by ='phylum') %>% full_join(laksh_phylum, by='phylum')  %>% full_join(ti_phylum, by='phylum') %>% full_join(tii_phylum, by='phylum') %>% full_join(tiii_phlyum, by='phylum') %>% full_join(tiv_phylum, by='phylum'))


for_phylumGenusDiversityPlot <- phylum_data[1:10,]

for_phylumGenusDiversityPlot[11,] <- c("Other",as.numeric(colSums(phylum_data[11:nrow(phylum_data),2:ncol(phylum_data)], na.rm = TRUE)))


colnames(for_phylumGenusDiversityPlot) <- c("Phylum","Bang","Mang","Laksh","TI","TII","TIII","TIV")

data_percentage <- apply(as.data.frame(sapply(for_phylumGenusDiversityPlot[,-1], as.numeric)), 2, function(x){x*100/sum(x,na.rm=T)})

rownames(data_percentage) <- for_phylumGenusDiversityPlot[,1]

data = melt(data_percentage)

ggplot(data=data, aes(x = variable, y = value, fill = taxon)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_igv()+  theme(legend.position="right") + labs(fill = "Phylum") + xlab("Population") + ylab("Diversity - Genus Level within Phyla (%)")


