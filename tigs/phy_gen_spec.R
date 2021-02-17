#setwd("D:/IBAB_Data/tigs/new_data/pooled")
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(reshape)
library(ggplot2)
library(ggsci)
bang <- read.csv("Phy_Gen_Spc/Bangalore_All_kaiju.PGS.table",header = T,sep="\t", stringsAsFactors = FALSE)
mang <- read.csv("Phy_Gen_Spc/Mangalore_All_kaiju.PGS.table",header = T,sep="\t", stringsAsFactors = FALSE)
laksh <- read.csv("Phy_Gen_Spc/Lakshadweep_All_kaiju.PGS.table",header = T,sep="\t", stringsAsFactors = FALSE)
TI <- read.csv("Phy_Gen_Spc/TI_All_kaiju.PGS.table",header = T,sep="\t", stringsAsFactors = FALSE)
TII <- read.csv("Phy_Gen_Spc/TII_All_kaiju.PGS.table",header = T,sep="\t", stringsAsFactors = FALSE)
TIII <- read.csv("Phy_Gen_Spc/TIII_All_kaiju.PGS.table",header = T,sep="\t", stringsAsFactors = FALSE)
TIV <- read.csv("Phy_Gen_Spc/TIV_All_kaiju.PGS.table",header = T,sep="\t", stringsAsFactors = FALSE)

kaiju_separate_pgs <- function(kaiju_dataframe){
    taxon_data <- separate(data = kaiju_dataframe, col = taxon_name, into = c("phylum","genus","species"), sep = ";")
    return(taxon_data)
}
kaiju_separate_pgs_comma <- function(kaiju_dataframe){
  taxon_data <- separate(data = kaiju_dataframe, col = taxon_name, into = c("phylum","genus","species"), sep = ",")
  return(taxon_data)
}

kaiju_taxon_stat <- function(kaiju_dataframe, taxon_level){
  taxon_level <- enquo(taxon_level)
  kaiju_dataframe <- kaiju_dataframe %>% select(phylum,species) %>% unique()
  kaiju_taxa_count <- kaiju_dataframe %>% group_by(!!taxon_level) %>%  summarise(total=n()) %>% arrange(desc(total))
  return(kaiju_taxa_count)
}

bang_sep <- kaiju_separate_pgs(bang)
mang_sep <- kaiju_separate_pgs(mang)
laksh_sep <- kaiju_separate_pgs(laksh)
TI_sep <- kaiju_separate_pgs(TI)
TII_sep <- kaiju_separate_pgs_comma(TII)
TIII_sep <- kaiju_separate_pgs(TIII)
TIV_sep <- kaiju_separate_pgs(TIV)

bang_pgs <- bang_sep[,c(3,5,6,7)]
mang_pgs <- mang_sep[,c(3,5,6,7)]
laksh_pgs <- laksh_sep[,c(3,5,6,7)]
ti_pgs <- TI_sep[,c(3,5,6,7)]
tii_pgs <- TII_sep[,c(3,5,6,7)]
tiii_pgs <- TIII_sep[,c(3,5,6,7)]
tiv_pgs <- TIV_sep[,c(3,5,6,7)]

bang_alpha <- kaiju_taxon_stat(bang_pgs, phylum)
mang_alpha <- kaiju_taxon_stat(mang_pgs, phylum)
laksh_alpha <- kaiju_taxon_stat(laksh_pgs, phylum)
TI_alpha <- kaiju_taxon_stat(ti_pgs, phylum)
TII_alpha <- kaiju_taxon_stat(tii_pgs, phylum)
TIII_alpha <- kaiju_taxon_stat(tiii_pgs, phylum)
TIV_alpha <- kaiju_taxon_stat(tiv_pgs, phylum)

alpha_diversity <- as.data.frame(full_join(bang_alpha, mang_alpha, by ='phylum') %>% full_join(laksh_alpha, by='phylum')%>% 
 full_join(TI_alpha, by='phylum') %>% full_join(TII_alpha, by='phylum') %>% full_join(TIII_alpha, by='phylum') %>% full_join(TIV_alpha, by='phylum'))

colnames(alpha_diversity) <- c("Phylum","Bang","Mang","Laksh","TI","TII","TIII","TIV")

data = melt(alpha_diversity[1:20,c(1,2,3,5,6,7,8)])
head(data)
ggplot(data=data, aes(x = variable, y = value, fill = Phylum)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_igv() +  theme(legend.position="right") + labs(fill = "Phylum") + xlab("Population") + ylab("Diversity - Speceis Level within Phyla - Non Scaled")

#---------------------------------------- Scaling

data_percentage <- apply(as.data.frame(sapply(alpha_diversity[,-1], as.numeric)), 2, function(x){x*100/sum(x,na.rm=T)})

rownames(data_percentage) <- alpha_diversity[,1]

data = melt(data_percentage[1:20,c(1,2,4,5,6,7)])
head(data)
ggplot(data=data, aes(x = X2, y = value, fill = X1)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_igv() +  theme(legend.position="right") + labs(fill = "Phylum") + xlab("Population") + ylab("Diversity - Speceis Level within Phyla - Scaled")



head(bang_groups)
bang_groups <- bang_pgs %>% group_by(species) %>%  summarise(across(reads, sum)) %>% arrange(desc(reads))
mang_groups <- mang_pgs %>% group_by(species) %>%  summarise(across(reads, sum)) %>% arrange(desc(reads))
laksh_groups <- laksh_pgs %>% group_by(species) %>%  summarise(across(reads, sum)) %>% arrange(desc(reads))
TI_groups <- ti_pgs %>% group_by(species) %>%  summarise(across(reads, sum)) %>% arrange(desc(reads))
TII_groups <- tii_pgs %>% group_by(species) %>%  summarise(across(reads, sum)) %>% arrange(desc(reads))
TIII_groups <- tiii_pgs %>% group_by(species) %>%  summarise(across(reads, sum)) %>% arrange(desc(reads))
TIV_groups <- tiv_pgs %>% group_by(species) %>%  summarise(across(reads, sum)) %>% arrange(desc(reads))

  beta_diversity <- as.data.frame(full_join(bang_groups, mang_groups, by ='species') %>% full_join(laksh_groups, by='species')%>% 
                                     full_join(TI_groups, by='species') %>% full_join(TII_groups, by='species') %>% full_join(TIII_groups, by='species') %>% full_join(TIV_groups, by='species'))
  beta_diversity <- beta_diversity %>% drop_na()
  colnames(beta_diversity) <- c("Species","Bang","Mang","Laksh","TI","TII","TIII","TIV")
  
  data = melt(beta_diversity[1:20,c(1,2,3,5,6,7,8)])
  head(data)
  ggplot(data=data, aes(x = variable, y = value, fill = Species)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_igv() +  theme(legend.position="right") + labs(fill = "Species") + xlab("Population") + ylab("Abundance - Speceis Level within Phyla - Non Scaled")
  

#---------------------------------------- Scaling

data_percentage <- apply(as.data.frame(sapply(beta_diversity[1:20,-1], as.numeric)), 2, function(x){x*100/sum(x,na.rm=T)})

rownames(data_percentage) <- beta_diversity[1:20,1]

data = melt(data_percentage[1:20,c(1,2,4,5,6,7)])
head(data)
ggplot(data=data, aes(x = X2, y = value, fill = X1)) +  geom_bar(stat = "identity") + theme_bw() +  scale_fill_igv() +  theme(legend.position="right") + labs(fill = "Species") + xlab("Population") + ylab("Diversity - Speceis Level within Phyla - Scaled")

