library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(ggsci)

delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}

kaiju_separate <- function(kaiju_dataframe, numCols){
  if(numCols == 12){
    taxon_data <- separate(data = kaiju_dataframe, col = taxon_name, into = c("Organism","kingdom","phylum","class","order","family","genus","species"), sep = ";")
  }else{
    taxon_data <- separate(data = kaiju_dataframe, col = taxon_name, into = c("phylum","class","order","family","genus","species"), sep = ";")
  }

 # taxon_data_filter <- delete.na(taxon_data, 3)
  return(taxon_data)
}

# kaiju_separate <- function(kaiju_dataframe){
#   taxon_data <- separate(data = kaiju_dataframe, col = taxon_name, into = c("Organisms","Kingdom","phylum","class","order","family","genus","species"), sep = ";")
#   #taxon_data_filter <- delete.na(taxon_data, 3)
#   return(taxon_data_filter)
# }

kaiju_groupby <- function(kaiju_dataframe, master_group, taxon_level, numCols){
  kaiju_dataframe <- kaiju_separate(kaiju_dataframe, numCols)
  taxon_level <- enquo(taxon_level)
  print(ncol(kaiju_dataframe))
  if(ncol(kaiju_dataframe) == 12){
    kaiju_master_group <- kaiju_dataframe[kaiju_dataframe$kingdom==master_group,]
    kaiju_groups <- kaiju_master_group %>% group_by(!!taxon_level) %>%  summarise(across(reads, sum)) %>% arrange(desc(reads))
  }
  else{
    kaiju_groups <- kaiju_dataframe %>% group_by(!!taxon_level) %>%  summarise(across(reads, sum)) %>% arrange(desc(reads))
  }
  return(kaiju_groups)
}

kaiju_join <- function(join_type,df1,df2,df3,df4,df5,df6,df7){
  print(paste0("Doing '",join_type,"' For givendataframe"),quote = F)
  warning("Replacing 'NA' with 1",call. = TRUE)
  if(join_type == "inner"){
    join_phylum <- inner_join(df1,df2,by="phylum") %>% inner_join(df3,by="phylum")%>% 
      inner_join(df4,by="phylum") %>% inner_join(df5,by="phylum") %>% inner_join(df6,by="phylum") %>% 
      inner_join(df7,by="phylum")
  } else if(join_type == "full"){
    join_phylum <- full_join(df1,df2,by="phylum") %>% full_join(df3,by="phylum")%>% 
      full_join(df4,by="phylum") %>% full_join(df5,by="phylum") %>% full_join(df6,by="phylum") %>% 
      full_join(df7,by="phylum")
  }else if(join_type == "left"){
    join_phylum <- left_join(df1,df2,by="phylum") %>% left_join(df3,by="phylum")%>% 
      left_join(df4,by="phylum") %>% left_join(df5,by="phylum") %>% left_join(df6,by="phylum") %>% 
      left_join(df7,by="phylum")
  }else if(join_type == "right"){
    join_phylum <- right_join(df1,df2,by="phylum") %>% right_join(df3,by="phylum")%>% 
      right_join(df4,by="phylum") %>% right_join(df5,by="phylum") %>% right_join(df6,by="phylum") %>% 
      right_join(df7,by="phylum")
  }else{
    print("Joint not found",quote = F)
    return(join_phylum=NULL)
  }
  join_phylum[is.na(join_phylum)] <- 1
  return(join_phylum)
}

kaiju_diversity <- function(kaiju_dataframe){
  print("Separting Taxon Column for extracting the Taxon Inforamtion",quote = F)
  kaiju_dataframe <- kaiju_separate(kaiju_dataframe)
  print("Txon Separtion done",quote = F)
  print(paste0("Number of Line in given dataframe are: ",nrow(taxon_data)),quote = F)
  print(paste0("Number of lines after removing bottom 3 rows : " , nrow(head(taxon_data, -3))),quote = F)
  taxon_data <- head(taxon_data, -3)
  print("Filtering data based on number of reads ( Number of Reads > 50)",quote = F)
  taxon_data <- taxon_data[taxon_data$reads>50,]
  print(paste0("Number of rows remaining after filtering : ",nrow(taxon_data)),quote = F)
  print("Selecting Phlyum + Genus column for removing the redundancy")
  unique_taxon <- unique(taxon_data[,c(5,9)])
  print("Counting number of genus present in given phlyum, and transposing the matrix")
  phylum_genus_df <- as.data.frame(t(unique_taxon %>%  pivot_wider(names_from = phylum, values_from = genus)))
  print("You have phlyum_genus_df which contains phylum Respective genus names")
  print("Creating Final Dataframe containign Genus Coutn + Phlyum Name")
  final_df <- as.data.frame(lengths(phylum_genus_df$V1))
  final_df$phylum <- rownames(final_df)
  print("Coulumn naming & rearrenging")
  colnames(final_df)<-c("GenusCount","Phlyum")
  final_df <- final_df[,c(2,1)]
  print("Sorting dataframe By Geneus Count")
  final_df <-final_df[order(-final_df$GenusCount),]
  return(final_df)
}

kaiju_taxon_stat <- function(kaiju_dataframe, taxon_level){
  taxon_level <- enquo(taxon_level)
  kaiju_dataframe <- kaiju_dataframe %>% select(phylum,genus) %>% unique()
  kaiju_taxa_count <- kaiju_dataframe %>% group_by(!!taxon_level) %>%  summarise(total=n()) %>% arrange(desc(total))
  return(kaiju_taxa_count)
}
