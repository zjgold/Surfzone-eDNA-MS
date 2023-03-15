# ---
# title: "Generalized Decontamination Script"
# author: "Zack Gold"
# date: "1/02/2023"
# output: html_document
# ---

# Load Libraries ----------------------------------------------------------------------------------------------------------------
library (tidyverse)
library (vegan)
library (proxy)
library(reshape2)
library(microDecon)
library(stringr)
library(knitr)
library(ranacapa)
library(dplyr)
library(tibble)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(plotly)
library(optparse)
library(fitdistrplus)
library(rstan)
library(shinystan)
library(bayesplot)
library(broom)
library(analyze.stuff)
library(here)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("/Users/zgold/Documents/UCLA/Surfzone/2018_analysis/analysis/decontam/decontamination_utilities.R")

#Test inputs ----------------------------------------------------------------------------------------------------------------
working_directory <- here("decontam")
input_anacapa_path_1 <- here("analysis","data","anacapa_20230101","ucsb_2018_elas_combo_20230101","c19_fishcard_taxonomy_tables","Summary_by_percent_confidence","60","c19_fishcard_ASV_raw_taxonomy_60.txt")
input_anacapa_path_2 <- here("analysis","data","anacapa_20230101","ucsb_2018_elas_fishcard_20230101","fishcard_12S_all_taxonomy_tables","Summary_by_percent_confidence","60","fishcard_12S_all_ASV_raw_taxonomy_60.txt")
input_anacapa_path_3 <- here("analysis","data","anacapa_20230101","ucsb_2018_miu_combo_20230101","c19_fishcard_taxonomy_tables","Summary_by_percent_confidence","60","c19_fishcard_ASV_raw_taxonomy_60.txt")
input_anacapa_path_4 <- here("analysis","data","anacapa_20230101","ucsb_2018_miu_fishcard_20230101","fishcard_12S_all_taxonomy_tables","Summary_by_percent_confidence","60","fishcard_12S_all_ASV_raw_taxonomy_60.txt")

input_meta_path <- here("analysis","data","metadata_ucsb_2018_032020.txt")
number_anacapa_tables <- "two"

read_type <- "merged_only"

low_dist_probability_cutoff <- 0.05
minimum_read_cutoff <- 1000
high_dist_probability_cutoff <- 1

som_level <- "double"
som_chains<-10
som_iterations<-1000
som_filter_choice <- "som_min_threshold_2"

max_dist_probability_cutoff <- 0.9
#---

#---

#Load Data ----------------------------------------------------------------------------------------------------------------
dir.create(here("analysis","decontam","Output_csv"))
dir.create(here("analysis","decontam","Output_R"))
dir.create(here("analysis","decontam","Output_plots"))

metadata <- read.table(input_meta_path, header = 1, sep = "\t", stringsAsFactors = F)

anacapa_table_1 <- read.table(input_anacapa_path_1, header = 1, sep = "\t", stringsAsFactors = F)
anacapa_table_2 <- read.table(input_anacapa_path_2, header = 1, sep = "\t", stringsAsFactors = F)
anacapa_table_3 <- read.table(input_anacapa_path_3, header = 1, sep = "\t", stringsAsFactors = F)
anacapa_table_4 <- read.table(input_anacapa_path_4, header = 1, sep = "\t", stringsAsFactors = F)

#Format Anacapa Table 1
anacapa_table_1[1,1] %>% str_remove(.,"forward_") %>% str_sub(., end=-3) -> barcode_1
barcode_1_use <- paste0("elas_",barcode_1)
anacapa_table_1$Miseq_run <- barcode_1_use
  
###Fix Names
anacapa_table_1_names <- colnames(anacapa_table_1)
anacapa_table_1_names %>% as.data.frame() -> anacapa_table_1_names
anacapa_table_1_names[1,1] <- "seq_number"
colnames(anacapa_table_1_names) <- c("Seq_number")
anacapa_table_1_names %>% mutate(.,  Seq_number =str_remove(Seq_number,paste0(barcode_1,"_"))) -> anacapa_table_1_names
  
left_join(anacapa_table_1_names, metadata) %>% dplyr::select(Seq_number, New_name) %>% 
    mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_1_names
colnames(anacapa_table_1) <- anacapa_table_1_names$final_names
anacapa_table_1 %>% 
  mutate(., seq_number = str_replace(seq_number, "c19_fishcard", barcode_1_use)) -> anacapa_table_1

#Format Anacapa Table 2
anacapa_table_2[1,1] %>% str_remove(.,"forward_") %>% str_sub(., end=-3) -> barcode_2
barcode_2_use <- paste0("elas_",barcode_2)
anacapa_table_2$Miseq_run <- barcode_2_use
  
###Fix Names
anacapa_table_2_names <- colnames(anacapa_table_2)
anacapa_table_2_names %>% as.data.frame() -> anacapa_table_2_names
anacapa_table_2_names[1,1] <- "seq_number"
colnames(anacapa_table_2_names) <- c("Seq_number")
anacapa_table_2_names %>% mutate(.,  Seq_number =str_remove(Seq_number,paste0(barcode_2,"_"))) -> anacapa_table_2_names
  
left_join(anacapa_table_2_names, metadata) %>% dplyr::select(Seq_number, New_name) %>% 
    mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_2_names
colnames(anacapa_table_2) <- anacapa_table_2_names$final_names
anacapa_table_2 %>% 
  mutate(., seq_number = str_replace(seq_number, "fishcard_12S_all", barcode_2_use)) -> anacapa_table_2

#Format Anacapa Table 3
anacapa_table_3[1,1] %>% str_remove(.,"forward_") %>% str_sub(., end=-3) -> barcode_3
barcode_3_use <- paste0("miu_",barcode_3)
anacapa_table_3$Miseq_run <- barcode_3_use
  
###Fix Names
anacapa_table_3_names <- colnames(anacapa_table_3)
anacapa_table_3_names %>% as.data.frame() -> anacapa_table_3_names
anacapa_table_3_names[1,1] <- "seq_number"
colnames(anacapa_table_3_names) <- c("Seq_number")
anacapa_table_3_names %>% mutate(.,  Seq_number =str_remove(Seq_number,paste0(barcode_3,"_"))) -> anacapa_table_3_names
  
left_join(anacapa_table_3_names, metadata) %>% dplyr::select(Seq_number, New_name) %>% 
    mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_3_names
colnames(anacapa_table_3) <- anacapa_table_3_names$final_names
anacapa_table_3 %>% 
  mutate(., seq_number = str_replace(seq_number, "c19_fishcard", barcode_3_use)) -> anacapa_table_3

#Format Anacapa Table 4
anacapa_table_4[1,1] %>% str_remove(.,"forward_") %>% str_sub(., end=-3) -> barcode_4
barcode_4_use <- paste0("miu_",barcode_4)
anacapa_table_4$Miseq_run <- barcode_4_use

###Fix Names
anacapa_table_4_names <- colnames(anacapa_table_4)
anacapa_table_4_names %>% as.data.frame() -> anacapa_table_4_names
anacapa_table_4_names[1,1] <- "seq_number"
colnames(anacapa_table_4_names) <- c("Seq_number")
anacapa_table_4_names %>% mutate(.,  Seq_number =str_remove(Seq_number,paste0(barcode_4,"_"))) -> anacapa_table_4_names

left_join(anacapa_table_4_names, metadata) %>% dplyr::select(Seq_number, New_name) %>% 
  mutate(final_names = coalesce(New_name, Seq_number)) %>%  dplyr::select(final_names) %>% as.list() -> anacapa_table_4_names
colnames(anacapa_table_4) <- anacapa_table_4_names$final_names
anacapa_table_4 %>% 
  mutate(., seq_number = str_replace(seq_number, "fishcard_12S_all", barcode_4_use)) -> anacapa_table_4

#Generate Hash.key
anacapa_table_1 %>% 
    dplyr::select(seq_number, sum.taxonomy) %>% 
    filter(., str_detect(seq_number, "merged")) %>% 
    filter(., !str_detect(seq_number, "unmerged")) %>% 
  dplyr::select(sum.taxonomy.c19= sum.taxonomy)-> hash.key_1_c19
  head(hash.key_1_c19)
anacapa_table_2 %>% 
    dplyr::select(seq_number, sum.taxonomy) %>% 
    filter(., str_detect(seq_number, "merged")) %>% 
    filter(., !str_detect(seq_number, "unmerged")) ->  hash.key_1_fish
head(hash.key_1_fish)
anacapa_table_3 %>% 
    filter(., str_detect(seq_number, "merged")) %>% 
    filter(., !str_detect(seq_number, "unmerged")) %>% 
    dplyr::select(sum.taxonomy.c19= sum.taxonomy)  -> hash.key_2_c19
head(hash.key_2_c19)
anacapa_table_4 %>% 
    filter(., str_detect(seq_number, "merged")) %>% 
    filter(., !str_detect(seq_number, "unmerged"))  -> hash.key_2_fish
head(hash.key_2_fish)
  dim(hash.key_1_c19)
  dim(hash.key_1_fish)
  
cbind(hash.key_1_fish, hash.key_1_c19) %>% 
    mutate(., mammal = str_detect(sum.taxonomy.c19,"Mammalia")) %>% 
    mutate(., bird = str_detect(sum.taxonomy.c19,"Aves")) %>% 
    mutate_if(is.character, list(~na_if(.,""))) %>%  
    mutate(., empty = is.na(sum.taxonomy)) %>% 
    mutate(., empty_2 = is.na(sum.taxonomy.c19)) %>% 
    mutate(., sum.taxonomy_final = case_when(mammal==TRUE ~ sum.taxonomy.c19,
                                             bird==TRUE ~ sum.taxonomy.c19,
                                             sum.taxonomy == "Eukaryota;Chordata;Actinopteri;;Embiotocidae;Amphistichus;Amphistichus argenteus" ~ "Eukaryota;Chordata;Actinopteri;;Embiotocidae;Amphistichus;",
                                             sum.taxonomy == "Eukaryota;Chordata;Actinopteri;;Embiotocidae;Amphistichus;Amphistichus koelzi" ~ "Eukaryota;Chordata;Actinopteri;;Embiotocidae;Amphistichus;",
                                             sum.taxonomy == "Eukaryota;Chordata;Actinopteri;Atheriniformes;Atherinopsidae;Atherinopsis;Atherinopsis californiensis" ~ "Eukaryota;Chordata;Actinopteri;Atheriniformes;Atherinopsidae;Atherinops;Atherinops affinis",
                                             sum.taxonomy == "Eukaryota;Chordata;Actinopteri;Syngnathiformes;Syngnathidae;Syngnathus;Syngnathus leptorhynchus" ~ "Eukaryota;Chordata;Actinopteri;Syngnathiformes;Syngnathidae;Syngnathus;",
                                             (empty==TRUE & empty_2 ==FALSE) ~ sum.taxonomy.c19,
                                             (empty==FALSE & empty_2 ==TRUE) ~ sum.taxonomy,
                                             (empty==FALSE & empty_2 ==FALSE) ~ sum.taxonomy,
                                             sum.taxonomy == "Eukaryota;Chordata;Actinopteri;;;;" & sum.taxonomy.c19 != "Eukaryota;Chordata;Actinopteri;;;;" & sum.taxonomy.c19 != "" ~ sum.taxonomy.c19,
                                             sum.taxonomy == "Eukaryota;Chordata;;;;;" & sum.taxonomy.c19 != "Eukaryota;Chordata;;;;;" ~ sum.taxonomy.c19,
                                             TRUE ~ sum.taxonomy)) %>% 
    dplyr::select(seq_number,sum.taxonomy =sum.taxonomy_final) -> hash.key_1_updated

cbind(hash.key_2_fish, hash.key_2_c19) %>% 
  mutate(., mammal = str_detect(sum.taxonomy.c19,"Mammalia")) %>% 
  mutate(., bird = str_detect(sum.taxonomy.c19,"Aves")) %>% 
  mutate_if(is.character, list(~na_if(.,""))) %>%  
  mutate(., empty = is.na(sum.taxonomy)) %>% 
  mutate(., empty_2 = is.na(sum.taxonomy.c19)) %>% 
  mutate(., sum.taxonomy_final = case_when(mammal==TRUE ~ sum.taxonomy.c19,
                                           bird==TRUE ~ sum.taxonomy.c19,
                                           sum.taxonomy == "Eukaryota;Chordata;Actinopteri;;Embiotocidae;Amphistichus;Amphistichus argenteus" ~ "Eukaryota;Chordata;Actinopteri;;Embiotocidae;Amphistichus;",
                                           sum.taxonomy == "Eukaryota;Chordata;Actinopteri;;Embiotocidae;Amphistichus;Amphistichus koelzi" ~ "Eukaryota;Chordata;Actinopteri;;Embiotocidae;Amphistichus;",
                                           sum.taxonomy == "Eukaryota;Chordata;Actinopteri;Atheriniformes;Atherinopsidae;Atherinopsis;Atherinopsis californiensis" ~ "Eukaryota;Chordata;Actinopteri;Atheriniformes;Atherinopsidae;Atherinops;Atherinops affinis",
                                           sum.taxonomy == "Eukaryota;Chordata;Actinopteri;Syngnathiformes;Syngnathidae;Syngnathus;Syngnathus leptorhynchus" ~ "Eukaryota;Chordata;Actinopteri;Syngnathiformes;Syngnathidae;Syngnathus;",
                                           (empty==TRUE & empty_2 ==FALSE) ~ sum.taxonomy.c19,
                                           (empty==FALSE & empty_2 ==TRUE) ~ sum.taxonomy,
                                           (empty==FALSE & empty_2 ==FALSE) ~ sum.taxonomy,
                                           sum.taxonomy == "Eukaryota;Chordata;Actinopteri;;;;" & sum.taxonomy.c19 != "Eukaryota;Chordata;Actinopteri;;;;" & sum.taxonomy.c19 != "" ~ sum.taxonomy.c19,
                                           sum.taxonomy == "Eukaryota;Chordata;;;;;" & sum.taxonomy.c19 != "Eukaryota;Chordata;;;;;" ~ sum.taxonomy.c19,
                                           TRUE ~ sum.taxonomy)) %>% 
  dplyr::select(seq_number,sum.taxonomy =sum.taxonomy_final) -> hash.key_2_updated

rbind(hash.key_1_updated,hash.key_2_updated) -> hash.key

saveRDS(hash.key,file= here("analysis","decontam","Output_R","hash.key.RDS"))
saveRDS(barcode_4_use,file= here("analysis","decontam","Output_R","barcode_4_use.RDS"))
saveRDS(barcode_2_use,file= here("analysis","decontam","Output_R","barcode_2_use.RDS"))

# Convert to Long Data

###Format for Long Data
anacapa_table_2$seq_number <- factor(anacapa_table_2$seq_number)
anacapa_table_2$Miseq_run <- factor(anacapa_table_2$Miseq_run)

columns <- colnames(anacapa_table_2)
remove <- c("seq_number","sum.taxonomy","Miseq_run")

gathercols <-  columns[! columns %in% remove] 

anacapa_table_2 %>% 
  pivot_longer(., cols=gathercols, names_to="sample", values_to="reads") -> anacapa_table_2_long

anacapa_table_2_long$reads <- as.numeric(anacapa_table_2_long$reads)

###Format for Long Data
anacapa_table_4$seq_number <- factor(anacapa_table_4$seq_number)
anacapa_table_4$Miseq_run <- factor(anacapa_table_4$Miseq_run)

columns <- colnames(anacapa_table_4)
remove <- c("seq_number","sum.taxonomy","Miseq_run")

gathercols <-  columns[! columns %in% remove] 

anacapa_table_4 %>% 
  pivot_longer(., cols=gathercols, names_to="sample", values_to="reads") -> anacapa_table_4_long

anacapa_table_4_long$reads <- as.numeric(anacapa_table_4_long$reads)
#Merge All Tables
ASV.table_long <- rbind(anacapa_table_2_long,anacapa_table_4_long)

ASV.table_long %>% 
  dplyr::select(-sum.taxonomy) %>% 
  left_join(hash.key) -> ASV.table_long
 
#---

#---

# Cleaning Process 0: Remove all Forward, Reverse, and Unmerged reads & Remove Singletons ----------------------------------------------------------------------------------------------------------------

#Filter Merged only reads

#Filter out forward, reverse, and unmerged reads
ASV.table_long %>% 
    filter(., str_detect(seq_number,"merged")) %>% 
    filter(., !str_detect(seq_number,"unmerged")) -> ASV.table_used

#Calculate % ASVs Kept
ASV.table_long %>%  dim() -> all_dim
ASV.table_used %>%  dim() -> used_only_dim

paste0(round(used_only_dim[[1]]/all_dim[[1]]*100,2),"% ASVs retained")

saveRDS(ASV.table_long, file=here("analysis","decontam","ASV.table_long.RDS"))
saveRDS(ASV.table_used, file=here("analysis","decontam","ASV.table_used.RDS"))

ASV.table_long %>% 
  group_by(Miseq_run) %>% 
  dplyr::summarise(n_distinct(seq_number), sum(reads))

ASV.table_used %>% 
  group_by(Miseq_run) %>% 
  dplyr::summarise(n_distinct(seq_number), sum(reads))

# Technical Replicate Exploration -------
ASV.table_used <- readRDS(file=here("analysis","decontam","ASV.table_used.RDS"))
saveRDS(hash.key,file= here("analysis","decontam","Output_R","hash.key.RDS"))
saveRDS(barcode_4_use,file= here("analysis","decontam","Output_R","barcode_4_use.RDS"))
saveRDS(barcode_2_use,file= here("analysis","decontam","Output_R","barcode_2_use.RDS"))

ASV.table_used %>% 
left_join(metadata, by=c("sample"="New_name")) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarise(sample_read_depth= sum(reads)) -> sample_read_depth

ASV.table_used %>% 
  filter(., !str_detect(seq_number,"elas")) %>% 
  left_join(metadata, by=c("sample"="New_name")) %>% 
  left_join(sample_read_depth) %>% 
  dplyr::group_by(Site, Bio_rep,seq_number) %>% 
  dplyr::summarise(mean_reads = mean(reads), 
                   sd_reads=sd(reads), 
                   max_reads=max(reads),
                   min_reads=min(reads),
                   mean_prop = mean(reads/sample_read_depth), 
                   sd_prop=sd(reads/sample_read_depth), 
                   max_prop=max(reads/sample_read_depth),
                   min_prop=min(reads/sample_read_depth),
                   mean_read_depth=mean(sample_read_depth)) %>% 
  ungroup() %>% 
  group_by(seq_number) %>%
  filter(., max_prop >0.01) %>% 
  filter(., sd_prop > 0) %>% 
  mutate( cv = sd_prop/mean_prop) -> ASV_prop_data


ASV_prop_data %>% 
  ggplot(., aes(x=Bio_rep, y= cv)) +  geom_boxplot() + facet_wrap(~Site)

ASV_prop_data %>% 
  tidyr::unite(., "Bio_rep",c("Site","Bio_rep")) %>% 
  ggplot(., aes(x=cv)) +  geom_histogram() + facet_wrap(~Bio_rep)

ASV_prop_data %>% 
  tidyr::unite(., "Bio_rep",c("Site","Bio_rep")) %>% 
  ggplot(., aes(x=cv)) +  geom_histogram() 

ASV_prop_data %>% 
  ggplot(., aes(x=mean_read_depth, y=mean_prop, color=seq_number)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90))+ 
  geom_errorbar(aes(ymin=mean_prop-sd_prop, ymax=mean_prop+sd_prop), width=.2,position=position_dodge(.9)) +
  theme(legend.position = "none")

ASV_prop_data %>% 
  ggplot(., aes(x=mean_read_depth, y=cv, color=seq_number)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90))+ 
  theme(legend.position = "none")

#---

#---


ASV.table_used <- readRDS(file=here("analysis","decontam","ASV.table_used.RDS"))
hash.key <- readRDS(file= here("analysis","decontam","Output_R","hash.key.RDS"))
barcode_4_use <- readRDS(file= here("analysis","decontam","Output_R","barcode_4_use.RDS"))
barcode_2_use <-readRDS(file= here("analysis","decontam","Output_R","barcode_2_use.RDS"))
metadata <- read.table(input_meta_path, header = 1, sep = "\t", stringsAsFactors = F)



# Cleaning Process 1: Estimation of *Tag-jumping* or sample *Cross-talk* ----------------------------------------------------------------------------------------------------------------

## Step 1: Nest the dataset by origin of ASVs

###Identify Positives, Negatives, and Samples

###Create list of control samples
metadata %>% 
  filter(Sample_Control=="Control") %>% 
  dplyr::select(New_name) %>% unique() -> controls
controls <- controls$New_name

metadata %>% 
  filter(Control_Type=="Pos") %>% 
  dplyr::select(New_name) %>% unique()-> pos_controls
pos_controls <- pos_controls$New_name

metadata %>% 
  filter(Control_Type=="Blank") %>% 
  dplyr::select(New_name) %>% unique() -> neg_controls
neg_controls <- neg_controls$New_name

metadata %>% 
  filter(Location=="Gold_mpawatch") %>% 
  dplyr::select(New_name) %>% unique() -> Gold_mpawatch
Gold_mpawatch <- Gold_mpawatch$New_name


###New column that labels each ASV as from Positive (control) or Sample
ASV.table_used %>% 
  mutate(source = case_when(sample %in% pos_controls~"Positives",
                            sample %in% neg_controls~"Blanks",
                            TRUE ~"Samples")) -> ASV.table_used


###Convert to tibble
ASV.table_used <- as_tibble(ASV.table_used)

###Remove empty sequences
ASV.table_used %>% 
  filter(reads != 0)  -> ASV.table_used

###Rename Columns and remove seq_number
ASV.table_used %>%
  mutate(sample = as.character(sample),
         nReads = reads)  -> ASV.table_used

###ASVs in Positive Controls
ASV.table_used %>% 
  filter (source == "Positives") %>%
  dplyr::group_by(seq_number) %>% 
  dplyr::summarise(tot = sum(reads)) %>% 
  arrange(desc(tot)) %>% 
  pull(seq_number) -> all.seqs.in.positives

hash.key %>% 
  filter(seq_number %in% all.seqs.in.positives) %>% as_tibble() -> pos.contam.species

write.csv(pos.contam.species, file=here("analysis","decontam","pos.contam.species.csv"))

###ASVs in Negative Controls

ASV.table_used %>% 
  filter (source == "Blanks") %>%
  dplyr::group_by(seq_number) %>% 
  dplyr::summarise(tot = sum(reads)) %>% 
  arrange(desc(tot)) %>% 
  pull(seq_number) -> all.seqs.in.blanks

hash.key %>% 
  filter(seq_number %in% all.seqs.in.blanks) %>% as_tibble() -> blank.contam.species

write.csv(blank.contam.species, file=here("analysis","decontam","blank.contam.species.csv"))

### Visualize Read Counts Across Samples for Barcode_1
ASV.table_used %>% 
  group_by(sample) %>%
  filter(., Miseq_run== barcode_2_use) %>% 
  mutate (TotalReadsperSample = sum(nReads)) %>%
  arrange(desc(TotalReadsperSample)) %>%
  ggplot(., aes(x=sample, y=TotalReadsperSample, color=source)) + 
  geom_point() +ggtitle("Read Count Across Samples") + 
  theme(axis.text.x = element_text(angle = 90)) -> plot_1

ggsave(plot=plot_1,here("analysis","decontam","Output_plots","Sample_Read_Depth_Barcode_1.png"), device = "png", width = 12, height = 8, units = "in")

### Visualize Read Counts Across Samples for Barcode_2

ASV.table_used %>% 
    group_by(sample) %>%
    filter(., Miseq_run== barcode_4_use) %>% 
    mutate (TotalReadsperSample = sum(nReads)) %>%
    arrange(desc(TotalReadsperSample)) %>%
    ggplot(., aes(x=sample, y=TotalReadsperSample, color=source)) + 
    geom_point() +ggtitle("Read Count Across Samples") + 
    theme(axis.text.x = element_text(angle = 90)) -> plot_2
  
ggsave(plot=plot_2,here("analysis","decontam","Output_plots","Sample_Read_Depth_Barcode_2.png"), device = "png", width = 12, height = 8, units = "in")

  
###Nesting the dataset
ASV.table_used %>% 
  dplyr::group_by(Miseq_run, source) %>% 
  nest() %>% 
  pivot_wider(names_from=source, values_from=data) -> ASV.nested


####Summary.file.1
ASV.nested %>% 
  ungroup() %>% 
  dplyr::transmute(.,Miseq_run,Summary = purrr::map(Samples, ~ how.many(ASVtable = ., round = 0)))  -> ASV.summary

ASV.summary$Summary

## Step 2: Model the composition of the positive controls of each run 

###Jumping vector

ASV.nested %>% 
  mutate (contam.tibble = purrr::map(Positives, 
                                     function(.x){
                                       .x %>%
                                         ungroup() %>% 
                                         group_by(sample) %>%
                                         mutate (TotalReadsperSample = sum(nReads)) %>%
                                         mutate (proportion = nReads/TotalReadsperSample) %>%
                                         group_by(seq_number) %>%
                                         dplyr::summarise (vector_contamination = max(proportion))
                                     }) ) -> ASV.nested

###Vector Contamination in Barcode_1
ASV.nested$contam.tibble[[1]] %>% as.data.frame() %>% 
  ggplot(aes(x= vector_contamination))+
  geom_histogram() -> vc_plot_1 # Check how it looks
ggsave(plot=vc_plot_1,here("analysis","decontam","Output_plots","Vector_Contamination_Barcode_1.png"), device = "png", width = 12, height = 8, units = "in")

###Vector Contamination in Barcode_2

ASV.nested$contam.tibble[[2]] %>% as.data.frame() %>% 
    ggplot(aes(x= vector_contamination))+
    geom_histogram() -> vc_plot_2 # Check how it looks
ggsave(plot=vc_plot_2,here("analysis","decontam","Output_plots","Vector_Contamination_Barcode_2.png"), device = "png", width = 12, height = 8, units = "in")


##Step 3: Substract the composition of the positive controls from the environment samples

ASV.nested %>% 
  ungroup() %>% 
  mutate(Step1.cleaned.tibble = map2(Samples, contam.tibble, function(.x,.y){ 
    .x %>%
      dplyr::group_by (sample) %>%
      mutate (TotalReadsperSample = sum (nReads)) %>%
      left_join(.y, by = "seq_number") %>%
      mutate (Updated_nReads = ifelse (!is.na(vector_contamination),  nReads - (ceiling(vector_contamination*TotalReadsperSample)), nReads)) %>%
      filter (Updated_nReads > 0) %>%
      ungroup() %>% 
      dplyr::select (sample, seq_number, nReads = Updated_nReads)
  })) -> ASV.nested

###Add this step to the summary table we were creating

####Summary.file.2
ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step1.cleaned.tibble, ~ how.many(ASVtable = .,round = "1.Jump"))) %>% 
  left_join(ASV.summary) %>% #use left join when there are many miseq runs to join
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary 

ASV.summary$Summary

#---

#---


# Cleaning Process 2: **Discarding PCR replicates with low number of reads** ----------------------------------------------------------------------------------------------------------------

###Pull out Sample Read Depth

ASV.nested %>% 
  dplyr::select(Miseq_run, Step1.cleaned.tibble) %>% 
  unnest(Step1.cleaned.tibble) %>% 
  group_by(Miseq_run,sample) %>%
  dplyr::summarise(tot = sum(nReads)) %>% 
  arrange(desc(tot))-> all.reps

all.reps %>%  
  filter(., Miseq_run==barcode_2_use) ->all.reps_1

all.reps %>%  
  filter(., Miseq_run==barcode_4_use) ->all.reps_2

all.reps_1 %>% 
  pull(tot) -> reads.per.sample_1

all.reps_2 %>% 
  pull(tot) -> reads.per.sample_2

names(reads.per.sample_1) <- all.reps_1 %>% pull(sample)  
names(reads.per.sample_2) <- all.reps_2 %>% pull(sample)  


### Fit Normal Distribution

fit_1 <- fitdist(reads.per.sample_1, "gamma", lower=c(0,0), start=list(scale=1,shape=1))
fit_2 <- fitdist(reads.per.sample_2, "gamma", lower=c(0,0), start=list(scale=1,shape=1))

all.reps_1 %>%  
  mutate(prob = pgamma(tot, shape = fit_1$estimate[[2]], scale = fit_1$estimate[[1]], lower.tail = TRUE,
                       log.p = FALSE)) -> all.reps_1

all.reps_2 %>%  
  mutate(prob = pgamma(tot, shape = fit_2$estimate[[2]], scale = fit_2$estimate[[1]], lower.tail = TRUE,
                       log.p = FALSE)) -> all.reps_2

### Remove Outlier Samples
outliers_1 <- all.reps_1 %>% 
  filter(prob < low_dist_probability_cutoff  | tot < minimum_read_cutoff | prob > high_dist_probability_cutoff) # changed to 0.05 to save the two samples
outliers_1 %>% View()

outliers_2 <- all.reps_2 %>% 
  filter(prob < low_dist_probability_cutoff  | tot < minimum_read_cutoff | prob > high_dist_probability_cutoff) # changed to 0.05 to save the two samples
outliers_2 %>%  view()

outliers <- rbind(outliers_1,outliers_2)
outliers %>% 
  unite(unique_sample,c("Miseq_run","sample"), remove=F) -> outliers

ASV.nested %>% 
  mutate (Step1.cleaned.tibble.edited = purrr::map(Step1.cleaned.tibble, 
                                           function(.x){
                                             .x %>%
                                               mutate(., barcode_1_col= str_detect(seq_number,barcode_2_use)) %>% 
                                               mutate(., barcode_2_col = str_detect(seq_number,barcode_4_use)) %>% 
                                               mutate(Miseq_run = if_else(barcode_2_col=="FALSE",barcode_2_use,barcode_4_use)) %>% 
                                               unite(unique_sample,c("Miseq_run","sample"), remove=F)
                                           }) ) -> ASV.nested

ASV.nested %>% 
  mutate(Step.2.low.read.depth = purrr::map (Step1.cleaned.tibble.edited, ~ filter(.,!unique_sample %in% outliers$unique_sample) %>% ungroup)) -> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step.2.low.read.depth, ~ how.many(ASVtable = .,round = "2.Low.Read.Depth"))) %>% 
  left_join(ASV.summary) %>% 
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary 

ASV.summary$Summary

#---

#---

# Cleaning Process 3: **Site Occupancy Modeling** ----------------------------------------------------------------------------------------------------------------

##STAN Model

#Metadata test run issue: 2 sequencing runs used different names for the same site so metadata had repeated data. 
#Needed to make a metadata file that was distinct(), otherwise each sample gets duplicated in the spread
metadata %>%  dplyr::select(-Seq_number) %>%  distinct() -> metadata

#Determine which SOM level
  
  ####Set up data with different hierarchical levels and reduce down to unique patterns.
  #Format 2 level hierarcjocal patterns
ASV.nested %>% 
    dplyr::select(Step.2.low.read.depth) %>% 
    unnest(Step.2.low.read.depth) %>% # This stage is a long table
    ungroup() %>%
    mutate(nReads = 1) %>% 
    left_join(metadata, by= c("sample"="New_name")) %>% 
    dplyr::select(seq_number, Site, Bio_rep, Tech_rep, nReads) %>% 
    group_by(Site) %>% 
    mutate(., K = n_distinct(Tech_rep)) -> data_som_int
  
data_som_int %>% 
    pivot_wider(values_from = nReads, names_from = Tech_rep, names_prefix = "Tech") %>% 
    ungroup() %>% 
    mutate(ndetections = rowSums(dplyr::select(., starts_with("Tech")),na.rm = TRUE)) %>% 
    dplyr::select(., !starts_with("Tech"))  %>% 
    pivot_wider(values_from = ndetections, names_from = Bio_rep, names_prefix = "Bio_rep") -> data_tech_summarized
  
  
  #Now need to try to reduce the number of repeated tests. 
  #Many species have identical patterns of presence (e.g. pattern of technical replicates a species was detected)
data_tech_summarized %>% 
    unite(repetition.level,starts_with("Bio_rep"), sep = ".", remove = FALSE) -> Pattern.of.presence_tech
  
  #Select unique patterns
Pattern.of.presence_tech %>% 
    group_by(repetition.level) %>% 
    slice(1) -> unique_patterns
  
  #Format for STAN Model
  # K is number of technical replicates for each Biological replicate
  # ndetections is the number of times a species was detected out of K trials
  
unique_patterns %>% 
    pivot_longer(., cols=starts_with("Bio_rep"), names_to="Bio_rep",values_to="ndetections") %>% 
    filter(., !is.na(ndetections)) %>% 
    mutate(Bio_rep =  str_replace_all(Bio_rep, "Bio_rep","")) -> unique_data
  
  #Create unique identifier for species; for use in hierarchical modeling
  #STAN model requires species names to be numbers.
SS_species <- unite(data = unique_data,
                      col = SS_species,
                      c("seq_number")
  ) %>% pull(SS_species)
  unique_data$Species <- match(SS_species, unique(SS_species)) #index for unique site-species combinations
  
  #Create unique identifier for combinations of site-biological_rep-species; for use in hierarchical modeling
SDS <- unite(data = unique_data,
               col = SDS,
               c("Site", "Bio_rep", "seq_number")
  ) %>% pull(SDS)
unique_data$SiteDateSpecies <- match(SDS, unique(SDS)) #index for unique site-date-species combinations
  
  #Create unique identifier for combinations of site-species; for use in hierarchical modeling
SS <- unite(data = unique_data,
              col = SS,
              c("Site", "seq_number")
  ) %>% pull(SS)
unique_data$SiteSpecies <- match(SS, unique(SS)) #index for unique site-species combinations
  
  #Fix names to match STAN model format
unique_data %>%
    mutate(., N= ndetections) %>%
    dplyr::select(-ndetections)-> unique_data
  
unique_data$N <- as.numeric(unique_data$N)
  ##Stan Model
sink(here("analysis","decontam","Stan_SOM_hierarchical.stan"))
  cat(
    "data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of samples (nrow)
    int<lower=1> Species[S];    // index of species, each of which will have a different value for p11 and p10
    int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
    int<lower=1> L[S];   // index of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> K[S];   // number of replicates per site (ncol)
    int<lower=0> N[S]; // number of detections among these replicates
    int z[S];   // integer flag to help estimate psi parameter
}
parameters{/////////////////////////////////////////////////////////////////////
real<lower=0,upper=1> psi[Nloc];  //commonness parameter
real<lower=0,upper=1> p11[Nspecies]; //true positive detection rate
real<lower=0,upper=1> p10[Nspecies]; //false positive detection rate
}
transformed parameters{/////////////////////////////////////////////////////////////////////
}
model{/////////////////////////////////////////////////////////////////////
real p[S];

for (i in 1:S){
z[i] ~ bernoulli(psi[L[i]]);
p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
N[i] ~ binomial(K[i], p[i]);
}; 

//priors
psi ~ beta(2,2); 
p11 ~ beta(2,2); 
p10 ~ beta(1,10);
}
generated quantities{
real<lower=0,upper=1> Occupancy_prob[S];    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015

for (i in 1:S){
Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
/ ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
+ (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
);
}
}

",
fill=TRUE)
sink()
  
  ##Run Stan model
  #note this will take a while the first time you run a particular model, because it needs to compile from C++
  
myHierarchicalModel <- stan(file = here("analysis","decontam","Stan_SOM_hierarchical.stan"), 
                              data = list(
                                S = nrow(unique_data),
                                Species = unique_data$Species,
                                Nspecies = length(unique(unique_data$Species)),
                                L = unique_data$SiteSpecies,
                                Nloc = length(unique(unique_data$SiteDateSpecies)),
                                K = unique_data$K,
                                N = unique_data$N,
                                z = ifelse(unique_data$N > 0, 1, 0)
                              ), 
                              chains = som_chains,   #number of chains
                              iter = som_iterations   #number of iterations per chain
)
  
  
  
myHierarchicalStanResults <- tidy(tibble(as.data.frame(myHierarchicalModel)))
  
  #See histogram of occupancy probabilities for unique patterns of presence
myHierarchicalStanResults %>% 
    filter(grepl("Occupancy_prob", column)) %>% 
    ggplot() + geom_histogram( aes(mean)) +
    ggtitle(label = "Histogram of Mean Occupancy Probabilities", subtitle = "Unique Patterns of Presence") +
    ylab("Count") + xlab("Mean Occupancy Probability")-> plot_3
  
  ggsave(plot=plot_3,here("analysis","decontam","Output_plots","Histogram_of_Mean_Occupancy_Probabilities_unique_patterns_of_presence.png"), device = "png", width = 12, height = 8, units = "in")
  
  #Now Re-attach SOM to unique patterns
  myHierarchicalStanResults %>% 
    filter(grepl("Occupancy_prob", column)) %>%
    separate(column, into=c("column","SiteDateSpecies"), sep="([\\[\\]])") %>% 
    pivot_wider(., names_from=column, values_from = c(mean, sd))-> myHierarchicalStanResults_wide
  myHierarchicalStanResults_wide$SiteDateSpecies <- as.numeric(myHierarchicalStanResults_wide$SiteDateSpecies)
  
  unique_data %>% 
    left_join(myHierarchicalStanResults_wide) %>% 
    group_by(seq_number,Site) %>% 
    summarise(., max_Occupancy_prob = max(mean_Occupancy_prob)) -> unique_data_SOM
  
  unique_patterns %>% 
    left_join(unique_data_SOM, by=c("seq_number","Site")) %>% 
    dplyr::select(-seq_number,-Site)-> unique_patterns_SOM
  
  Pattern.of.presence_tech %>% 
    left_join(unique_patterns_SOM, by=c("K","repetition.level")) %>% 
    left_join(hash.key)-> SOM_results
  
  #See histogram of occupancy probabilities for unique patterns of presence
  SOM_results %>% 
    ggplot() + geom_histogram( aes(max_Occupancy_prob)) +
    ggtitle(label = "Histogram of Mean Occupancy Probabilities", subtitle = "All ASVs") +
    ylab("Count") + xlab("Mean Occupancy Probability")-> plot_4
  
  ggsave(plot=plot_4,here("analysis","decontam","Output_plots","Histogram_of_Mean_Occupancy_Probabilities_all_ASVs.png"), device = "png", width = 12, height = 8, units = "in")
  
  
  unique_patterns_SOM %>% 
    ungroup() %>% 
    mutate (ndetections = rowSums(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
    filter(., ndetections ==1) %>% 
    summarise(max(max_Occupancy_prob)) -> som_min_threshold_1
  
  unique_patterns_SOM %>% 
    ungroup() %>% 
    mutate(maxdetections = rowMaxs(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
    filter(., maxdetections == 1) %>% 
    summarise(max(max_Occupancy_prob)) -> som_min_threshold_2
  
  unique_patterns_SOM %>% 
    ungroup() %>% 
    mutate(maxdetections = rowMaxs(dplyr::select(., starts_with("Bio_rep")),na.rm = TRUE)) %>% 
    filter(., maxdetections == 3) %>% 
    summarise(min(max_Occupancy_prob)) -> som_min_threshold_3
  
  
  list(som_min_threshold_1,som_min_threshold_2,som_min_threshold_3) -> som_min_thresh
  names(som_min_thresh) <- c("som_min_threshold_1","som_min_threshold_2","som_min_threshold_3")
  
  #Summarize
  SOM_results %>% 
    mutate(., to_keep = case_when(max_Occupancy_prob < som_min_threshold_1$`max(max_Occupancy_prob)` ~ "Toss",
                                  max_Occupancy_prob > som_min_threshold_2$`max(max_Occupancy_prob)` ~ "Keep",
                                  TRUE ~ "Borderline")) %>% 
    count(to_keep) %>% 
    mutate(., total=sum(n)) %>% 
    mutate(per=paste0(round(100*n/total,2),'%')) %>% 
    dplyr::select(-total)

#Save Data
SOM_results %>% 
  dplyr::select(-K, -repetition.level)-> SOM_results_to_keep_all

saveRDS(SOM_results_to_keep_all,here("analysis","decontam","occupancy_results_all.RDS"))

SOM_results_to_keep_all %>% 
  write.csv(here("analysis","decontam",paste0("Occupancy_all_",Sys.Date(),".csv")), row.names = F)

SOM_results_to_keep_all %>% 
  filter(., max_Occupancy_prob >= som_min_thresh[[som_filter_choice]][[1]]) -> SOM_results_to_keep_filtered

saveRDS(SOM_results_to_keep_filtered,here("analysis","decontam","occupancy_results_filtered.RDS"))
SOM_results_to_keep_filtered <- readRDS(here("analysis","decontam","occupancy_results_filtered.RDS"))

SOM_results_to_keep_filtered %>% 
  write.csv(here("analysis","decontam",paste0("Occupancy_all_",Sys.Date(),".csv")), row.names = F)


ASV.nested %>% 
  mutate (Step2.tibble.edited = purrr::map(Step.2.low.read.depth, 
                                           function(.x){
                                             .x %>%
                                               left_join(hash.key) %>% 
                                               dplyr::select(sample, seq_number, nReads, sum.taxonomy)
                                           }) ) -> ASV.nested
ASV.nested %>% 
  mutate(Step3.tibble = purrr::map (Step2.tibble.edited, ~ filter(.,sum.taxonomy %in% SOM_results_to_keep_filtered$sum.taxonomy) %>% ungroup))-> ASV.nested

ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step3.tibble, ~ how.many(ASVtable = .,round ="3.Occupancy"))) %>% 
  left_join(ASV.summary) %>%
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary

ASV.summary$Summary

#---

#---

# Cleaning Process 4: **Dissimilarity between PCR (biological) replicates** ----------------------------------------------------------------------------------------------------------------

ASV.nested %>% 
  dplyr::select(Miseq_run,Step3.tibble) %>% 
  unnest(Step3.tibble) %>% 
  as.data.frame() %>% 
  ungroup() %>% 
  left_join(metadata, by=c("sample"="New_name"))-> cleaned.tibble.post_occ


ASV.nested %>% 
  dplyr::select(Miseq_run,Step2.tibble.edited) %>% 
  unnest(Step2.tibble.edited) %>% 
  as.data.frame() %>% 
  ungroup() %>% 
  left_join(metadata, by=c("sample"="New_name"))-> cleaned.tibble.pre_occ


## How many samples, how many ASVs
cleaned.tibble.post_occ %>% 
  dplyr::summarise(n_distinct(sample),
                   n_distinct(seq_number))

cleaned.tibble.pre_occ %>% 
  dplyr::summarise(n_distinct(sample),
                   n_distinct(seq_number))


## eDNA Index
cleaned.tibble.post_occ %>%
  dplyr::group_by(sample, Miseq_run) %>%
  mutate (Tot = sum(nReads),
          Normalized.reads = nReads / Tot) -> cleaned.tibble.post_occ #transforms raw number of reads to eDNA index


cleaned.tibble.post_occ %>%
  filter(., !str_detect(Miseq_run,"elas")) %>% 
  ungroup() %>% 
  dplyr::select( -Normalized.reads,-Tot,-sample,-Technical_rep) %>% distinct() %>% 
  group_by(sum.taxonomy,Site, Bio_rep,Tech_rep) %>% 
  dplyr::summarise(nReads=sum(nReads)) %>% 
  unite("Replicate",Site,Bio_rep,Tech_rep) %>% ungroup() %>% dplyr::select(Replicate) %>% unique()-> sample_check 


cleaned.tibble.post_occ %>%
  ungroup() %>% 
  dplyr::select( -Normalized.reads,-Tot,-sample,-Technical_rep) %>% distinct() %>% 
  filter(., !str_detect(Miseq_run,"elas")) -> cleaned.tibble.post_occ_v2

cleaned.tibble.post_occ_v2 %>% 
  unite("Replicate",Site,Bio_rep,Tech_rep) %>% 
  group_by(Replicate) %>% 
  dplyr::summarise(tot_reads=sum(nReads)) -> tot_reads

cleaned.tibble.post_occ %>%
  ungroup() %>% 
  dplyr::select( -Normalized.reads,-Tot,-sample,-Technical_rep) %>% distinct() -> cleaned.tibble.post_occ_v3

saveRDS(cleaned.tibble.post_occ_v3, file =here("analysis","decontam","surf_zone_data.RDS"))


cleaned.tibble.post_occ_v2 %>% 
  group_by(sum.taxonomy,Site, Bio_rep,Tech_rep) %>% 
  dplyr::summarise(nReads=sum(nReads)) %>% 
  unite("Replicate",Site,Bio_rep,Tech_rep, remove=F) %>% 
  left_join(tot_reads) %>% 
  mutate(., prop=nReads/tot_reads) %>% 
  dplyr::select(-Replicate,-nReads,-tot_reads) %>% 
  pivot_wider(names_from="Tech_rep", values_from="prop", values_fill=0) %>%
  pivot_longer(cols=`1`:`3`,names_to="Tech_rep", values_to="prop") %>% 
  unite("Replicate",Site,Bio_rep,Tech_rep, remove=F) %>% ungroup() %>% 
  filter(., Replicate %in% sample_check$Replicate) %>% 
  left_join(tot_reads) %>% 
  ggplot(aes(x=tot_reads, y=prop,color=Bio_rep, shape=Tech_rep)) +geom_point() +facet_wrap(sum.taxonomy~Site)

cleaned.tibble.post_occ_v2 %>% 
  group_by(sum.taxonomy,Site, Bio_rep,Tech_rep) %>% 
  dplyr::summarise(nReads=sum(nReads)) %>% 
  unite("Replicate",Site,Bio_rep,Tech_rep, remove=F) %>% 
  left_join(tot_reads) %>% 
  mutate(., prop=nReads/tot_reads) %>% 
  dplyr::select(-Replicate,-nReads,-tot_reads) %>% 
  pivot_wider(names_from="Tech_rep", values_from="prop", values_fill=0) %>%
  pivot_longer(cols=`1`:`3`,names_to="Tech_rep", values_to="prop") %>% 
  unite("Replicate",Site,Bio_rep,Tech_rep, remove=F) %>% ungroup() %>% 
  filter(., Replicate %in% sample_check$Replicate) %>% 
  left_join(tot_reads) -> holder_1

holder_1 %>% 
  group_by(sum.taxonomy,Site, Bio_rep) %>% 
  dplyr::summarise(mean_prop = mean(prop),
         max_prop=max(prop),
         min_prop=min(prop)) ->  avg_data

holder_1 %>% 
  left_join(avg_data) %>% 
  mutate(zero_inf = case_when(min_prop > 0 & max_prop > 0 ~"Three reps above zero",
                              min_prop == 0 & max_prop > 0 ~"One rep is zero",
                              min_prop == 0 & max_prop == 0 ~"All zero")) %>% 
  ggplot(aes(x=tot_reads, y=prop, color=zero_inf, shape=Bio_rep)) +geom_point() +facet_wrap(sum.taxonomy~.,scales = "free") -> plot_prop

 ggsave(plot=plot_prop,here("analysis","decontam","species_prop_plot.png"), device = "png", width = 30, height = 30, units = "in")
  
 
 holder_1 %>% 
   left_join(avg_data) %>% 
   separate(sum.taxonomy, into = c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species"), remove =F, sep=";") %>% 
   dplyr::mutate(name = ifelse(is.na(Species)| Species == "", Genus, Species)) %>% 
   dplyr::mutate(name = ifelse(is.na(Genus) | Genus == "", Family, name)) %>%
   dplyr::mutate(name = ifelse(is.na(Family) | Family == "", Order, name)) %>%
   filter(., !is.na(name)) %>% 
   mutate(zero_inf = case_when(min_prop > 0 & max_prop > 0 ~"Three reps above zero",
                               min_prop == 0 & max_prop > 0 ~"One rep is zero",
                               min_prop == 0 & max_prop == 0 ~"All zero")) %>% 
   mutate(read_depth = case_when(tot_reads > 100000 ~">100K",
                                 tot_reads >= 50000  ~"50K-100K",
                               TRUE ~"<50K")) %>% 
   filter(., prop >0) %>% 
   ggplot(aes(x=zero_inf, y=prop)) +
   geom_violin()  -> plot_prop2
 
 ggsave(plot=plot_prop2,here("analysis","decontam","species_prop_plot2.png"), device = "png", width = 12, height = 8, units = "in", limitsize = FALSE)
 
 
cleaned.tibble.pre_occ %>%
  dplyr::group_by(sample, Miseq_run) %>%
  mutate (Tot = sum(nReads),
          Normalized.reads = nReads / Tot) -> cleaned.tibble.pre_occ #transforms raw number of reads to eDNA index

#Calculate Vegan Distances
tibble_to_vegdist_all (cleaned.tibble.post_occ) -> all.distances.full.post
tibble_to_vegdist_all (cleaned.tibble.pre_occ) -> all.distances.full.pre


as.tibble(subset(melt(as.matrix(all.distances.full.post)))) -> all.distances.melted.post
as.tibble(subset(melt(as.matrix(all.distances.full.pre)))) -> all.distances.melted.pre


  # Now, create a three variables for all distances, they could be PCR replicates, BIOL replicates, or from the same site
all.distances.melted.pre %>%
    separate(Var1, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
    left_join(metadata) %>% 
    mutate(., Site1=Site,
           Bio_rep1=Bio_rep,
           Tech_rep1=Tech_rep,
           Miseq_run1=Miseq_run) %>% 
    dplyr::select(Var1,Miseq_run1,Site1,Bio_rep1,Tech_rep1,Var2,value) %>% 
    separate(Var2, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
    left_join(metadata) %>% 
    mutate(., Site2=Site,
           Bio_rep2=Bio_rep,
           Tech_rep2=Tech_rep,
           Miseq_run2=Miseq_run) %>% 
    dplyr::select(Var1,Var2,Site1,Miseq_run1,Bio_rep1,Tech_rep1,Miseq_run2,Site2,Bio_rep2,Tech_rep2,value) %>%
    unite( Site1, Bio_rep1, col= "station1", remove=F) %>% 
    unite( Site2, Bio_rep2, col= "station2", remove=F) %>% 
    unite( Miseq_run1, station1, col= "typers1", remove=F) %>% 
    unite( Miseq_run2, station2, col= "typers2", remove=F) %>% 
    unite( Miseq_run1, Site1, col= "type_site1", remove=F) %>% 
    unite( Miseq_run2, Site2, col= "type_site2", remove=F) %>% 
    ungroup() %>% 
    mutate(Distance.type = case_when(typers1 == typers2 ~ "PCR.replicates",
                                     type_site1 == type_site2 ~ "Biol.replicates",
                                     Miseq_run1 == Miseq_run2 ~ "Same Barcode Different Site", 
                                     TRUE ~ "Different Barcode")) %>%
    dplyr::select(Sample1 = Var1, Sample2 = Var2 , value , Distance.type) %>%
    filter (Sample1 != Sample2) -> all.distances.to.plot.pre
  
  
  
  # Checking all went well
sapply(all.distances.to.plot.pre, function(x) summary(is.na(x)))
  
all.distances.to.plot.pre$Distance.type <- all.distances.to.plot.pre$Distance.type  %>% fct_relevel( "PCR.replicates", "Biol.replicates", "Same Barcode Different Site", "Different Barcode")
all.distances.to.plot.pre %>% 
  filter(., Distance.type != "Different Barcode") %>% distinct() %>% 
  filter(., !str_detect(Sample1,"elas_")) %>% 
  filter(., !str_detect(Sample2,"elas_"))-> diff_remove

diff_remove %>% 
  filter(., value > .5 ) %>% 
  filter(., Distance.type == "PCR.replicates")  %>% View()

plot_5 <- ggplot (diff_remove , aes (fill = Distance.type, x = value,after_stat(density))) +
    geom_histogram(stat = 'bin', alpha = 0.9, binwidth = 0.05) + #xlim(0, 1) +
    facet_wrap( ~ Distance.type) +
    labs (x = "Pairwise Dissimilarity", y = "Density" ,
          fill = "Groups", title = "eDNA Pairwise Dissimilarity Between Samples", subtitle = "Pre Occupancy") +theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                  panel.background = element_blank())
ggsave(plot=plot_5,here("miu_eDNA_Pairwise_Dissimilarity_Between_Samples_pre_occupancy.png"), device = "png", width = 12, height = 8, units = "in")
  
  
all.distances.melted.post %>%
    separate(Var1, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
    left_join(metadata) %>% 
    mutate(., Site1=Site,
           Bio_rep1=Bio_rep,
           Tech_rep1=Tech_rep,
           Miseq_run1=Miseq_run) %>% 
    dplyr::select(Var1,Miseq_run1,Site1,Bio_rep1,Tech_rep1,Var2,value) %>% 
    separate(Var2, into = c("Miseq_run", "New_name"), sep=":", remove=F) %>% 
    left_join(metadata) %>% 
    mutate(., Site2=Site,
           Bio_rep2=Bio_rep,
           Tech_rep2=Tech_rep,
           Miseq_run2=Miseq_run) %>% 
    dplyr::select(Var1,Var2,Site1,Miseq_run1,Bio_rep1,Tech_rep1,Miseq_run2,Site2,Bio_rep2,Tech_rep2,value) %>%
    unite( Site1, Bio_rep1, col= "station1", remove=F) %>% 
    unite( Site2, Bio_rep2, col= "station2", remove=F) %>% 
    unite( Miseq_run1, station1, col= "typers1", remove=F) %>% 
    unite( Miseq_run2, station2, col= "typers2", remove=F) %>% 
    unite( Miseq_run1, Site1, col= "type_site1", remove=F) %>% 
    unite( Miseq_run2, Site2, col= "type_site2", remove=F) %>% 
    ungroup() %>% 
    mutate(Distance.type = case_when(typers1 == typers2 ~ "PCR.replicates",
                                     type_site1 == type_site2 ~ "Biol.replicates",
                                     Miseq_run1 == Miseq_run2 ~ "Same Barcode Different Site", 
                                     TRUE ~ "Different Barcode")) %>%
    dplyr::select(Sample1 = Var1, Sample2 = Var2 , value , Distance.type) %>%
    filter (Sample1 != Sample2) -> all.distances.to.plot.post
  
  
  # Checking all went well
sapply(all.distances.to.plot.post, function(x) summary(is.na(x)))
  
all.distances.to.plot.post$Distance.type <- all.distances.to.plot.post$Distance.type  %>% fct_relevel( "PCR.replicates", "Biol.replicates", "Same Barcode Different Site", "Different Barcode")
  
plot_6 <- ggplot (all.distances.to.plot.post , aes (fill = Distance.type, x = value,after_stat(density))) +
    geom_histogram(stat = 'bin', alpha = 0.9, binwidth = 0.05) + #xlim(0, 1) +
    facet_wrap( ~ Distance.type) +
    labs (x = "Pairwise Dissimilarity", y = "Density" ,
          fill = "Groups", title = "eDNA Pairwise Dissimilarity Between Samples", subtitle = "Post Site Occupancy") +theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                        panel.background = element_blank())
  
ggsave(plot=plot_6,here("analysis","decontam","Output_plots","eDNA_Pairwise_Dissimilarity_Between_Samples_post_occupancy.png"), device = "png", width = 12, height = 8, units = "in")
  
# Instead of chosing based on the pairwise distances, we will use distance to centroid

  #Need to Remove Single Sampels
cleaned.tibble.post_occ  %>% 
    group_by(Miseq_run,Site,Bio_rep) %>% 
    dplyr::summarise(cases = n_distinct(Tech_rep)) %>% 
    filter(., cases == 1) %>% 
    unite(col="location_id", Miseq_run,Site,Bio_rep, remove = FALSE, sep=":") %>% ungroup() %>% 
    dplyr::select(location_id) -> singles
  #Currently this method does not filter out these single tech reps. There are only 3
  
  #Nest Data
cleaned.tibble.post_occ  %>% 
    unite(col="location_id", Miseq_run,Site,Bio_rep, remove = FALSE, sep=":") %>% 
    filter(., !location_id %in% singles$location_id) %>% 
    group_by(location_id) %>% nest() -> nested.cleaning
  
nested.cleaning %>% 
    mutate(matrix = map(data, tibble_to_vegdist_bio_rep)) -> nested.cleaning
  
  # Convert
nested.cleaning %>% mutate(ncomparisons = map(matrix, length)) -> nested.cleaning
  
  #Calculate distances to centroid
nested.cleaning <- nested.cleaning %>% mutate(distances = map2(matrix, location_id, dist_to_centroid))
  nested.cleaning %>% 
    separate(location_id, c("Miseq_run","Site","Bio_rep"), sep=":") %>% 
    unnest_longer(distances) %>% 
    dplyr::select(Miseq_run,Site,Bio_rep, distances_id,distances) -> all_distances.groups #unnest data


#Calculate normal distribution of distances to centroid
normparams.step5 <- MASS::fitdistr(all_distances.groups$distances, "normal")$estimate

#Calculate Probability
probs.step5 <- pnorm(all_distances.groups$distances, normparams.step5[1], normparams.step5[2])

#Determine Outliers
outliers.step5 <- which(probs.step5>max_dist_probability_cutoff)

#Remove Outliers
discard.step5 <-names(all_distances.groups$distances[outliers.step5])

to_write_discarded.step5 <- tibble(distances_id = discard.step5,
                                   distance = all_distances.groups$distances[outliers.step5])

to_write_discarded.step5 <- to_write_discarded.step5 %>% left_join(tibble(distances_id = discard.step5,
                                                                          probs = probs.step5[outliers.step5]))
write_csv(to_write_discarded.step5 ,here("analysis","decontam","step5.discarded_samples.csv"))

## Plot Final Replication Levels
  cleaned.tibble.post_occ %>% 
    unite(col="distances_id", Miseq_run,Site,Bio_rep,Tech_rep, remove = FALSE, sep=":") ->cleaned.tibble.post_occ_for_plotting
  all_distances.groups %>%
    dplyr::select(-Miseq_run,-Site,-Bio_rep) -> all_distances.groups_for_plotting
  
  cleaned.tibble.post_occ_for_plotting %>%
    filter(., !distances_id %in% discard.step5) %>% 
    dplyr::group_by(Miseq_run,Site,Bio_rep) %>% 
    dplyr::summarise(cases = n_distinct(distances_id)) %>% 
    ggplot()+
    geom_raster(aes(x= Site, y = Bio_rep, fill = cases))+
    geom_text(aes(x= Site, y = Bio_rep, label = cases),color="white") +
    ggtitle("Sample Replication Level") +facet_wrap(~Miseq_run) +
    theme(axis.text.x = element_text(angle = 90))-> plot_7
  
  ggsave(plot=plot_7,here("analysis","decontam","Output_plots","Final_sample_replication_level.png"), device = "png", width = 20, height = 8, units = "in")


##Remove Samples
ASV.nested %>% 
    mutate (Step3.tibble_edited = purrr::map(Step3.tibble, 
                                             function(.x){
                                               .x %>%
                                                 mutate(., barcode_1_col= str_detect(seq_number,barcode_2_use)) %>% 
                                                 mutate(., barcode_2_col = str_detect(seq_number,barcode_4_use)) %>% 
                                                 mutate(Miseq_run = if_else(barcode_2_col=="FALSE",barcode_2_use,barcode_4_use)) %>% 
                                                 left_join(metadata, by=c("sample"="New_name")) %>% 
                                                 unite(col="distances_id", Miseq_run,Site,Bio_rep,Tech_rep, remove = FALSE, sep=":")
                                             }) ) -> ASV.nested
  
  
#Filter Sample
ASV.nested %>% 
    mutate(Step4.tibble = map (Step3.tibble_edited,  ~ filter(.,! distances_id %in% to_write_discarded.step5$distances_id))) -> ASV.nested


#Visualize Results of Clearence Process 4
ASV.nested %>% 
  transmute(Miseq_run, Summary.1 = purrr::map(Step4.tibble, ~ how.many(ASVtable = .,round ="4.Dissimilarity"))) %>% 
  left_join(ASV.summary) %>%
  mutate(Summary   = map2(Summary, Summary.1, bind_rows)) %>%
  dplyr::select(-Summary.1) -> ASV.summary


ASV.summary$Summary

##Save Data
saveRDS(ASV.nested,file=here("analysis","decontam","ASV.nested_final.RDS"))
saveRDS(ASV.summary,file=here("analysis","decontam","ASV.summary_final.RDS"))


#---

#---

#Generate Final Outputs ----------------------------------------------------------------------------------------------------------------

##Code for merging ASV tables

###Identify Unique Species for Sum.taxonomy merging

hash.key %>% 
  distinct(.,sum.taxonomy) -> hashes_unique

hashes_unique$number <- row.names(hashes_unique)
hashes_unique$number <- paste0("taxon_",hashes_unique$number)
row.names(hashes_unique)<-hashes_unique$number

hash.key %>% 
  left_join(hashes_unique, by="sum.taxonomy") -> hash.key.updated

#---

#---

###Pre Occupancy Sum by Taxonomy, All PCR Tech Reps Separate Samples
  
hash.key.updated$number %>% unique() -> total_taxa
  
  ASV.nested$Step3.tibble[[1]] %>% 
    mutate(miseq = ASV.nested$Miseq_run[[1]]) %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    unite(miseq,sample, col="Sample") %>% 
    left_join(hash.key.updated) %>% 
    dplyr::group_by(number,Sample) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    spread(., Sample, nReads) %>% #convert to wide data format
    replace(is.na(.), 0) -> barcode_1_pre_df
  barcode_1_pre_df$number %>%  unique() -> barcode_1_pre_df_taxa
  
  ASV.nested$Step3.tibble[[2]] %>% 
    mutate(miseq = ASV.nested$Miseq_run[[2]]) %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    unite(miseq,sample, col="Sample") %>% 
    left_join(hash.key.updated) %>% 
    dplyr::group_by(number,Sample) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    spread(., Sample, nReads) %>% #convert to wide data format
    replace(is.na(.), 0) -> barcode_2_pre_df
  barcode_2_pre_df$number %>%  unique() -> barcode_2_pre_df_taxa
  
  total_kept_taxa <- (append(barcode_1_pre_df_taxa,barcode_2_pre_df_taxa)) %>% unique()
  
  if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_pre_df_taxa))) {
    barcode_1_pre_df %>% ungroup() %>% 
      arrange(number) -> barcode_1_pre_df
  } else {
    barcode_1_pre_df %>% ungroup() %>% 
      add_row(number=setdiff(total_kept_taxa,barcode_1_pre_df_taxa)) %>% 
      arrange(number) %>% 
      replace(is.na(.), 0) -> barcode_1_pre_df
  }
  
  if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_pre_df_taxa))) {
    barcode_2_pre_df %>% ungroup() %>% 
      arrange(number) -> barcode_2_pre_df
  } else {
    barcode_2_pre_df %>% ungroup() %>% 
      add_row(number=setdiff(total_kept_taxa,barcode_2_pre_df_taxa)) %>% 
      arrange(number) %>% 
      replace(is.na(.), 0) -> barcode_2_pre_df
  }
  
  barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
  row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
  barcode_1_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_pre_df
  
  barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)
  row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
  barcode_2_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_pre_df
  
  ####First, we want to create proportions by dividing by the rowsums:
  ####we could do this with sweep() or mutate_all() or other ways, but using vegan:
  
  barcode_1_prop <- decostand(barcode_1_pre_df, method = "total", MARGIN = 2)
  barcode_2_prop <- decostand(barcode_2_pre_df, method = "total", MARGIN = 2)
  
  ####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
  ####We do this by scaling everything to the max observed in each row. 
  
  ####to do this WITHIN a dataset, we could just do (again, using vegan):
  barcode_1_prop_index <- decostand(barcode_1_prop, method = "max", MARGIN = 1)
  barcode_2_prop_index <- decostand(barcode_2_prop, method = "max", MARGIN = 1)
  
  ####This gives us an index between 0 and 1 for each species in each dataset.  
  
  ####But if we want to combine datasets, this second step has to happen in the combined dataset, so it all gets scaled to 0-1.  
  ####easy enough:
  
  combined_index <- decostand(cbind(barcode_1_prop,barcode_2_prop), method = "max", MARGIN = 1)
  ####How both datasets are combined, on a common, comparable scale.
  
  ### Output Read Count Data
  pre_results_reads = cbind(barcode_1_prop,barcode_2_prop)
  
  hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]
  
  pre_results_reads$number <- rownames(pre_results_reads)
  
  pre_results_reads %>% 
    left_join(hash.key.updated.2, by="number") %>% 
    dplyr::select(-number,-seq_number) -> pre_results_reads
  
  saveRDS(pre_results_reads,file=here("analysis","decontam","Output_R","pre_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.RDS"))
  write_csv(pre_results_reads ,here("analysis","decontam","Output_csv","pre_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.csv"))
  
  ### Output eDNA Index Data
  
  combined_index$number <- rownames(combined_index)
  
  combined_index %>% 
    left_join(hash.key.updated.2, by="number") %>% 
    dplyr::select(-number,-seq_number) -> combined_index
  
  saveRDS(combined_index,file=here("analysis","decontam","Output_R","pre_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.RDS"))
  write_csv(combined_index ,here("analysis","decontam","Output_csv","pre_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.csv"))
  
  #---
  
  #---
  
  ###Pre Occupancy Sum by Taxonomy, All PCR Tech Reps Separate Samples
  
  hash.key.updated$number %>% unique() -> total_taxa
  
  ASV.nested$Step4.tibble[[1]] %>% 
    mutate(miseq = ASV.nested$Miseq_run[[1]]) %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    unite(miseq,sample, col="Sample") %>% 
    left_join(hash.key.updated) %>% 
    dplyr::group_by(number,Sample) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    spread(., Sample, nReads) %>% #convert to wide data format
    replace(is.na(.), 0) -> barcode_1_post_occ_df
  barcode_1_post_occ_df$number %>%  unique() -> barcode_1_post_occ_df_taxa
  
  
  ASV.nested$Step4.tibble[[2]] %>% 
    mutate(miseq = ASV.nested$Miseq_run[[2]]) %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    unite(miseq,sample, col="Sample") %>% 
    left_join(hash.key.updated) %>% 
    dplyr::group_by(number,Sample) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    spread(., Sample, nReads) %>% #convert to wide data format
    replace(is.na(.), 0) -> barcode_2_post_occ_df
  barcode_2_post_occ_df$number %>%  unique() -> barcode_2_post_occ_df_taxa
  
  total_kept_taxa <- append(barcode_1_post_occ_df_taxa,barcode_2_post_occ_df_taxa) %>% unique()
  
  if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_post_occ_df_taxa))) {
    barcode_1_post_occ_df %>% ungroup() %>% 
      arrange(number) -> barcode_1_post_occ_df
  } else {
    barcode_1_post_occ_df %>% ungroup() %>% 
      add_row(number=setdiff(total_kept_taxa,barcode_1_post_occ_df_taxa)) %>%
      arrange(number) %>% 
      replace(is.na(.), 0) -> barcode_1_post_occ_df
  }
  
  if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_post_occ_df_taxa))) {
    barcode_2_post_occ_df %>% ungroup() %>% 
      arrange(number) -> barcode_2_post_occ_df
  } else {
    barcode_2_post_occ_df %>% ungroup() %>% 
      add_row(number=setdiff(total_kept_taxa,barcode_2_post_occ_df_taxa)) %>% 
      arrange(number) %>% 
      replace(is.na(.), 0) -> barcode_2_post_occ_df
  }
  
  
  barcode_1_post_occ_df <- as.data.frame(barcode_1_post_occ_df)
  row.names(barcode_1_post_occ_df) <- barcode_1_post_occ_df$number
  barcode_1_post_occ_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_post_occ_df
  
  barcode_2_post_occ_df <- as.data.frame(barcode_2_post_occ_df)
  row.names(barcode_2_post_occ_df) <- barcode_2_post_occ_df$number
  barcode_2_post_occ_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_post_occ_df
  
  ####first, we want to create proportions by dividing by the rowsums:
  ####we could do this with sweep() or mutate_all() or other ways, but using vegan:
  
  barcode_1_post_occ_df_prop <- decostand(barcode_1_post_occ_df, method = "total", MARGIN = 2)
  barcode_2_post_occ_df_prop <- decostand(barcode_2_post_occ_df, method = "total", MARGIN = 2)
  
  ####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
  ####We do this by scaling everything to the max observed in each row. 
  
  ####to do this WITHIN a dataset, we could just do (again, using vegan):
  barcode_1_post_occ_df_index <- decostand(barcode_1_post_occ_df_prop, method = "max", MARGIN = 1)
  barcode_2_post_occ_df_index <- decostand(barcode_2_post_occ_df_prop, method = "max", MARGIN = 1)
  
  ####this gives us an index between 0 and 1 for each species in each dataset.  
  
  ####But if we want to combine datasets, this second step has to happen in the combined dataset, so it all gets scaled to 0-1.  
  ####easy enough:
  
  combined_index <- decostand(cbind(barcode_1_post_occ_df_index,barcode_2_post_occ_df_index), method = "max", MARGIN = 1)
  ####now both datasets are combined, on a common, comparable scale.
  
  ### Output Post Occupancy Read Data
  results_post_occ_reads = cbind(barcode_1_post_occ_df,barcode_2_post_occ_df)
  
  hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]
  
  results_post_occ_reads$number <- rownames(results_post_occ_reads)
  
  results_post_occ_reads %>% 
    left_join(hash.key.updated.2, by="number") %>% 
    dplyr::select(-number,-seq_number) -> results_post_occ_reads
  
  saveRDS(results_post_occ_reads,file=here("analysis","decontam","Output_R","post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.RDS"))
  write_csv(results_post_occ_reads ,here("analysis","decontam","Output_csv","post_occupancy_results_sum.taxonomy_tech_reps_separate_read_counts.csv"))
  
  ### Output Post Occupancy eDNA Index Data
  combined_index$number <- rownames(combined_index)
  
  combined_index %>% 
    left_join(hash.key.updated.2, by="number") %>% 
    dplyr::select(-number,-seq_number) -> combined_index
  
  saveRDS(combined_index,file=here("analysis","decontam","Output_R","post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.RDS"))
  write_csv(combined_index ,here("analysis","decontam","Output_csv","post_occupancy_results_sum.taxonomy_tech_reps_separate_eDNA_index.csv"))
  
  #Code for Merging Tech Reps
  
  ###Pre Occupancy Sum by Taxonomy, Biological Replicates Separate
  
  ASV.nested$Step3.tibble[[1]] %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    mutate(Sample = sample) %>% 
    left_join(hash.key.updated) %>% 
    dplyr::group_by(number,Sample) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    spread(., Sample, nReads) %>% #convert to wide data format
    replace(is.na(.), 0)  -> barcode_1_pre_df
  barcode_1_pre_df$number %>%  unique() -> barcode_1_pre_df_taxa
  
  metadata %>% 
    filter(., !(New_name %in% controls)) %>% 
    filter(.,!(New_name %in% colnames(barcode_1_pre_df))) %>% 
    filter(., !New_name %in% Gold_mpawatch) %>% 
    pull(New_name) %>%  unique()-> columns2add
  
  barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
  
  barcode_1_pre_df %>% 
    tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
    replace(is.na(.), 0) %>% 
    dplyr::select(sort(tidyselect::peek_vars()))-> barcode_1_pre_df
  
  row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
  barcode_1_pre_df %>% ungroup() -> barcode_1_pre_df
  
  ASV.nested$Step3.tibble[[2]] %>% 
    mutate(Sample = sample) %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    left_join(hash.key.updated) %>% 
    dplyr::group_by(number,Sample) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    spread(., Sample, nReads) %>% #convert to wide data format
    replace(is.na(.), 0)  -> barcode_2_pre_df
  barcode_2_pre_df$number %>%  unique() -> barcode_2_pre_df_taxa
  
  metadata %>% 
    filter(., !(New_name %in% controls)) %>% 
    filter(.,!(New_name %in% colnames(barcode_2_pre_df))) %>% 
    filter(., !New_name %in% Gold_mpawatch) %>% 
    pull(New_name) %>%  unique()-> columns2add
  
  barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)
  
  barcode_2_pre_df %>% 
    tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
    replace(is.na(.), 0) %>% 
    dplyr::select(sort(tidyselect::peek_vars())) -> barcode_2_pre_df
  
  row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
  barcode_2_pre_df %>% ungroup() -> barcode_2_pre_df
  
  total_kept_taxa <- (append(barcode_1_pre_df_taxa,barcode_2_pre_df_taxa)) %>% unique()
  
  if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_pre_df_taxa))) {
    barcode_1_pre_df %>% ungroup() %>% 
      arrange(number) -> barcode_1_pre_df
  } else {
    barcode_1_pre_df %>% ungroup() %>% 
      add_row(number=setdiff(total_kept_taxa,barcode_1_pre_df_taxa)) %>% 
      arrange(number) %>% 
      replace(is.na(.), 0) -> barcode_1_pre_df
  }
  
  if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_pre_df_taxa))) {
    barcode_2_pre_df %>% ungroup() %>% 
      arrange(number) -> barcode_2_pre_df
  } else {
    barcode_2_pre_df %>% ungroup() %>% 
      add_row(number=setdiff(total_kept_taxa,barcode_2_pre_df_taxa)) %>% 
      arrange(number) %>% 
      replace(is.na(.), 0) -> barcode_2_pre_df
  }
  
  barcode_1_pre_df <- as.data.frame(barcode_1_pre_df)
  row.names(barcode_1_pre_df) <- barcode_1_pre_df$number
  barcode_1_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_pre_df
  
  barcode_2_pre_df <- as.data.frame(barcode_2_pre_df)
  row.names(barcode_2_pre_df) <- barcode_2_pre_df$number
  barcode_2_pre_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_pre_df
  
  barcode_1_pre_df+barcode_2_pre_df -> combined_df
  
  ####First, we want to create proportions by dividing by the rowsums:
  ####we could do this with sweep() or mutate_all() or other ways, but using vegan:
  
  combined_df_prop <- decostand(combined_df, method = "total", MARGIN = 2)
  
  ####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
  ####We do this by scaling everything to the max observed in each row. 
  
  ####To do this WITHIN a dataset, we could just do (again, using vegan):
  ####eDNA Index between 0-1 made by straight adding reads across each decontaminated datatable
  combined_df_index <- decostand(combined_df_prop, method = "max", MARGIN = 1)
  
  #Output read counts
  hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]
  
  combined_df$number <- rownames(combined_df)
  
  combined_df %>% 
    left_join(hash.key.updated.2, by="number") %>% 
    dplyr::select(-number,-seq_number) -> combined_df
  
  saveRDS(combined_df,file=here("analysis","decontam","Output_R","pre_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.RDS"))
  write_csv(combined_df,here("analysis","decontam","Output_csv","pre_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.csv"))
  
  ### Output eDNA Index
  combined_df_index$number <- rownames(combined_df_index)
  
  combined_df_index %>% 
    left_join(hash.key.updated.2, by="number") %>% 
    dplyr::select(-number,-seq_number) -> combined_df_index
  
  saveRDS(combined_df_index,file=here("analysis","decontam","Output_R","pre_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.RDS"))
  write_csv(combined_df_index,here("analysis","decontam","Output_csv","pre_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.csv"))
  
  
  ###Post Occupancy Merge, Bio reps sep
  
  ASV.nested$Step4.tibble[[1]] %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    mutate(Sample = sample) %>% 
    left_join(hash.key.updated) %>% 
    dplyr::group_by(number,Sample) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    spread(., Sample, nReads) %>% #convert to wide data format
    replace(is.na(.), 0)  -> barcode_1_post_df
  barcode_1_post_df$number %>%  unique() -> barcode_1_post_df_taxa
  
  metadata %>% 
    filter(., !(New_name %in% controls)) %>% 
    filter(.,!(New_name %in% colnames(barcode_1_post_df))) %>% 
    filter(., !New_name %in% Gold_mpawatch) %>% 
    pull(New_name) %>%  unique()-> columns2add
  
  barcode_1_post_df <- as.data.frame(barcode_1_post_df)
  
  barcode_1_post_df %>% 
    tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
    replace(is.na(.), 0) %>% 
    dplyr::select(sort(tidyselect::peek_vars()))-> barcode_1_post_df
  
  row.names(barcode_1_post_df) <- barcode_1_post_df$number
  barcode_1_post_df %>% ungroup() -> barcode_1_post_df
  
  ASV.nested$Step4.tibble[[2]] %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    mutate(Sample = sample) %>% 
    left_join(hash.key.updated) %>% 
    dplyr::group_by(number,Sample) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    spread(., Sample, nReads) %>% #convert to wide data format
    replace(is.na(.), 0)  -> barcode_2_post_df
  barcode_2_post_df$number %>%  unique() -> barcode_2_post_df_taxa
  
  metadata %>% 
    filter(., !(New_name %in% controls)) %>% 
    filter(.,!(New_name %in% colnames(barcode_2_post_df))) %>% 
    filter(., !New_name %in% Gold_mpawatch) %>% 
    pull(New_name) %>%  unique()-> columns2add
  
  barcode_2_post_df <- as.data.frame(barcode_2_post_df)
  
  barcode_2_post_df %>% 
    tibble::add_column(!!!set_names(as.list(rep(NA, length(columns2add))),nm=columns2add)) %>% 
    replace(is.na(.), 0) %>% 
    dplyr::select(sort(tidyselect::peek_vars())) -> barcode_2_post_df
  
  row.names(barcode_2_post_df) <- barcode_2_post_df$number
  barcode_2_post_df %>% ungroup() -> barcode_2_post_df
  
  total_kept_taxa <- (append(barcode_1_post_df_taxa,barcode_2_post_df_taxa)) %>% unique()
  
  if (rlang::is_empty(setdiff(total_kept_taxa,barcode_1_post_df_taxa))) {
    barcode_1_post_df %>% ungroup() %>% 
      arrange(number) -> barcode_1_post_df
  } else {
    barcode_1_post_df %>% ungroup() %>% 
      add_row(number=setdiff(total_kept_taxa,barcode_1_post_df_taxa)) %>% 
      arrange(number) %>% 
      replace(is.na(.), 0) -> barcode_1_post_df
  }
  
  if (rlang::is_empty(setdiff(total_kept_taxa,barcode_2_post_df_taxa))) {
    barcode_2_post_df %>% ungroup() %>% 
      arrange(number) -> barcode_2_post_df
  } else {
    barcode_2_post_df %>% ungroup() %>% 
      add_row(number=setdiff(total_kept_taxa,barcode_2_post_df_taxa)) %>% 
      arrange(number) %>% 
      replace(is.na(.), 0) -> barcode_2_post_df
  }
  
  barcode_1_post_df <- as.data.frame(barcode_1_post_df)
  row.names(barcode_1_post_df) <- barcode_1_post_df$number
  barcode_1_post_df %>% ungroup() %>% dplyr::select(-number) -> barcode_1_post_df
  
  barcode_2_post_df <- as.data.frame(barcode_2_post_df)
  row.names(barcode_2_post_df) <- barcode_2_post_df$number
  barcode_2_post_df %>% ungroup() %>% dplyr::select(-number) -> barcode_2_post_df
  
  barcode_1_post_df+barcode_2_post_df -> combined_df
  
  ####First, we want to create proportions by dividing by the rowsums:
  ####we could do this with sweep() or mutate_all() or other ways, but using vegan:
  
  combined_df_prop <- decostand(combined_df, method = "total", MARGIN = 2)
  
  ####Second, we want to ask how the proprortion for each species has changed across columns (samples). 
  ####We do this by scaling everything to the max observed in each row. 
  
  ####To do this WITHIN a dataset, we could just do (again, using vegan):
  ####eDNA Index between 0-1 made by straight adding reads across each decontaminated datatable
  combined_df_index <- decostand(combined_df_prop, method = "max", MARGIN = 1)
  
  #Output read counts
  hash.key.updated.2 <- hash.key.updated[!duplicated(hash.key.updated$number), ]
  
  combined_df$number <- rownames(combined_df)
  
  combined_df %>% 
    left_join(hash.key.updated.2, by="number") %>% 
    dplyr::select(-number,-seq_number) -> combined_df
  
  saveRDS(combined_df,file=here("analysis","decontam","Output_R","post_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.RDS"))
  write_csv(combined_df,here("analysis","decontam","Output_csv","post_occupancy_results_sum.taxonomy_tech_reps_summed_read_counts.csv"))
  
  ### Output eDNA Index
  combined_df_index$number <- rownames(combined_df_index)
  
  combined_df_index %>% 
    left_join(hash.key.updated.2, by="number") %>% 
    dplyr::select(-number,-seq_number) -> combined_df_index
  
  saveRDS(combined_df_index,file=here("analysis","decontam","Output_R","post_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.RDS"))
  write_csv(combined_df_index,here("analysis","decontam","Output_csv","post_occupancy_results_sum.taxonomy_tech_reps_summed_eDNA_index.csv"))
  
  #---
  
  #---
  
  #Code for Merging Sites
  
  ###Pre Occupancy Merge by Site
  ASV.nested %>% 
    dplyr::select(Step3.tibble,Miseq_run) %>% 
    unnest(Step3.tibble) %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    left_join(hash.key.updated) %>% 
    left_join(metadata, by=c("sample"="New_name")) %>% 
    dplyr::group_by(sum.taxonomy,Site) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    pivot_wider(., names_from=Site, values_from=nReads, values_fill=0)  -> pre_wide_reads
  
  saveRDS(pre_wide_reads,file=here("analysis","decontam","Output_R","pre_occ_site_sum.taxonomy_reads_summed.RDS"))
  write_csv(pre_wide_reads ,here("analysis","decontam","Output_csv","pre_occ_site_sum.taxonomy_reads_summed.csv"))
  
  ASV.nested %>% 
    dplyr::select(Step3.tibble,Miseq_run) %>% 
    unnest(Step3.tibble) %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    left_join(hash.key.updated) %>% 
    left_join(metadata, by=c("sample"="New_name")) %>% 
    ungroup() %>% 
    group_by(Site, sum.taxonomy) %>% 
    summarise(meanreads = mean(nReads)) %>% 
    dplyr::group_by(Site) %>%
    mutate (Tot = sum(meanreads),
            Row.sums = meanreads / Tot) %>% 
    dplyr::group_by (sum.taxonomy) %>%
    mutate (Colmax = max (Row.sums),
            Normalized.reads = Row.sums / Colmax) %>% 
    ungroup() %>% 
    dplyr::select(Site,Normalized.reads, sum.taxonomy) %>% 
    pivot_wider(., names_from=Site, values_from=Normalized.reads, values_fill=0)  -> pre_wide
  
  saveRDS(pre_wide,file=here("analysis","decontam","Output_R","pre_occ_site_averaged_sum.taxonomy_e_index.RDS"))
  write_csv(pre_wide, here("analysis","decontam","Output_csv","pre_occ_site_averaged_sum.taxonomy_e_index.csv"))
  
  ###Post Occupancy Merge by Site
  ASV.nested %>% 
    dplyr::select(Step4.tibble) %>% 
    unnest(Step4.tibble) %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    ungroup() %>% 
    dplyr::group_by(sum.taxonomy,Site) %>%
    dplyr::summarise(nReads=sum(nReads)) %>% 
    pivot_wider(., names_from=Site, values_from=nReads, values_fill=0) -> post_wide_reads
  
  saveRDS(post_wide_reads,file=here("analysis","decontam","Output_R","post_occ_site_sum.taxonomy_reads_summed.RDS"))
  write_csv(post_wide_reads ,here("analysis","decontam","Output_csv","post_occ_site_sum.taxonomy_reads_summed.csv"))
  
  ASV.nested %>% 
    dplyr::select(Step4.tibble) %>% 
    unnest(Step4.tibble) %>% 
    filter(., !sample %in% Gold_mpawatch) %>% 
    ungroup() %>% 
    group_by(Site, sum.taxonomy) %>% 
    summarise(meanreads = mean(nReads)) %>% 
    dplyr::group_by(Site) %>%
    mutate (Tot = sum(meanreads),
            Row.sums = meanreads / Tot) %>% 
    dplyr::group_by (sum.taxonomy) %>%
    mutate (Colmax = max(Row.sums),
            Normalized.reads = Row.sums / Colmax) %>% 
    dplyr::select(Site,Normalized.reads, sum.taxonomy) %>% 
    pivot_wider(., names_from=Site, values_from=Normalized.reads, values_fill=0)   -> post_wide
  
  
saveRDS(post_wide,file=here("analysis","decontam","Output_R","post_occ_site_averaged_sum.taxonomy_e_index.RDS"))
write_csv(post_wide,here("analysis","decontam","Output_csv","post_occ_site_averaged_sum.taxonomy_e_index.csv"))

post_wide$sum.taxonomy %>% sort()
