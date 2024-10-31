#####################################################################################################################
##The aim of this code is to help to compare the diversity of fungi in biopsies within the different type of groups##
#####################################################################################################################

#First of all, the working directory must be changed to where all the data is, so it is not necessary to set the complete path when loading other objects
setwd("path_to_where_your_files_are")

#Then, a phyloseq object needs to be created, in doing so the next webpage has been very useful
#############################################################################################
#    Title: Phyloseq Introduction and Import: Phyloseq and Microbiome analysis in R         #
#    Author: Siobhon Egan during an internship program at Pawsey Supercomputing Centre      #
#    Date: 2017/2018                                                                        #
#    Availability: https://cryptick-lab.github.io/NGS-Analysis/_site/R-PhyloseqIntro.html   #
#############################################################################################

#Next, the csv files, that have the information related to the otus, will be read , taxonomy and the metadata of the experiment
library(readr) #Library needed to read the csv files
library(readxl)#Library needed to read the excel files

#Dataset of otus, the rows must be the otus and the columns must be the samples
uparse_otus <- read_excel("otu_mat.xlsx",col_types = c("text",rep("numeric",25)))
#As explained in Siobhon Egab's page, both the uparse_otus dataset and the taxonomy matrix need to be ordered by the same column (the one containing the otu names)
sorting<-order(uparse_otus$id)
uparse_otus<-uparse_otus[sorting,]
#The first column, which is the one containing the otus, needs to be deleted because the phyloseq function needs a particular format,  then it will be transformed into a matrix 
#and set both sample and rownames (it is very important that the order of the columns in the otu matrix is the same that the order of the rows in the taxonomy matrix)
otumat <- as.matrix(uparse_otus[,-1])
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
head(otumat)

#Dataset with the information of taxonomy, must be ordered by the otus, so that it has the same order that the otu matrix, delete the column of otus and convert it into a matrix
taxonomy <- read_delim("tax_mat.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
taxonomy<-taxonomy[sorting,]
taxmat <- as.matrix(taxonomy[,-1])
rownames(taxmat) <- rownames(otumat)
#The names of the columns will be set in the way it is needed by phyloseq
colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
head(taxmat)

#By doing this it is ensured that the datasets are indeed converted into a matrix
class(otumat)
class(taxmat)

#The library that is needed to create the phyloseq object is loaded
library("phyloseq")
#The physeq object with both the otu table and tax table is created
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
TAX
physeq = phyloseq(OTU, TAX)
physeq
sample_names(physeq)

#Dataset with the the meta data information, it is converted so that the columns are factors
#IBD_type: CONTROL, ACTIVE UC, QUIESCENT UC, ACTIVE CROHN, QUIESCENT CROHN
meta_data <- read_excel("metadatafungi.xlsx")
sampledata = sample_data(data.frame(meta_data, row.names=sample_names(physeq), stringsAsFactors=FALSE))
head(sampledata)
relation_names_id<-relation_names_id%>%left_join(meta_data,by = join_by("sample_id"=="sample-id"))

#The random tree associated to the physeq object is created, one created with qiime2 or any other tool could have been loaded instead of creating a new one
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

#Create final phyloseq object
physeq1 = merge_phyloseq(physeq, sampledata, random_tree)
physeq1

#If you do not have the package BiocManager installed it will be installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#The library microbiome from BiocManager will be installed and then loaded
BiocManager::install("microbiome") 
library(microbiome)


#This is just to save the physeq object, it is not strictly needed
saveRDS(physeq1,"phyloseq_fungi.RDS")

#The aim of this section of code is to create a box and whiskers plot with the alpha diversity with several diversity measure measures, in doing so the next webpage has been of great use
############################################################################################
#    Title: Alpha Diversity and Sequencing Depth                                           #
#    Author: Siobhon Egan during an internship program at Pawsey Supercomputing Centre     #
#    Date: 2017/2018                                                                       #
#    Availability: https://cryptick-lab.github.io/NGS-Analysis/_site/R-AlphaDiversity.html #
############################################################################################
library(ggplot2)
physeq_prune <- prune_species(speciesSums(physeq1) > 0, physeq1)#prune the physeq object

#Create the plots of alpha diversity
#First, the diversity distances for both Shannon's and Simpson's measures are calculated
measures<-estimate_richness(physeq_prune,measures = c("Simpson","Shannon"))
measures$names_samples<-rownames(measures)
measures<-measures%>%left_join(relation_names_id,by=join_by("names_samples"=="code"))

#The following library is loaded to use fct_relevel() which will reorder the levels of IBD_type factor
library(forcats)
#The following library is loaded to modify the theme and title of the plots when they are set together
library(patchwork)

#Then, the plots will be ploted and then joined to have them side by side
Simpson<-measures%>%
  mutate(name = fct_relevel(IBD_type,"CONTROL","ACTIVE CROHN","QUIESCENT CROHN","ACTIVE UC","QUIESCENT UC"))%>%
  ggplot(aes(x=name,y=Simpson,fill=name,alpha=0.3))+
  geom_boxplot()+geom_jitter(aes(color=name),width = 0.2,size=2,alpha=1)+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=20,face="bold"),axis.title.x = element_blank())
Simpson

Shannon<-measures%>%
  mutate(name = fct_relevel(IBD_type,"CONTROL","ACTIVE CROHN","QUIESCENT CROHN","ACTIVE UC","QUIESCENT UC"))%>%
  ggplot(aes(x=name,y=Shannon,fill=name,alpha=0.3))+
  geom_boxplot()+geom_jitter(aes(color=name),width = 0.2,size=2,alpha=1)+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=20,face="bold"),axis.title.x = element_blank())
Shannon

x11()
Simpson+Shannon+plot_annotation(title = "Alpha diversity in Fungi Biopsies")&theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))

#Now, the Kruskal-Wallis test will be performed
measures<-measures%>%mutate(IBD_type = fct_relevel(IBD_type,"CONTROL","ACTIVE CROHN","QUIESCENT CROHN","ACTIVE UC","QUIESCENT UC"))
kruskal.test(Simpson ~ IBD_type, data = measures)
kruskal.test(Shannon ~ IBD_type, data = measures)

#Now, lefse is going to be used to represent some plots that will help to determine which species are augmented or decreased between groups of IBD_type
#Install microbiomeMaker and load it
BiocManager::install("microbiomeMarker") 
library(microbiomeMarker)

#Subset the physeq object into controls and sick
controls<- subset_samples(physeq1, IBD_type=="CONTROL")
nsamples(physeq)
nsamples(controls)

samples_only<-subset_samples(physeq1, IBD_type=="ACTIVE UC"|IBD_type=="QUIESCENT CROHN"|IBD_type=="ACTIVE CROHN"|IBD_type=="QUIESCENT UC")
nsamples(samples_only)

#The objects are pruned 
controls_prune <- prune_species(speciesSums(controls) > 0, controls)
samples_only_prune <- prune_species(speciesSums(samples_only)>0, samples_only)

#Trying to aply lefse between active and quiescent colitis
colitis<- subset_samples(physeq1, IBD_type=="QUIESCENT UC"|IBD_type=="ACTIVE UC")
colitis_prune <- prune_species(speciesSums(colitis) > 0, colitis)

#Lefse can use different normalization methods, all were tried and none identified biomarkers. It is probably due to the fact that there are only 5 patients in each group
lefse_colitis <- run_lefse(
  colitis_prune,
  group = "IBD_type",
  subgroup = NULL,
  taxa_rank = "Genus",
  transform = "identity",
  norm = "CPM",
  kw_cutoff = 0.05, # default =0.05
  lda_cutoff = 2, # default=2
  bootstrap_n = 30, # integer, the number of bootstrap iteration for LDA, default 30
  bootstrap_fraction = 2/3, # numeric, the subsampling fraction value for each bootstrap iteration, default 2/3
  wilcoxon_cutoff = 0.05, # default =0.05
  multigrp_strat = FALSE, # logical, for multiple group tasks, whether the test is performed in a one-against one (more strict) or in a one-against all setting, default FALSE.
  strict = "0", # multiple testing options, 0 for no correction (default), 1 for independent comparisons, 2 for independent comparison.
  sample_min = 2, # integer, minimum number of samples per subclass for performing wilcoxon test, default 10
  only_same_subgrp = FALSE, # logical, whether perform the wilcoxon test only among the subgroups with the same name, default FALSE
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_colitis
#As no biomarkers were identified, the 'plot_ef_bar' plot can not be done

#The same is done with Crohn
crohn<- subset_samples(physeq1, IBD_type=="QUIESCENT CROHN"|IBD_type=="ACTIVE CROHN")
crohn_prune <- prune_species(speciesSums(crohn) > 0, crohn)
#Tried different options of normalization and transformation but no biomarkers were identified either
lefse_crohn <- run_lefse(
  crohn_prune,
  group = "IBD_type",
  subgroup = NULL,
  taxa_rank = "Genus",
  transform = "identity",
  norm = "CPM",
  kw_cutoff = 0.05, # default =0.05
  lda_cutoff = 2, # default=2
  bootstrap_n = 30, # integer, the number of bootstrap iteration for LDA, default 30
  bootstrap_fraction = 2/3, # numeric, the subsampling fraction value for each bootstrap iteration, default 2/3
  wilcoxon_cutoff = 0.05, # default =0.05
  multigrp_strat = FALSE, # logical, for multiple group tasks, whether the test is performed in a one-against one (more strict) or in a one-against all setting, default FALSE.
  strict = "0", # multiple testing options, 0 for no correction (default), 1 for independent comparisons, 2 for independent comparison.
  sample_min = 2, # integer, minimum number of samples per subclass for performing wilcoxon test, default 10
  only_same_subgrp = FALSE, # logical, whether perform the wilcoxon test only among the subgroups with the same name, default FALSE
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_crohn

#Quiescent IBD (Crohn and Colitis) vs Active IBD (Crohn and Colitis) is now tried. Same as before, different options were tried but no biomarkers were identified
library (dplyr)
samples_only@sam_data@.Data[[2]]<-ifelse(samples_only@sam_data@.Data[[2]] %in% c("ACTIVE UC","ACTIVE CROHN"),"active", "quiescent")
samples_only@sam_data
samples_only<- subset_samples(samples_only, IBD_type=="active"|IBD_type=="quiescent")
samples_only_prune <- prune_species(speciesSums(samples_only)>0, samples_only)
samples_only_prune@tax_table<-samples_only_prune@tax_table[!is.na(samples_only_prune@tax_table[,6]),]
lefse_activity <- run_lefse(
  samples_only_prune,
  group = "IBD_type",
  subgroup = NULL,
  taxa_rank = "Genus",
  transform = "log10",
  norm = "CPM",
  kw_cutoff = 0.05, # default =0.05
  lda_cutoff = 2, # default=2
  bootstrap_n = 30, # integer, the number of bootstrap iteration for LDA, default 30
  bootstrap_fraction = 2/3, # numeric, the subsampling fraction value for each bootstrap iteration, default 2/3
  wilcoxon_cutoff = 0.05, # default =0.05
  multigrp_strat = FALSE, # logical, for multiple group tasks, whether the test is performed in a one-against one (more strict) or in a one-against all setting, default FALSE.
  strict = "0", # multiple testing options, 0 for no correction (default), 1 for independent comparisons, 2 for independent comparison.
  sample_min = 2, # integer, minimum number of samples per subclass for performing wilcoxon test, default 10
  only_same_subgrp = FALSE, # logical, whether perform the wilcoxon test only among the subgroups with the same name, default FALSE
  curv = FALSE) # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE

lefse_activity

#No biomarkers were identifed in none of the different tests that were tried. It is most probably due to the fact that there are just 5 patients in each group and they are quite similar (as shown in the pcoa plot)
