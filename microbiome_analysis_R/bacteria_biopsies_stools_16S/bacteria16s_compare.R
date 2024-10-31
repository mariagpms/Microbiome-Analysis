##################################################################################################################################################
##The aim of this code is to help to compare the diversity of the bacteria in both stools and biopsies within the different type of IBD patients##
##################################################################################################################################################

#First of all, the working directory must be changed to where all the data is, so it is not necessary to set the complete path when loading other objects
setwd("path_to_where_your_files_are")

#Then, a phyloseq object needs to be created, in doing so the next webpage has been very useful
#############################################################################################
#    Title: Phyloseq Introduction and Import: Phyloseq and Microbiome analysis in R         #
#    Author: Siobhon Egan during an internship program at Pawsey Supercomputing Centre      #
#    Date: 2017/2018                                                                        #
#    Availability: https://cryptick-lab.github.io/NGS-Analysis/_site/R-PhyloseqIntro.html   #
#############################################################################################

#Next, the csv files that have the information related to the otus, taxonomy and the metadata of the experiment will be read 
library(readr) #Library needed to read the csv files

#Dataset of otus, the rows must be the otus and the columns must be the samples
uparse_otus <- read_csv("otu_mat.csv")
#As explained in Siobhon Egab's page, both the uparse_otus dataset and the taxonomy matrix need to be ordered by the same column (the one containing the otu names)
sorting<-order(uparse_otus$`sample-id`)
uparse_otus<-uparse_otus[sorting,]
#The first column, which is the one containing the otus, needs to be deleted because the phyloseq function needs a particular format. Then it will be transformed into a matrix 
#and set both sample and rownames (it is very important that the order of the columns in the otu matrix is the same that the order of the rows in the taxonomy matrix)
otumat <- as.matrix(uparse_otus[,-1])
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
head(otumat)

#Dataset with the information of taxonomy, must be ordered by the otus so that it has the same order as the otu matrix. Then delete the column of otus and convert it into a matrix
taxonomy <- read_delim("tax-mat.csv", delim = ";",escape_double = FALSE, trim_ws = TRUE)
taxonomy<-taxonomy[sorting,]
taxmat <- as.matrix(taxonomy[,-1])
rownames(taxmat) <- rownames(otumat)
#The names of the columns will be set in the way it is needed by phyloseq function
colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
head(taxmat)

#By doing this it is ensured that the datasets are indeed converted into a matrix
class(otumat)
class(taxmat)

#The library that is needed to create the phyloseq object is loaded
library(phyloseq)
#The physeq object with both the otu table and tax table is created
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
TAX
physeq = phyloseq(OTU, TAX)
physeq
sample_names(physeq)

#Dataset with the the meta data information is converted so that the columns are factors
#IBD_type: CONTROL, ACTIVE UC, QUIESCENT UC, ACTIVE CROHN, QUIESCENT CROHN
#Status: IBD, healthy
#Sample_type: Biopsies, Stools
meta_data <- read_delim("metadatamerge.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
sampledata = sample_data(data.frame(meta_data, row.names=sample_names(physeq), stringsAsFactors=FALSE))
head(sampledata)

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
saveRDS(physeq1,"phyloseq_merge.RDS")

#The aim of this section of code is to create a box and whiskers plot with the alpha diversity with several diversity measures, in doing so the next webpage has been of great use
############################################################################################
#    Title: Alpha Diversity and Sequencing Depth                                           #
#    Author: Siobhon Egan during an internship program at Pawsey Supercomputing Centre     #
#    Date: 2017/2018                                                                       #
#    Availability: https://cryptick-lab.github.io/NGS-Analysis/_site/R-AlphaDiversity.html #
############################################################################################
library(ggplot2)
physeq_prune <- prune_species(speciesSums(physeq1) > 0, physeq1)#prune the physeq object
#First, it needs the pruned physeq object, in here it is both stools and biopsies which does not make sense but it was just for the purpose of learning how to do it
pwhisk<-plot_richness(physeq_prune, measures=c("Simpson","Shannon"),x = "IBD_type",color = "IBD_type")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
pwhisk+geom_boxplot() 

pwhisk_sample<-plot_richness(physeq_prune, measures=c("Simpson","Shannon"),x = "Sample_type",color = "Sample_type")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
pwhisk_sample+geom_boxplot()

#Now a plot for biopsies and another for stools will be created
#The physeq object has to be subset to have one for biopsies and another for stools. Then, both objects will be pruned and a box and whiskers plot will be done with Simpson's and Shannon's measure
physeq1_biop<-subset_samples(physeq1,Sample_type=="Biopsies")
physeq1_st<-subset_samples(physeq1,Sample_type=="Stools")
#Warning: when using prune_species and speciesSum a warning will be shown on the console, it is because it will be deprecated as long as the object is created (check the environment) there is no problema in still using it
#Warning message:
#  1: In prune_species(speciesSums(quies_colitis) > 0, quies_colitis) : 'prune_species' will be deprecated in the future.
#Use 'prune_taxa' instead.
#See help("Deprecated") and help("phyloseq-deprecated").
#2: In speciesSums(quies_colitis) : 'speciesSums' will be deprecated in the future.
#Use 'taxa_sums' instead.
#See help("Deprecated") and help("phyloseq-deprecated").
biop_prune<-prune_species(speciesSums(physeq1_biop) > 0, physeq1_biop)
st_prune<-prune_species(speciesSums(physeq1_st)>0,physeq1_st)

#The alpha diversity plots for biopsies will be created
#First, calculate the diversity distances for both Simpson's and Shannon's measures
measures_biop<-estimate_richness(biop_prune,measures = c("Simpson","Shannon"))
measures_biop$names_samples<-rownames(measures_biop)
measures_biop<-measures_biop%>%left_join(relation_names_id,by=join_by("names_samples"=="code"))
#The following library is loaded to use fct_relevel() which will reorder the levels of IBD_type factor
library(forcats)
#Then, the boxplots of alpha diversity for Simpson's and Shannon's measure are created
biop_Simpson<-measures_biop%>%
  mutate(name = fct_relevel(IBD_type,"CONTROL","ACTIVE CROHN","QUIESCENT CROHN","ACTIVE UC","QUIESCENT UC"))%>%
  ggplot(aes(x=name,y=Simpson,fill=name,alpha=0.3))+
  geom_boxplot()+geom_jitter(aes(color=name),width = 0.2,size=2,alpha=1)+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=20,face="bold"),axis.title.x = element_blank())
biop_Simpson

biop_Shannon<-measures_biop%>%
  mutate(name = fct_relevel(IBD_type,"CONTROL","ACTIVE CROHN","QUIESCENT CROHN","ACTIVE UC","QUIESCENT UC"))%>%
  ggplot(aes(x=name,y=Shannon,fill=name,alpha=0.3))+
  geom_boxplot()+geom_jitter(aes(color=name),width = 0.2,size=2,alpha=1)+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=20,face="bold"),axis.title.x = element_blank())
biop_Shannon
#The following library is loaded to modify the theme and title of the plots when they are set together
library(patchwork)
x11()
biop_Simpson+biop_Shannon+plot_annotation(title = "Alpha diversity in Biopsies")&theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))

#Now, it will be done for stools
measures_st<-estimate_richness(st_prune,measures = c("Simpson","Shannon"))
measures_st$names_samples<-rownames(measures_st)
measures_st<-measures_st%>%left_join(relation_names_id,by=join_by("names_samples"=="code"))

st_Simpson<-measures_st%>%
  mutate(name = fct_relevel(IBD_type,"CONTROL","ACTIVE CROHN","QUIESCENT CROHN","ACTIVE UC","QUIESCENT UC"))%>%
  ggplot(aes(x=name,y=Simpson,fill=name,alpha=0.3))+
  geom_boxplot()+geom_jitter(aes(color=name),width = 0.2,size=2,alpha=1)+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=20,face="bold"),axis.title.x = element_blank())
st_Simpson

st_Shannon<-measures_st%>%
  mutate(name = fct_relevel(IBD_type,"CONTROL","ACTIVE CROHN","QUIESCENT CROHN","ACTIVE UC","QUIESCENT UC"))%>%
  ggplot(aes(x=name,y=Shannon,fill=name,alpha=0.3))+
  geom_boxplot()+geom_jitter(aes(color=name),width = 0.2,size=2,alpha=1)+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5,size=20,face="bold"),axis.title.x = element_blank())
st_Shannon

x11()
st_Simpson+st_Shannon+plot_annotation(title = "Alpha diversity in Stool")&theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))

#Now, the Kruskal-Wallis test for both biopsies and stools will be performed
measures_biop<-measures_biop%>%mutate(IBD_type = fct_relevel(IBD_type,"CONTROL","ACTIVE CROHN","QUIESCENT CROHN","ACTIVE UC","QUIESCENT UC"))
kruskal.test(Simpson ~ IBD_type, data = measures_biop)
kruskal.test(Shannon ~ IBD_type, data = measures_biop)

measures_st<-measures_st%>%mutate(IBD_type = fct_relevel(IBD_type,"CONTROL","ACTIVE CROHN","QUIESCENT CROHN","ACTIVE UC","QUIESCENT UC"))
kruskal.test(Simpson ~ IBD_type, data = measures_st)
kruskal.test(Shannon ~ IBD_type, data = measures_st)

#Now, lefse is going to be used to represent some plots that will help to determine which species are augmented or decreased between biopsies and stools
#Install microbiomeMaker and load it
BiocManager::install("microbiomeMarker") 
library(microbiomeMarker)

#All types together (Control, active and quiescent both crohn and uc) evaluating biopsies vs stools
lefse_merge <- run_lefse(
  physeq_prune,
  group = "Sample_type",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_merge #there were 112 microbiome markers identified
x11()
plot_ef_bar(lefse_merge)+ggtitle("Stools-Biopsies (all): Genus")

#Do it within family instead of Genus as done before
lefse_family <- run_lefse(
  physeq_prune,
  group = "Sample_type",
  subgroup = NULL,
  taxa_rank = "Family",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_family
x11()
plot_ef_bar(lefse_family)+ggtitle("Stools-Biopsies (all): Family")

#Now with species despite with our data it is not advised to do so
lefse_species <- run_lefse(
  physeq_prune,
  group = "Sample_type",
  subgroup = NULL,
  taxa_rank = "Species",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_species
x11()
plot_ef_bar(lefse_species)

#The physeq object will be subset to see the differences in the groups

#Control: Stools-Biopsies
controls<-subset_samples(physeq1, status=="healthy")
controls_prune <- prune_species(speciesSums(controls) > 0, controls)
lefse_contols <- run_lefse(
  controls_prune,
  group = "Sample_type",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_contols #26 microbiome markers
p5<-plot_ef_bar(lefse_contols)+ggtitle("Control")
x11()
p5

#Check the differences between healthy and IBD
lefse_status <- run_lefse(
  physeq_prune,
  group = "status",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_status
plot_ef_bar(lefse_status)+xlim(c(0,5))+ggtitle("IBD-Healthy (all)")

#Within the patients with IBD
sick<-subset_samples(physeq1, status=="IBD")
sick_prune <- prune_species(speciesSums(sick) > 0, sick)
lefse_sick <- run_lefse(
  sick_prune,
  group = "Sample_type",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_sick
x11()
plot_ef_bar(lefse_sick)+ggtitle("IBD: Stools-Biopsies")

#Subset it by patients with Crohn's disease
crohn<-subset_samples(sick, IBD_type=="ACTIVE CROHN"|IBD_type=="QUIESCENT CROHN")
crohn_prune<-prune_species(speciesSums(crohn) > 0, crohn)
lefse_crohn<- run_lefse(
  crohn_prune,
  group = "Sample_type",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_crohn
x11()
plot_ef_bar(lefse_crohn)+ggtitle("Crohn: Stools-Biopsies")

#Subset it between active and quiescent crohn to do the lefse analysis between stools and biopsies in within each group
active_crohn<-subset_samples(crohn,IBD_type=="ACTIVE CROHN")
quies_crohn<-subset_samples(crohn,IBD_type=="QUIESCENT CROHN")
prune_accrohn<-prune_species(speciesSums(active_crohn) > 0, active_crohn)
prune_quicrohn<-prune_species(speciesSums(quies_crohn) > 0, quies_crohn)

lefse_accrohn<- run_lefse(
  prune_accrohn,
  group = "Sample_type",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_accrohn
lefse_quicrohn<- run_lefse(
  prune_quicrohn,
  group = "Sample_type",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_quicrohn
p3<-plot_ef_bar(lefse_accrohn)+ggtitle("Active Crohn")
p4<-plot_ef_bar(lefse_quicrohn)+ggtitle("Quiescent Crohn")
x11()
p4+p3

#Do the same with colitis
colitis<-subset_samples(sick, IBD_type=="ACTIVE UC"|IBD_type=="QUIESCENT UC")
colitis_prune<-prune_species(speciesSums(colitis) > 0, colitis)
lefse_colitis<- run_lefse(
  colitis_prune,
  group = "Sample_type",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_colitis
x11()
plot_ef_bar(lefse_colitis)+ggtitle("Ulcerative Colitis: Stools-Biopsies")

active_colitis<-subset_samples(colitis,IBD_type=="ACTIVE UC")
quies_colitis<-subset_samples(colitis,IBD_type=="QUIESCENT UC")
prune_acuc<-prune_species(speciesSums(active_colitis) > 0, active_colitis)
prune_quiuc<-prune_species(speciesSums(quies_colitis) > 0, quies_colitis)
lefse_acuc<- run_lefse(
  prune_acuc,
  group = "Sample_type",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_acuc
lefse_quiuc<- run_lefse(
  prune_quiuc,
  group = "Sample_type",
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
  curv = FALSE # logical, whether perform the wilcoxon test using the Curtis’s approach, defalt FALSE
)
lefse_quiuc
p1<-plot_ef_bar(lefse_acuc)+ggtitle("Active UC")
p2<-plot_ef_bar(lefse_quiuc)+ggtitle("Quiescent UC")
x11()
p1+p2

#Active
x11()
p1+p3

#Quiescent
x11()
p2+p4

#All
x11()
p1+p2+p5+p3+p4
