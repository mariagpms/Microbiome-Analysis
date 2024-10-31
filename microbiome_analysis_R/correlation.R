#######################################################################################################
#The aim of this code is to calculate the correlations between bacteria and fungi                     #                     
#in biopsies for each group. Then, a heatmap will be ploted                                           #
#######################################################################################################
#First of all, the required libraries will be loaded: to read data, to process it and to plot it
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)

#First of all, the working directory must be changed to where the fungi data is, so it is not necessary to set the complete path when loading the objects
setwd("Path/microbione_analysis_R/fungi_biopsies_18S")

#The otu dataset will be loaded, and changed the column name to the same as the taxonomy
uparse_otus_fungi <- read_excel("otu_mat.xlsx",col_types = c("text",rep("numeric",25)))
colnames(uparse_otus_fungi)[1]<-"Feature ID"
#The taxonomy dataset is loaded
taxonomy_fungi <- read_delim("tax_mat.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
#The metadata dataset is loaded
meta_data_fungi <- read_excel("metadatafungi.xlsx")

#In the otu's dataset, the column with the otu ids will be changed to the Genus that correspond to each otu
#It will be done my joining the otu dataset with the taxonomy one by the column of otus (Feature ID). Then, all the columns except for the one corresponding to Genus, the number 7, will be deleted from the final dataset
otu_tax_fungi<-uparse_otus_fungi%>%left_join(taxonomy_fungi,by="Feature ID")%>%select(-colnames(taxonomy_fungi)[-7])

#With the next 5 lines of code the names of the samples for each IBD_type are obtained, they are the same for bacteria and fungi
control<-meta_data_fungi[meta_data_fungi[,2]=="CONTROL",1]
acuc<-meta_data_fungi[meta_data_fungi[,2]=="ACTIVE UC",1]
quiuc<-meta_data_fungi[meta_data_fungi[,2]=="QUIESCENT UC",1]
accrohn<-meta_data_fungi[meta_data_fungi[,2]=="ACTIVE CROHN",1]
quicrohn<-meta_data_fungi[meta_data_fungi[,2]=="QUIESCENT CROHN",1]

#For the 5 levels of IBD_type the dataset will be summarized by the Genus
#Then, those genus which a value of 0 for all the patients will be deleted. 
#This is done by creating a new column that sums the frequencies of all the patients and then the rows that have 0 for this column will be deleted
fungi_control<-otu_tax_fungi%>%select(Genus,control$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(control$`sample-id`),sum)
fungi_control$sum_all<-rowSums(subset(fungi_control,select = control$`sample-id`))
fungi_control<-fungi_control%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))

fungi_acuc<-otu_tax_fungi%>%select(Genus,acuc$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(acuc$`sample-id`),sum)
fungi_acuc$sum_all<-rowSums(subset(fungi_acuc,select = acuc$`sample-id`))
fungi_acuc<-fungi_acuc%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))

fungi_quiuc<-otu_tax_fungi%>%select(Genus,quiuc$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(quiuc$`sample-id`),sum)
fungi_quiuc$sum_all<-rowSums(subset(fungi_quiuc,select = quiuc$`sample-id`))
fungi_quiuc<-fungi_quiuc%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))

fungi_accrohn<-otu_tax_fungi%>%select(Genus,accrohn$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(accrohn$`sample-id`),sum)
fungi_accrohn$sum_all<-rowSums(subset(fungi_accrohn,select = accrohn$`sample-id`))
fungi_accrohn<-fungi_accrohn%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))

fungi_quicrohn<-otu_tax_fungi%>%select(Genus,quicrohn$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(quicrohn$`sample-id`),sum)
fungi_quicrohn$sum_all<-rowSums(subset(fungi_quicrohn,select = quicrohn$`sample-id`))
fungi_quicrohn<-fungi_quicrohn%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))

#Second of all, the working directory must be changed to where the bacteria data is, so it is not necessary to set the complete path when loading the objects
setwd("Path/microbione_analysis_R/bacteria_biopsies_stools_16S")

#Both the metadata file as well as the taxonomy matrix are loaded
meta_data_bac <- read_delim("metadatamerge.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
tax_mat_bac <- read_delim("tax-mat.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)

#Once the otu matrix is in a csv it is loaded
otu_mat_bac <- read_csv("otu_mat.csv")
colnames(otu_mat_bac)[1]<-"Feature ID" #the first column does not have a name, it is set to be the same as in the taxonomy matrix, it refers to the column containing the otu identifiers

#The otu matrix is split into a biopsies otu matrix and a stools otu matrix. The first column is needed as it contains the otus and the rest of the columns that are being kept correspond to the samples of biopsies or stools
#Keep in mind that this could have been done differently, for example extracting the sample names for biopsies/stools from the metadata matrix and then giving that as well as the first column name ("Feature ID") to get the same output
otu_mat_biop<-otu_mat_bac[,1:26]
otu_biop_tax<-otu_mat_biop%>%left_join(tax_mat_bac,by="Feature ID")%>%select(-colnames(tax_mat_bac)[-7])
colnames(otu_biop_tax)<-colnames(otu_tax_fungi)#Although the sample names are the same, the ones in bacteria have a suffix so we replace then with the others

#For the 5 levels of IBD_type the dataset will be summarized by the Genus, the procedure will be the same as before but then a column with the means will be calculated
#Then, the top 20 frequencies will be kept, to have the top 20 bacteria 
bac_control<-otu_biop_tax%>%select(Genus,control$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(control$`sample-id`),sum)
bac_control$sum_all<-rowSums(subset(bac_control,select = control$`sample-id`))
bac_control<-bac_control%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))
bac_control$mean_top<-rowMeans(subset(bac_control,select = control$`sample-id`))
bac_control<-bac_control%>%arrange(desc(mean_top))%>%head(20)%>%select(-mean_top)

bac_acuc<-otu_biop_tax%>%select(Genus,acuc$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(acuc$`sample-id`),sum)
bac_acuc$sum_all<-rowSums(subset(bac_acuc,select = acuc$`sample-id`))
bac_acuc<-bac_acuc%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))
bac_acuc$mean_top<-rowMeans(subset(bac_acuc,select = acuc$`sample-id`))
bac_acuc<-bac_acuc%>%arrange(desc(mean_top))%>%head(20)%>%select(-mean_top)

bac_quiuc<-otu_biop_tax%>%select(Genus,quiuc$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(quiuc$`sample-id`),sum)
bac_quiuc$sum_all<-rowSums(subset(bac_quiuc,select = quiuc$`sample-id`))
bac_quiuc<-bac_quiuc%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))
bac_quiuc$mean_top<-rowMeans(subset(bac_quiuc,select = quiuc$`sample-id`))
bac_quiuc<-bac_quiuc%>%arrange(desc(mean_top))%>%head(20)%>%select(-mean_top)

bac_accrohn<-otu_biop_tax%>%select(Genus,accrohn$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(accrohn$`sample-id`),sum)
bac_accrohn$sum_all<-rowSums(subset(bac_accrohn,select = accrohn$`sample-id`))
bac_accrohn<-bac_accrohn%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))
bac_accrohn$mean_top<-rowMeans(subset(bac_accrohn,select = accrohn$`sample-id`))
bac_accrohn<-bac_accrohn%>%arrange(desc(mean_top))%>%head(20)%>%select(-mean_top)

bac_quicrohn<-otu_biop_tax%>%select(Genus,quicrohn$`sample-id`)%>%group_by(Genus)%>%summarise_at(vars(quicrohn$`sample-id`),sum)
bac_quicrohn$sum_all<-rowSums(subset(bac_quicrohn,select = quicrohn$`sample-id`))
bac_quicrohn<-bac_quicrohn%>%filter(sum_all!=0)%>%select(-sum_all)%>%filter(!is.na(Genus))
bac_quicrohn$mean_top<-rowMeans(subset(bac_quicrohn,select = quicrohn$`sample-id`))
bac_quicrohn<-bac_quicrohn%>%arrange(desc(mean_top))%>%head(20)%>%select(-mean_top)


#Union of the dataframes for bacteria and fungi
bac_fungi_control<-union(bac_control,fungi_control)
control_genus<-bac_fungi_control%>%select(Genus)
bac_fungi_control<-bac_fungi_control%>%select(-Genus)
bac_fungi_control<-t(bac_fungi_control)
colnames(bac_fungi_control)<-control_genus$Genus

bac_fungi_acuc<-union(bac_acuc,fungi_acuc)
acuc_genus<-bac_fungi_acuc%>%select(Genus)
bac_fungi_acuc<-bac_fungi_acuc%>%select(-Genus)
bac_fungi_acuc<-t(bac_fungi_acuc)
colnames(bac_fungi_acuc)<-acuc_genus$Genus

bac_fungi_quiuc<-union(bac_quiuc,fungi_quiuc)
quiuc_genus<-bac_fungi_quiuc%>%select(Genus)
bac_fungi_quiuc<-bac_fungi_quiuc%>%select(-Genus)
bac_fungi_quiuc<-t(bac_fungi_quiuc)
colnames(bac_fungi_quiuc)<-quiuc_genus$Genus

bac_fungi_accrohn<-union(bac_accrohn,fungi_accrohn)
accrohn_genus<-bac_fungi_accrohn%>%select(Genus)
bac_fungi_accrohn<-bac_fungi_accrohn%>%select(-Genus)
bac_fungi_accrohn<-t(bac_fungi_accrohn)
colnames(bac_fungi_accrohn)<-accrohn_genus$Genus

bac_fungi_quicrohn<-union(bac_quicrohn,fungi_quicrohn)
quicrohn_genus<-bac_fungi_quicrohn%>%select(Genus)
bac_fungi_quicrohn<-bac_fungi_quicrohn%>%select(-Genus)
bac_fungi_quicrohn<-t(bac_fungi_quicrohn)
colnames(bac_fungi_quicrohn)<-quicrohn_genus$Genus

#The following library will be loaded to calculate the Spearman's correlation and then perform a significance test
library(rstatix)
#Now, the Spearman's correlations between the frequencies of bacteria and fungi will be calculated and a significance test performed
#Then, a heatmap will be ploted. The correlations that are significant at a 5% level will be marked on the plot
#This will be done for the 5 groups of IBD_type
#Control
cor_control_test<-as.data.frame(bac_fungi_control)%>%cor_test(method = "spearman")%>%
    filter(!((var1%in%bac_control$Genus)&(var2%in%bac_control$Genus)))%>%
    filter(!((var1%in%fungi_control$Genus)&(var2%in%fungi_control$Genus)))%>%
    mutate(var1_factor = factor(var1, levels = unique(c(var1, var2))),
         var2_factor = factor(var2, levels = unique(c(var1, var2))))%>%
    filter(as.numeric(var1_factor) < as.numeric(var2_factor))

cor_control_test%>%ggplot(aes(var1, var2, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
            limit = c(-1,1), space = "Lab", name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
  ggtitle("Heatmap - Control")+xlab("Bacteria")+ylab("Fungi")+
  coord_fixed()

cor_control_sig<-mutate(cor_control_test, sig = ifelse(p <0.05, "Sig.", "Non Sig."))
ggplot()+geom_tile(data=cor_control_sig,aes(var1, var2, fill = cor),color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
        limit = c(-1,1), space = "Lab", name="Correlation") +
  geom_tile(data = filter(cor_control_sig, sig == "Sig."), aes(var1, var2),
        linewidth = 1, colour = "black", fill = "transpa rent")+
  geom_text(data = cor_control_sig, aes(var1, var2, label = round(cor, 2),
                fontface = ifelse(sig == "Sig.", "bold", "plain")))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
    ggtitle("Heatmap - Control")+xlab("Bacteria")+ylab("Fungi")+
    coord_fixed()

#Active UC
cor_acuc_test<-as.data.frame(bac_fungi_acuc)%>%cor_test(method = "spearman")%>%
    filter(!((var1%in%bac_acuc$Genus)&(var2%in%bac_acuc$Genus)))%>%
    filter(!((var1%in%fungi_acuc$Genus)&(var2%in%fungi_acuc$Genus)))%>%
    mutate(var1_factor = factor(var1, levels = unique(c(var1, var2))),
         var2_factor = factor(var2, levels = unique(c(var1, var2))))%>%
    filter(as.numeric(var1_factor) < as.numeric(var2_factor))

cor_acuc_test%>%ggplot(aes(var1, var2, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
            limit = c(-1,1), space = "Lab", name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
  ggtitle("Heatmap - Active UC")+xlab("Bacteria")+ylab("Fungi")+
  coord_fixed()

cor_acuc_sig<-mutate(cor_acuc_test, sig = ifelse(p <0.05, "Sig.", "Non Sig."))
ggplot()+geom_tile(data=cor_acuc_sig,aes(var1, var2, fill = cor),color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
        limit = c(-1,1), space = "Lab", name="Correlation") +
  geom_tile(data = filter(cor_acuc_sig, sig == "Sig."), aes(var1, var2),
        linewidth = 1, colour = "black", fill = "transpa rent")+
  geom_text(data = cor_acuc_sig, aes(var1, var2, label = round(cor, 2),
                fontface = ifelse(sig == "Sig.", "bold", "plain")))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
    ggtitle("Heatmap - Active UC")+xlab("Bacteria")+ylab("Fungi")+
    coord_fixed()

#Quiescent UC
cor_quiuc_test<-as.data.frame(bac_fungi_quiuc)%>%cor_test(method = "spearman")%>%
    filter(!((var1%in%bac_quiuc$Genus)&(var2%in%bac_quiuc$Genus)))%>%
    filter(!((var1%in%fungi_quiuc$Genus)&(var2%in%fungi_quiuc$Genus)))%>%
    mutate(var1_factor = factor(var1, levels = unique(c(var1, var2))),
         var2_factor = factor(var2, levels = unique(c(var1, var2))))%>%
    filter(as.numeric(var1_factor) < as.numeric(var2_factor))

cor_quiuc_test%>%ggplot(aes(var1, var2, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
            limit = c(-1,1), space = "Lab", name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
  ggtitle("Heatmap - Quiescent UC")+xlab("Bacteria")+ylab("Fungi")+
  coord_fixed()

cor_quiuc_sig<-mutate(cor_quiuc_test, sig = ifelse(p <0.05, "Sig.", "Non Sig."))
ggplot()+geom_tile(data=cor_quiuc_sig,aes(var1, var2, fill = cor),color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
        limit = c(-1,1), space = "Lab", name="Correlation") +
  geom_tile(data = filter(cor_quiuc_sig, sig == "Sig."), aes(var1, var2),
        linewidth = 1, colour = "black", fill = "transpa rent")+
  geom_text(data = cor_quiuc_sig, aes(var1, var2, label = round(cor, 2),
                fontface = ifelse(sig == "Sig.", "bold", "plain")))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
    ggtitle("Heatmap - Quiescent UC")+xlab("Bacteria")+ylab("Fungi")+
    coord_fixed()

#Active Crohn
cor_accrohn_test<-as.data.frame(bac_fungi_accrohn)%>%cor_test(method = "spearman")%>%
    filter(!((var1%in%bac_accrohn$Genus)&(var2%in%bac_accrohn$Genus)))%>%
    filter(!((var1%in%fungi_accrohn$Genus)&(var2%in%fungi_accrohn$Genus)))%>%
    mutate(var1_factor = factor(var1, levels = unique(c(var1, var2))),
         var2_factor = factor(var2, levels = unique(c(var1, var2))))%>%
    filter(as.numeric(var1_factor) < as.numeric(var2_factor))

cor_accrohn_test%>%ggplot(aes(var1, var2, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
            limit = c(-1,1), space = "Lab", name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
  ggtitle("Heatmap - Active Crohn")+xlab("Bacteria")+ylab("Fungi")+
  coord_fixed()

cor_accrohn_sig<-mutate(cor_accrohn_test, sig = ifelse(p <0.05, "Sig.", "Non Sig."))
ggplot()+geom_tile(data=cor_accrohn_sig,aes(var1, var2, fill = cor),color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
        limit = c(-1,1), space = "Lab", name="Correlation") +
  geom_tile(data = filter(cor_accrohn_sig, sig == "Sig."), aes(var1, var2),
        linewidth = 1, colour = "black", fill = "transpa rent")+
  geom_text(data = cor_accrohn_sig, aes(var1, var2, label = round(cor, 2),
                fontface = ifelse(sig == "Sig.", "bold", "plain")))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
    ggtitle("Heatmap - Active Crohn")+xlab("Bacteria")+ylab("Fungi")+
    coord_fixed()

#Quiescent Crohn
cor_quicrohn_test<-as.data.frame(bac_fungi_quicrohn)%>%cor_test(method = "spearman")%>%
    filter(!((var1%in%bac_quicrohn$Genus)&(var2%in%bac_quicrohn$Genus)))%>%
    filter(!((var1%in%fungi_quicrohn$Genus)&(var2%in%fungi_quicrohn$Genus)))%>%
    mutate(var1_factor = factor(var1, levels = unique(c(var1, var2))),
         var2_factor = factor(var2, levels = unique(c(var1, var2))))%>%
    filter(as.numeric(var1_factor) < as.numeric(var2_factor))

cor_quicrohn_test%>%ggplot(aes(var1, var2, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
            limit = c(-1,1), space = "Lab", name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
  ggtitle("Heatmap - Quiescent Crohn")+xlab("Bacteria")+ylab("Fungi")+
  coord_fixed()

cor_quicrohn_sig<-mutate(cor_quicrohn_test, sig = ifelse(p <0.05, "Sig.", "Non Sig."))
ggplot()+geom_tile(data=cor_quicrohn_sig,aes(var1, var2, fill = cor),color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
        limit = c(-1,1), space = "Lab", name="Correlation") +
  geom_tile(data = filter(cor_quicrohn_sig, sig == "Sig."), aes(var1, var2),
        linewidth = 1, colour = "black", fill = "transpa rent")+
  geom_text(data = cor_quicrohn_sig, aes(var1, var2, label = round(cor, 2),
                fontface = ifelse(sig == "Sig.", "bold", "plain")))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    plot.title = element_text(hjust = 0.5,size=20,face = "bold")) +
    ggtitle("Heatmap - Quiescent Crohn")+xlab("Bacteria")+ylab("Fungi")+
    coord_fixed()