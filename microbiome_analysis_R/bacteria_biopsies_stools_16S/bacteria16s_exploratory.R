##################################################################################################################
#The aim of this code is to create a stacked bar plot with the 15 more represented genus of bacterias in our data#
#############Also, at the end of this file both ANOSIM and PERMANOVA tests will be performed######################
##################################################################################################################
#First of all, the required libraries will be loaded: to read data, to process it and to plot it
library(readr)
library(dplyr)
library(ggplot2)

#The working directory where all of our datasets are will be set. This way it is not necessary to give the full path when loading the data
setwd("path_to_where_your_files_are")

#Both the metadata file as well as the taxonomy matrix are loaded
meta_data <- read_delim("metadatamerge.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
tax_mat <- read_delim("tax-mat.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)

#If by any chance you delete and can not recover your otu matrix as a csv file, it can be redone using the otu matrix downloaded from qiime2
#metadata_3_otu <- read_delim("otu_mat_qiime.tsv",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#metadata_3_otu<-metadata_3_otu[-1,] #the first row is deleted as it just contains information about the columns
#metadata_3_otut<-t(metadata_3_otu) #we need the columns as rows so we transpose it
#colotu<-metadata_3_otut[1,] #we save the first row as it will be our column names
#metadata_3_otut<-metadata_3_otut[-1,] #we delete it as it has been saved
#colnames(metadata_3_otut)<-colotu #we give the column names the saved vector
#write.csv(metadata_3_otut,file = "otu_mat.csv") #The object that was created us saved as a csv file for being able to use it for this and the other code

#Once the otu matrix is in a csv it is loaded
otu_mat <- read_csv("otu_mat.csv")
colnames(otu_mat)[1]<-"Feature ID" #the first column does not have a name, it is set to be the same as in the taxonomy matrix, it refers to the column containing the otu identifiers
#The otu matrix is split into a biopsies otu matrix and a stools otu matrix. The first column is needed as it contains the otus and the rest of the columns that are being kept as they correspond to the samples of biopsies or stools
#Keep in mind that this could have been done differently, for example extracting the sample names for biopsies/stools from the metadata matrix and then giving that as well as the first column name ("Feature ID") to get the same output
otu_mat_biop<-otu_mat[,1:26]
otu_mat_st<-otu_mat[,c(1,27:56)]

#I would like to add a disclaimer from now on: There is probably another way to get to the same output, maybe even easier or faster but this is how it was done

#With these 5 lines we will have the names of the samples for each IBD_type and a column to say if they correspond to biopsies or stools
control<-meta_data%>%filter(IBD_type=="CONTROL")%>%select(`sample-id`,Sample_type)
acuc<-meta_data%>%filter(IBD_type=="ACTIVE UC")%>%select(`sample-id`,Sample_type)
quiuc<-meta_data%>%filter(IBD_type=="QUIESCENT UC")%>%select(`sample-id`,Sample_type)
accrohn<-meta_data%>%filter(IBD_type=="ACTIVE CROHN")%>%select(`sample-id`,Sample_type)
quicrohn<-meta_data%>%filter(IBD_type=="QUIESCENT CROHN")%>%select(`sample-id`,Sample_type)

#It will be first done for Biopsies, then for stools

#By doing this the biopsies otu matrix will be merged with the taxa one by the column "Feature ID". 
#Then, using select the columns that are not needed will be deleted, all the ones present in the taxonomy matrix minus Genus, which is the 7th column in the taxonomy matrix, if you wanted to use Family instead of Genus you would have put a 6 instead of the 7
otu_biop_tax<-otu_mat_biop%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
#A new variable for the mean within columns for each type of sample will be created. That is why the columns corresponding to the name of the samples calculated within lines 38 to 42 are selected
otu_biop_tax$mean_control<-rowMeans(select(otu_biop_tax,select = control$`sample-id`[control$Sample_type=="Biopsies"]))
otu_biop_tax$mean_acuc<-rowMeans(subset(otu_biop_tax,select = acuc$`sample-id`[acuc$Sample_type=="Biopsies"]))
otu_biop_tax$mean_quiuc<-rowMeans(subset(otu_biop_tax,select = quiuc$`sample-id`[quiuc$Sample_type=="Biopsies"]))
otu_biop_tax$mean_accrohn<-rowMeans(subset(otu_biop_tax,select = accrohn$`sample-id`[accrohn$Sample_type=="Biopsies"]))
otu_biop_tax$mean_quicrohn<-rowMeans(subset(otu_biop_tax,select = quicrohn$`sample-id`[quicrohn$Sample_type=="Biopsies"]))

#The columns that begin with "GEI" will be deleted as they correspond to the names of the samples and they are no longer needed
otu_biop_tax<-otu_biop_tax%>%select(-starts_with("GEI"))

#A new dataframe containing 3 columns is needed: IBD_type, Genus and the mean value. Therefore, our previous dataframe will be divided into 5 different dataframes, one for each level of IBD_type so that the column with the IBD_type can be created
#The Genus column will be kept and the one corresponding to the mean of the group as well. Then, the dataset will be grouped by Genus, adding the ones that are the same and deleting the ones with 0 as the result value
biop_control<-otu_biop_tax%>%select(Genus,mean_control)%>%group_by(Genus)%>%summarise_at(vars(mean_control),list(mean_value=sum))%>%filter(mean_value!=0)
biop_acuc<-otu_biop_tax%>%select(Genus,mean_acuc)%>%group_by(Genus)%>%summarise_at(vars(mean_acuc),list(mean_value=sum))%>%filter(mean_value!=0)
biop_quiuc<-otu_biop_tax%>%select(Genus,mean_quiuc)%>%group_by(Genus)%>%summarise_at(vars(mean_quiuc),list(mean_value=sum))%>%filter(mean_value!=0)
biop_accrohn<-otu_biop_tax%>%select(Genus,mean_accrohn)%>%group_by(Genus)%>%summarise_at(vars(mean_accrohn),list(mean_value=sum))%>%filter(mean_value!=0)
biop_quicrohn<-otu_biop_tax%>%select(Genus,mean_quicrohn)%>%group_by(Genus)%>%summarise_at(vars(mean_quicrohn),list(mean_value=sum))%>%filter(mean_value!=0)

#A new column for the type is created, by giving just the text it puts it in all of the observations of the dataset
biop_control$IBD_type<-"CONTROL"
biop_acuc$IBD_type<-"ACTIVE UC"
biop_quiuc$IBD_type<-"QUIESCENT UC"
biop_accrohn$IBD_type<-"ACTIVE CROHN"
biop_quicrohn$IBD_type<-"QUIESCENT CROHN"

#This is just to store the percentages that are being represented and the percentage of NA's that are not being considered when getting the top 15
percentage_biop<-matrix(rep(0,10),ncol = 5)
colnames(percentage_biop)<-c("CONTROL","ACTIVE UC","QUIESCENT UC","ACTIVE CROHN","QUIESCENT CROHN")
rownames(percentage_biop)<-c("Represented","NA not used")

#The dataset will be ordered from highest frequency to lowest and the top 15 without taking into account the NA's will be kept
#Warning, we are keeping the top 15, standardizing the values to 0-1 that is why is so important to say the percentage that is being represented.

#For control
biop_control_top<-biop_control%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_biop["Represented","CONTROL"]<-sum(biop_control_top$mean_value)
percentage_biop["NA not used","CONTROL"]<-as.numeric(biop_control%>%filter(is.na(Genus))%>%select(mean_value))
biop_control_top$mean_value<-biop_control_top$mean_value/percentage_biop["Represented","CONTROL"]
#For active uc
biop_acuc_top<-biop_acuc%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_biop["Represented","ACTIVE UC"]<-sum(biop_acuc_top$mean_value)
percentage_biop["NA not used","ACTIVE UC"]<-as.numeric(biop_acuc%>%filter(is.na(Genus))%>%select(mean_value))
biop_acuc_top$mean_value<-biop_acuc_top$mean_value/percentage_biop["Represented","ACTIVE UC"]
#For quiescent uc
biop_quiuc_top<-biop_quiuc%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_biop["Represented","QUIESCENT UC"]<-sum(biop_quiuc_top$mean_value)
percentage_biop["NA not used","QUIESCENT UC"]<-as.numeric(biop_quiuc%>%filter(is.na(Genus))%>%select(mean_value))
biop_quiuc_top$mean_value<-biop_quiuc_top$mean_value/percentage_biop["Represented","QUIESCENT UC"]
#For active crohn
biop_accrohn_top<-biop_accrohn%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_biop["Represented","ACTIVE CROHN"]<-sum(biop_accrohn_top$mean_value)
percentage_biop["NA not used","ACTIVE CROHN"]<-as.numeric(biop_accrohn%>%filter(is.na(Genus))%>%select(mean_value))
biop_accrohn_top$mean_value<-biop_accrohn_top$mean_value/percentage_biop["Represented","ACTIVE CROHN"]
#For quiescent crohn
biop_quicrohn_top<-biop_quicrohn%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_biop["Represented","QUIESCENT CROHN"]<-sum(biop_quicrohn_top$mean_value)
percentage_biop["NA not used","QUIESCENT CROHN"]<-as.numeric(biop_quicrohn%>%filter(is.na(Genus))%>%select(mean_value))
biop_quicrohn_top$mean_value<-biop_quicrohn_top$mean_value/percentage_biop["Represented","QUIESCENT CROHN"]

#The datasets are merged, one below the other
biop_represent<-union(biop_control_top,biop_acuc_top)
biop_represent<-union(biop_represent,biop_quiuc_top)
biop_represent<-union(biop_represent,biop_accrohn_top)
biop_represent<-union(biop_represent,biop_quicrohn_top)
biop_represent<-biop_represent%>%arrange(Genus)#It is ordered from least to most

#It is used to check how many different Genus are there in the dataset that is used to represent, 36. 
#Therefore, 36 colours that can be differentiated are chosen. There are libraries that can provide them but we chose to do it manually as it was hard to see the differences with the automatic ones 
length(unique(biop_represent$Genus))
our_color_scale<-c("#000000","#FFCC00","#0099FF","#FF00CC","#99CC00","#FF3300",
             "#000099","#FFCCFF","#666666","#FFFF00","#006633","#CCFFFF",
             "#330000","#9900FF","#99FF66","#6666FF","#CCCC33","#00CC99",
             "#FFFFCC","#996699","#CC0000","#66FFFF","#FF9933","#003366",
             "#F2F2F2","#663300","#0066CC","#FF6699","#E5C494","#6A3D9A",
             "#B3E2CD","#FC8D62","#A6CEE3","#4DAF4A","#8DA0CB","#F781BF")
#From left to right, all the parameters that are given: the dataset containing 3 columns (genus, mean_value and IBD_type), 
#aes is a function containing fill=Genus so that it is filled with a different color according to the Genus, with y=mean_value it is being said that the y-axis represents the relative frequency calculated and with x=IBD_type that it has the a different column for each IBD_type
#With geom_bar it is made a stacked barplot from the data given before, with position=fill it makes it stacked and with stat=identity it makes it into 0-1, with color=black the lines between the stacks are black
#With ggtitle the title is set
#With theme some parametres of the text of the title are being modified. hjust=0.5 centers the text, size=20 gives the title a size of 20 and face=bold makes the text bold 
biopplot<-ggplot(biop_represent,aes(fill=Genus,y=mean_value,x=IBD_type))+
  geom_bar(position = "fill",stat="identity",color="black")+
  scale_fill_manual(values=our_color_scale)+
  ggtitle("Top 15 Genus of Bacterias in Biopsies")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"))
x11()
biopplot

#It is now done for stools, it is the same
otu_stools_tax<-otu_mat_st%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
#Use Sample_type==Stools instead of Biopsies
otu_stools_tax$mean_control<-rowMeans(select(otu_stools_tax,select = control$`sample-id`[control$Sample_type=="Stools"]))
otu_stools_tax$mean_acuc<-rowMeans(subset(otu_stools_tax,select = acuc$`sample-id`[acuc$Sample_type=="Stools"]))
otu_stools_tax$mean_quiuc<-rowMeans(subset(otu_stools_tax,select = quiuc$`sample-id`[quiuc$Sample_type=="Stools"]))
otu_stools_tax$mean_accrohn<-rowMeans(subset(otu_stools_tax,select = accrohn$`sample-id`[accrohn$Sample_type=="Stools"]))
otu_stools_tax$mean_quicrohn<-rowMeans(subset(otu_stools_tax,select = quicrohn$`sample-id`[quicrohn$Sample_type=="Stools"]))
#To delete the sample names of stools the ones starting with MS are deleted
otu_stools_tax<-otu_stools_tax%>%select(-starts_with("MS"))
#Same
stools_control<-otu_stools_tax%>%select(Genus,mean_control)%>%group_by(Genus)%>%summarise_at(vars(mean_control),list(mean_value=sum))%>%filter(mean_value!=0)
stools_acuc<-otu_stools_tax%>%select(Genus,mean_acuc)%>%group_by(Genus)%>%summarise_at(vars(mean_acuc),list(mean_value=sum))%>%filter(mean_value!=0)
stools_quiuc<-otu_stools_tax%>%select(Genus,mean_quiuc)%>%group_by(Genus)%>%summarise_at(vars(mean_quiuc),list(mean_value=sum))%>%filter(mean_value!=0)
stools_accrohn<-otu_stools_tax%>%select(Genus,mean_accrohn)%>%group_by(Genus)%>%summarise_at(vars(mean_accrohn),list(mean_value=sum))%>%filter(mean_value!=0)
stools_quicrohn<-otu_stools_tax%>%select(Genus,mean_quicrohn)%>%group_by(Genus)%>%summarise_at(vars(mean_quicrohn),list(mean_value=sum))%>%filter(mean_value!=0)

percentage_st<-matrix(rep(0,10),ncol = 5)
colnames(percentage_st)<-c("CONTROL","ACTIVE UC","QUIESCENT UC","ACTIVE CROHN","QUIESCENT CROHN")
rownames(percentage_st)<-c("Represented","NA not used")

#Order it most to less frequency and get the top 15
#Warning, same as before
stools_control_top<-stools_control%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_st["Represented","CONTROL"]<-sum(stools_control_top$mean_value)
percentage_st["NA not used","CONTROL"]<-as.numeric(stools_control%>%filter(is.na(Genus))%>%select(mean_value))
stools_control_top$mean_value<-stools_control_top$mean_value/percentage_st["Represented","CONTROL"]

stools_acuc_top<-stools_acuc%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_st["Represented","ACTIVE UC"]<-sum(stools_acuc_top$mean_value)
percentage_st["NA not used","ACTIVE UC"]<-as.numeric(stools_acuc%>%filter(is.na(Genus))%>%select(mean_value))
stools_acuc_top$mean_value<-stools_acuc_top$mean_value/percentage_st["Represented","ACTIVE UC"]

stools_quiuc_top<-stools_quiuc%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_st["Represented","QUIESCENT UC"]<-sum(stools_quiuc_top$mean_value)
percentage_st["NA not used","QUIESCENT UC"]<-as.numeric(stools_quiuc%>%filter(is.na(Genus))%>%select(mean_value))
stools_quiuc_top$mean_value<-stools_quiuc_top$mean_value/percentage_st["Represented","QUIESCENT UC"]

stools_accrohn_top<-stools_accrohn%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_st["Represented","ACTIVE CROHN"]<-sum(stools_accrohn_top$mean_value)
percentage_st["NA not used","ACTIVE CROHN"]<-as.numeric(stools_accrohn%>%filter(is.na(Genus))%>%select(mean_value))
stools_accrohn_top$mean_value<-stools_accrohn_top$mean_value/percentage_st["Represented","ACTIVE CROHN"]

stools_quicrohn_top<-stools_quicrohn%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentage_st["Represented","QUIESCENT CROHN"]<-sum(stools_quicrohn_top$mean_value)
percentage_st["NA not used","QUIESCENT CROHN"]<-as.numeric(stools_quicrohn%>%filter(is.na(Genus))%>%select(mean_value))
stools_quicrohn_top$mean_value<-stools_quicrohn_top$mean_value/percentage_st["Represented","QUIESCENT CROHN"]

#Give the column with the type of IBD
stools_control_top$IBD_type<-"CONTROL"
stools_acuc_top$IBD_type<-"ACTIVE UC"
stools_quiuc_top$IBD_type<-"QUIESCENT UC"
stools_accrohn_top$IBD_type<-"ACTIVE CROHN"
stools_quicrohn_top$IBD_type<-"QUIESCENT CROHN"

#Union of all the datasets into one to create the stacked bar plot
stools_represent<-union(stools_control_top,stools_acuc_top)
stools_represent<-union(stools_represent,stools_quiuc_top)
stools_represent<-union(stools_represent,stools_accrohn_top)
stools_represent<-union(stools_represent,stools_quicrohn_top)
#Exactly the same as before
stoolsplot<-ggplot(stools_represent,aes(fill=Genus,y=mean_value,x=IBD_type))+
  geom_bar(position = "fill",stat="identity",color="black")+
  scale_fill_manual(values=our_color_scale) + 
  ggtitle("Top 15 Genus of Bacterias in Stools") +
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"))
x11()
stoolsplot


#####################################################################
####              ANOSIM  &  PERMANOVA                            ###
#####################################################################
#To do the ANOSIM test, the library vegan needs to be loaded
library(vegan)
#For the PERMANOVA test, the library ggvegan is needed
library(ggvegan)

#To run both the anosim as well as permanova tests the absolute frequencies matrix is needed, therefore it is read into the environment
freq_abs_otu<- read_table("freq_abs_16S_biop_stools.csv")
colnames(freq_abs_otu)[1]<-"Feature ID"

###################################################################################################################################
#The dataset will be divided into the 5 different levels of IBD_type and then check within that group, between biopsies and stools#
###################################################################################################################################
#It will be done for each level of IBD_type, it will be explained for control and for the rest is the same procedure

#Control
#The otus column and the ones corresponding to control samples are selected. Then the result of that is joined with the taxonomy matrix by the column of otus and
#the columns corresponding to the taxonomy matrix, except for the one of Genus (7) are deleted as they are no longer needed
otu_control_genus<-freq_abs_otu%>%select(c(1,control$`sample-id`))%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
#The datatset is grouped by Genus and the columns that refer to samples are summarized, as done in the first lines of this same code for the plots
otu_control_genus<-otu_control_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-3JCRR-16SM`:MS12350),list(sum))
#The dataset is transposed. When doing this some of the structure is affected, that is why we save the column names into a vector, then delete the row corresponding to it and set them back as column names
#Also, the row names are preserved as they are a column itself, that will then be added as "sample-id" in line 239
#It is important to note that the colunm names only correspond to the first 349 columns as it is the number of Genus that are in this dataset, check the size of controlcolsnames to know how many you have if you try to replicate this code
#Moreover, when setting the rownames to just the number they correspond to, be aware that 10 is because there are 10 samples of control between stools and biopsies, check the size of controlfirstcol
otu_control_genust<-t(otu_control_genus)
controlcolsnames<-otu_control_genust[1,]
otu_control_genust<-otu_control_genust[-1,]
controlfirstcol<-rownames(otu_control_genust)
rownames(otu_control_genust)<-1:10
otu_control_genust<-as.data.frame(otu_control_genust)
otu_control_genust$`sample-id`<-controlfirstcol
colnames(otu_control_genust)[1:349]<-controlcolsnames
control_genus_anosim<-otu_control_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL",],by="sample-id")

#Both ANOSIM and PERMANOVA need a matrix containing just the absolute frequencies, that is why all rows are kept as well as all the columns corresponding to the Genus
matrixcontrol<-as.matrix(control_genus_anosim[,1:349])
matrixcontrol<-matrix(as.numeric(matrixcontrol),ncol = 349,byrow = FALSE)#When converting it, the rows changed format to character, they are converted back to numericc
#Anosim with several distances, from now on only bray will be performed as the variation in the results is almost none
anosim(matrixcontrol,control_genus_anosim$Sample_type,distance = "bray")
anosim(matrixcontrol,control_genus_anosim$Sample_type,distance = "jaccard")
anosim(matrixcontrol,control_genus_anosim$Sample_type,distance = "chao")
#Permanova, the distance matrix, using the measure of bray will be calculated to perform the permanova test
dist_bray_control<-vegdist(matrixcontrol, method='bray')
adonis2(dist_bray_control~as.factor(control_genus_anosim$Sample_type))

#Active UC
otu_acuc_genus<-freq_abs_otu%>%select(c(1,acuc$`sample-id`))%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_acuc_genus<-otu_acuc_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HURH-4JAAR-16SM`:MS12353),list(sum))
otu_acuc_genust<-t(otu_acuc_genus)
acuccolsnames<-otu_acuc_genust[1,]
otu_acuc_genust<-otu_acuc_genust[-1,]
acucfirstcol<-rownames(otu_acuc_genust)
rownames(otu_acuc_genust)<-1:10
otu_acuc_genust<-as.data.frame(otu_acuc_genust)
otu_acuc_genust$`sample-id`<-acucfirstcol
colnames(otu_acuc_genust)[1:349]<-acuccolsnames
acuc_genus_anosim<-otu_acuc_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC",],by="sample-id")

matrixacuc<-as.matrix(acuc_genus_anosim[,1:349])
matrixacuc<-matrix(as.numeric(matrixacuc),ncol = 349,byrow = FALSE)
#ANOSIM
anosim(matrixacuc,acuc_genus_anosim$Sample_type,distance = "bray")
#Permanova
dist_bray_acuc<-vegdist(matrixacuc,method="bray")
adonis2(dist_bray_acuc~as.factor(acuc_genus_anosim$Sample_type))

#Quiescent UC
otu_quiuc<-freq_abs_otu%>%select(c(1,quiuc$`sample-id`))
otu_quiuc_genus<-otu_quiuc%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_quiuc_genus<-otu_quiuc_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-5YPM-16SM`:MS12357),list(sum))
otu_quiuc_genust<-t(otu_quiuc_genus)
quiuccolsnames<-otu_quiuc_genust[1,]
otu_quiuc_genust<-otu_quiuc_genust[-1,]
quiucfirstcol<-rownames(otu_quiuc_genust)
rownames(otu_quiuc_genust)<-1:12
otu_quiuc_genust<-as.data.frame(otu_quiuc_genust)
otu_quiuc_genust$`sample-id`<-quiucfirstcol
colnames(otu_quiuc_genust)[1:349]<-quiuccolsnames
quiuc_genus_anosim<-otu_quiuc_genust%>%left_join(meta_data[meta_data$IBD_type=="QUIESCENT UC",],by="sample-id")

matrixquiuc<-as.matrix(quiuc_genus_anosim[,1:349])
matrixquiuc<-matrix(as.numeric(matrixquiuc),ncol = 349,byrow = FALSE)
#Anosim
anosim(matrixquiuc,quiuc_genus_anosim$Sample_type,distance = "bray")
#Permanova
dist_bray_quiuc<-vegdist(matrixquiuc,method="bray")
adonis2(dist_bray_quiuc~as.factor(quiuc_genus_anosim$Sample_type))

#Active Crohn
otu_accrohn<-freq_abs_otu%>%select(c(1,accrohn$`sample-id`))
otu_accrohn_genus<-otu_accrohn%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_accrohn_genus<-otu_accrohn_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HURH-8JAMF-16SM`:MS12359),list(sum))
otu_accrohn_genust<-t(otu_accrohn_genus)
accrohncolsnames<-otu_accrohn_genust[1,]
otu_accrohn_genust<-otu_accrohn_genust[-1,]
accrohnfirstcol<-rownames(otu_accrohn_genust)
rownames(otu_accrohn_genust)<-1:11
otu_accrohn_genust<-as.data.frame(otu_accrohn_genust)
otu_accrohn_genust$`sample-id`<-accrohnfirstcol
colnames(otu_accrohn_genust)[1:349]<-accrohncolsnames
accrohn_genus_anosim<-otu_accrohn_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrixaccrohn<-as.matrix(accrohn_genus_anosim[,1:349])
matrixaccrohn<-matrix(as.numeric(matrixaccrohn),ncol = 349,byrow = FALSE)
#Anosim
anosim(matrixaccrohn,accrohn_genus_anosim$Sample_type,distance = "bray")
#Permanova
dist_bray_accrohn<-vegdist(matrixaccrohn,method="bray")
adonis2(dist_bray_accrohn~as.factor(accrohn_genus_anosim$Sample_type))

#Quiescent Crohn
otu_quicrohn<-freq_abs_otu%>%select(c(1,quicrohn$`sample-id`))
otu_quicrohn_genus<-otu_quicrohn%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_quicrohn_genus<-otu_quicrohn_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HURH-10MGB-16SM`:MS12352),list(sum))
otu_quicrohn_genust<-t(otu_quicrohn_genus)
quicrohncolsnames<-otu_quicrohn_genust[1,]
otu_quicrohn_genust<-otu_quicrohn_genust[-1,]
quicrohnfirstcol<-rownames(otu_quicrohn_genust)
rownames(otu_quicrohn_genust)<-1:12
otu_quicrohn_genust<-as.data.frame(otu_quicrohn_genust)
otu_quicrohn_genust$`sample-id`<-quicrohnfirstcol
colnames(otu_quicrohn_genust)[1:349]<-quicrohncolsnames
quicrohn_genus_anosim<-otu_quicrohn_genust%>%left_join(meta_data[meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrixquicrohn<-as.matrix(quicrohn_genus_anosim[,1:349])
matrixquicrohn<-matrix(as.numeric(matrixquicrohn),ncol = 349,byrow = FALSE)
#Anosim
anosim(matrixquicrohn,quicrohn_genus_anosim$Sample_type,distance = "bray")
#Permanova
dist_bray_quicrohn<-vegdist(matrixquicrohn,method="bray")
adonis2(dist_bray_quicrohn~as.factor(quicrohn_genus_anosim$Sample_type))

########################################################################################################################################
#The tests will now be performed dividing the dataset into stools and biopsies, and checking all the different combinations of IBD_type#
########################################################################################################################################
#This 10 lines of code will give the sample names for the different comparisons between the IBD_type levels. 
#Alongside this code some objects will be named using a number, it corresponds to the comparison that is being made in the same order as they can be seen in the following lines
#Also, there are going to be little to no explanations as it is the same procedure as before but for different groups
control_acuc<-meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="ACTIVE UC",1]
control_quiuc<-meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="QUIESCENT UC",1]
control_accrohn<-meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="ACTIVE CROHN",1]
control_quicrohn<-meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="QUIESCENT CROHN",1]
acuc_quiuc<-meta_data[meta_data$IBD_type=="QUIESCENT UC"|meta_data$IBD_type=="ACTIVE UC",1]
acuc_accrohn<-meta_data[meta_data$IBD_type=="ACTIVE CROHN"|meta_data$IBD_type=="ACTIVE UC",1]
acuc_quicrohn<-meta_data[meta_data$IBD_type=="QUIESCENT CROHN"|meta_data$IBD_type=="ACTIVE UC",1]
quiuc_accrohn<-meta_data[meta_data$IBD_type=="QUIESCENT UC"|meta_data$IBD_type=="ACTIVE CROHN",1]
quiuc_quicrohn<-meta_data[meta_data$IBD_type=="QUIESCENT UC"|meta_data$IBD_type=="QUIESCENT CROHN",1]
accrohn_quicrohn<-meta_data[meta_data$IBD_type=="ACTIVE CROHN"|meta_data$IBD_type=="QUIESCENT CROHN",1]

################################################################################
#It will be first done for biopsies

#Control-Active UC
otu_1<-freq_abs_otu%>%select(c(1,control_acuc$`sample-id`))
otu_1_genus<-otu_1%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_1_genus<-otu_1_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-3JCRR-16SM`:`GEICYL-HCU-16-MJCT-16SM`),list(sum))
otu_1_genust<-t(otu_1_genus)
columnsnames1<-otu_1_genust[1,]
otu_1_genust<-otu_1_genust[-1,]
firstcol1<-rownames(otu_1_genust)
rownames(otu_1_genust)<-1:10
otu_1_genust<-as.data.frame(otu_1_genust)
otu_1_genust$`sample-id`<-firstcol1
colnames(otu_1_genust)[1:349]<-columnsnames1
genus_1_anosim<-otu_1_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="ACTIVE UC",],by="sample-id")

matrix1<-as.matrix(genus_1_anosim[,1:349])
matrix1<-matrix(as.numeric(matrix1),ncol = 349,byrow = FALSE)
anosim(matrix1,genus_1_anosim$IBD_type,distance = "bray")
dist_bray_1<-vegdist(matrix1, method='bray')
adonis2(dist_bray_1~as.factor(genus_1_anosim$IBD_type))

#Control-Quiescent UC
otu_2<-freq_abs_otu%>%select(c(1,control_quiuc$`sample-id`))
otu_2_genus<-otu_2%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_2_genus<-otu_2_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-3JCRR-16SM`:`GEICYL-HCU-8AGR-16SM`),list(sum))
otu_2_genust<-t(otu_2_genus)
columnsnames2<-otu_2_genust[1,]
otu_2_genust<-otu_2_genust[-1,]
firstcol2<-rownames(otu_2_genust)
rownames(otu_2_genust)<-1:10
otu_2_genust<-as.data.frame(otu_2_genust)
otu_2_genust$`sample-id`<-firstcol2
colnames(otu_2_genust)[1:349]<-columnsnames2
genus_2_anosim<-otu_2_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="QUIESCENT UC",],by="sample-id")

matrix2<-as.matrix(genus_2_anosim[,1:349])
matrix2<-matrix(as.numeric(matrix2),ncol = 349,byrow = FALSE)
anosim(matrix2,genus_2_anosim$IBD_type,distance = "bray")
dist_bray_2<-vegdist(matrix2, method='bray')
adonis2(dist_bray_2~as.factor(genus_2_anosim$IBD_type))

#Control-Active Crohn
otu_3<-freq_abs_otu%>%select(c(1,control_accrohn$`sample-id`))
otu_3_genus<-otu_3%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_3_genus<-otu_3_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-3JCRR-16SM`:`GEICYL-HURH-18JLGSJ-16SM`),list(sum))
otu_3_genust<-t(otu_3_genus)
columnsnames3<-otu_3_genust[1,]
otu_3_genust<-otu_3_genust[-1,]
firstcol3<-rownames(otu_3_genust)
rownames(otu_3_genust)<-1:10
otu_3_genust<-as.data.frame(otu_3_genust)
otu_3_genust$`sample-id`<-firstcol3
colnames(otu_3_genust)[1:349]<-columnsnames3
genus_3_anosim<-otu_3_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrix3<-as.matrix(genus_3_anosim[,1:349])
matrix3<-matrix(as.numeric(matrix3),ncol = 349,byrow = FALSE)
anosim(matrix3,genus_3_anosim$IBD_type,distance = "bray")
dist_bray_3<-vegdist(matrix3, method='bray')
adonis2(dist_bray_3~as.factor(genus_3_anosim$IBD_type))

#Control-Quiescent Crohn
otu_4<-freq_abs_otu%>%select(c(1,control_quicrohn$`sample-id`))
otu_4_genus<-otu_4%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_4_genus<-otu_4_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HURH-10MGB-16SM`:`GEICYL-HURH-9CAG-16SM`),list(sum))
otu_4_genust<-t(otu_4_genus)
columnsnames4<-otu_4_genust[1,]
otu_4_genust<-otu_4_genust[-1,]
firstcol4<-rownames(otu_4_genust)
rownames(otu_4_genust)<-1:10
otu_4_genust<-as.data.frame(otu_4_genust)
otu_4_genust$`sample-id`<-firstcol4
colnames(otu_4_genust)[1:349]<-columnsnames4
genus_4_anosim<-otu_4_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix4<-as.matrix(genus_4_anosim[,1:349])
matrix4<-matrix(as.numeric(matrix4),ncol = 349,byrow = FALSE)
anosim(matrix4,genus_4_anosim$IBD_type,distance = "bray")
dist_bray_4<-vegdist(matrix4, method='bray')
adonis2(dist_bray_4~as.factor(genus_4_anosim$IBD_type))

#Active UC-Quiescent UC
otu_5<-freq_abs_otu%>%select(c(1,acuc_quiuc$`sample-id`))
otu_5_genus<-otu_5%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_5_genus<-otu_5_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HURH-4JAAR-16SM`:`GEICYL-HCU-8AGR-16SM`),list(sum))
otu_5_genust<-t(otu_5_genus)
columnsnames5<-otu_5_genust[1,]
otu_5_genust<-otu_5_genust[-1,]
firstcol5<-rownames(otu_5_genust)
rownames(otu_5_genust)<-1:10
otu_5_genust<-as.data.frame(otu_5_genust)
otu_5_genust$`sample-id`<-firstcol5
colnames(otu_5_genust)[1:349]<-columnsnames5
genus_5_anosim<-otu_5_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC"|meta_data$IBD_type=="QUIESCENT UC",],by="sample-id")

matrix5<-as.matrix(genus_5_anosim[,1:349])
matrix5<-matrix(as.numeric(matrix5),ncol = 349,byrow = FALSE)
anosim(matrix5,genus_5_anosim$IBD_type,distance = "bray")
dist_bray_5<-vegdist(matrix5, method='bray')
adonis2(dist_bray_5~as.factor(genus_5_anosim$IBD_type))

#Active UC-Active Crohn
otu_6<-freq_abs_otu%>%select(c(1,acuc_accrohn$`sample-id`))
otu_6_genus<-otu_6%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_6_genus<-otu_6_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HURH-4JAAR-16SM`:`GEICYL-HURH-18JLGSJ-16SM`),list(sum))
otu_6_genust<-t(otu_6_genus)
columnsnames6<-otu_6_genust[1,]
otu_6_genust<-otu_6_genust[-1,]
firstcol6<-rownames(otu_6_genust)
rownames(otu_6_genust)<-1:10
otu_6_genust<-as.data.frame(otu_6_genust)
otu_6_genust$`sample-id`<-firstcol6
colnames(otu_6_genust)[1:349]<-columnsnames6
genus_6_anosim<-otu_6_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC"|meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrix6<-as.matrix(genus_6_anosim[,1:349])
matrix6<-matrix(as.numeric(matrix6),ncol = 349,byrow = FALSE)
anosim(matrix6,genus_6_anosim$IBD_type,distance = "bray")
dist_bray_6<-vegdist(matrix6, method='bray')
adonis2(dist_bray_6~as.factor(genus_6_anosim$IBD_type))

#Active UC-Quiescent Crohn
otu_7<-freq_abs_otu%>%select(c(1,acuc_quicrohn$`sample-id`))
otu_7_genus<-otu_7%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_7_genus<-otu_7_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HURH-10MGB-16SM`:`GEICYL-HURH-9CAG-16SM`),list(sum))
otu_7_genust<-t(otu_7_genus)
columnsnames7<-otu_7_genust[1,]
otu_7_genust<-otu_7_genust[-1,]
firstcol7<-rownames(otu_7_genust)
rownames(otu_7_genust)<-1:10
otu_7_genust<-as.data.frame(otu_7_genust)
otu_7_genust$`sample-id`<-firstcol7
colnames(otu_7_genust)[1:349]<-columnsnames7
genus_7_anosim<-otu_7_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix7<-as.matrix(genus_7_anosim[,1:349])
matrix7<-matrix(as.numeric(matrix7),ncol = 349,byrow = FALSE)
anosim(matrix7,genus_7_anosim$IBD_type,distance = "bray")
dist_bray_7<-vegdist(matrix7, method='bray')
adonis2(dist_bray_7~as.factor(genus_7_anosim$IBD_type))

#Quiescent UC-Active Crohn
otu_8<-freq_abs_otu%>%select(c(1,quiuc_accrohn$`sample-id`))
otu_8_genus<-otu_8%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_8_genus<-otu_8_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-5YPM-16SM`:`GEICYL-HURH-18JLGSJ-16SM`),list(sum))
otu_8_genust<-t(otu_8_genus)
columnsnames8<-otu_8_genust[1,]
otu_8_genust<-otu_8_genust[-1,]
firstcol8<-rownames(otu_8_genust)
rownames(otu_8_genust)<-1:10
otu_8_genust<-as.data.frame(otu_8_genust)
otu_8_genust$`sample-id`<-firstcol8
colnames(otu_8_genust)[1:349]<-columnsnames8
genus_8_anosim<-otu_8_genust%>%left_join(meta_data[meta_data$IBD_type=="QUIESCENT UC"|meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrix8<-as.matrix(genus_8_anosim[,1:349])
matrix8<-matrix(as.numeric(matrix8),ncol = 349,byrow = FALSE)
anosim(matrix8,genus_8_anosim$IBD_type,distance = "bray")
dist_bray_8<-vegdist(matrix8, method='bray')
adonis2(dist_bray_8~as.factor(genus_8_anosim$IBD_type))

#Quiescent UC-Quiescent Crohn
otu_9<-freq_abs_otu%>%select(c(1,quiuc_quicrohn$`sample-id`))
otu_9_genus<-otu_9%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_9_genus<-otu_9_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HURH-10MGB-16SM`:`GEICYL-HURH-9CAG-16SM`),list(sum))
otu_9_genust<-t(otu_9_genus)
columnsnames9<-otu_9_genust[1,]
otu_9_genust<-otu_9_genust[-1,]
firstcol9<-rownames(otu_9_genust)
rownames(otu_9_genust)<-1:10
otu_9_genust<-as.data.frame(otu_9_genust)
otu_9_genust$`sample-id`<-firstcol9
colnames(otu_9_genust)[1:349]<-columnsnames9
genus_9_anosim<-otu_9_genust%>%left_join(meta_data[meta_data$IBD_type=="QUIESCENT UC"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix9<-as.matrix(genus_9_anosim[,1:349])
matrix9<-matrix(as.numeric(matrix9),ncol = 349,byrow = FALSE)
anosim(matrix9,genus_9_anosim$IBD_type,distance = "bray")
dist_bray_9<-vegdist(matrix9, method='bray')
adonis2(dist_bray_9~as.factor(genus_9_anosim$IBD_type))

#Active Crohn-Quiescent Crohn
otu_10<-freq_abs_otu%>%select(c(1,accrohn_quicrohn$`sample-id`))
otu_10_genus<-otu_10%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_10_genus<-otu_10_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HURH-10MGB-16SM`:`GEICYL-HURH-9CAG-16SM`),list(sum))
otu_10_genust<-t(otu_10_genus)
columnsnames10<-otu_10_genust[1,]
otu_10_genust<-otu_10_genust[-1,]
firstcol10<-rownames(otu_10_genust)
rownames(otu_10_genust)<-1:10
otu_10_genust<-as.data.frame(otu_10_genust)
otu_10_genust$`sample-id`<-firstcol10
colnames(otu_10_genust)[1:349]<-columnsnames10
genus_10_anosim<-otu_10_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE CROHN"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix10<-as.matrix(genus_10_anosim[,1:349])
matrix10<-matrix(as.numeric(matrix10),ncol = 349,byrow = FALSE)
anosim(matrix10,genus_10_anosim$IBD_type,distance = "bray")
dist_bray_10<-vegdist(matrix10, method='bray')
adonis2(dist_bray_10~as.factor(genus_10_anosim$IBD_type))

################################################################################
#Now for stools,the names of the objects will be repeated so beware of it before executing this section of code if you want to save the above objects

#Control-Active UC
otu_1<-freq_abs_otu%>%select(c(1,control_acuc$`sample-id`))
otu_1_genus<-otu_1%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_1_genus<-otu_1_genus%>%group_by(Genus)%>%summarise_at(vars(MS12330:MS12353),list(sum))
otu_1_genust<-t(otu_1_genus)
columnsnames1<-otu_1_genust[1,]
otu_1_genust<-otu_1_genust[-1,]
firstcol1<-rownames(otu_1_genust)
rownames(otu_1_genust)<-1:10
otu_1_genust<-as.data.frame(otu_1_genust)
otu_1_genust$`sample-id`<-firstcol1
colnames(otu_1_genust)[1:349]<-columnsnames1
genus_1_anosim<-otu_1_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="ACTIVE UC",],by="sample-id")

matrix1<-as.matrix(genus_1_anosim[,1:349])
matrix1<-matrix(as.numeric(matrix1),ncol = 349,byrow = FALSE)
anosim(matrix1,genus_1_anosim$IBD_type,distance = "bray")
dist_bray_1<-vegdist(matrix1, method='bray')
adonis2(dist_bray_1~as.factor(genus_1_anosim$IBD_type))

#Control-Quiescent UC
otu_2<-freq_abs_otu%>%select(c(1,control_quiuc$`sample-id`))
otu_2_genus<-otu_2%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_2_genus<-otu_2_genus%>%group_by(Genus)%>%summarise_at(vars(MS12331:MS12357),list(sum))
otu_2_genust<-t(otu_2_genus)
columnsnames2<-otu_2_genust[1,]
otu_2_genust<-otu_2_genust[-1,]
firstcol2<-rownames(otu_2_genust)
rownames(otu_2_genust)<-1:12
otu_2_genust<-as.data.frame(otu_2_genust)
otu_2_genust$`sample-id`<-firstcol2
colnames(otu_2_genust)[1:349]<-columnsnames2
genus_2_anosim<-otu_2_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="QUIESCENT UC",],by="sample-id")

matrix2<-as.matrix(genus_2_anosim[,1:349])
matrix2<-matrix(as.numeric(matrix2),ncol = 349,byrow = FALSE)
anosim(matrix2,genus_2_anosim$IBD_type,distance = "bray")
dist_bray_2<-vegdist(matrix2, method='bray')
adonis2(dist_bray_2~as.factor(genus_2_anosim$IBD_type))

#Control-Active Crohn
otu_3<-freq_abs_otu%>%select(c(1,control_accrohn$`sample-id`))
otu_3_genus<-otu_3%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_3_genus<-otu_3_genus%>%group_by(Genus)%>%summarise_at(vars(MS12331:MS12359),list(sum))
otu_3_genust<-t(otu_3_genus)
columnsnames3<-otu_3_genust[1,]
otu_3_genust<-otu_3_genust[-1,]
firstcol3<-rownames(otu_3_genust)
rownames(otu_3_genust)<-1:11
otu_3_genust<-as.data.frame(otu_3_genust)
otu_3_genust$`sample-id`<-firstcol3
colnames(otu_3_genust)[1:349]<-columnsnames3
genus_3_anosim<-otu_3_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrix3<-as.matrix(genus_3_anosim[,1:349])
matrix3<-matrix(as.numeric(matrix3),ncol = 349,byrow = FALSE)
anosim(matrix3,genus_3_anosim$IBD_type,distance = "bray")
dist_bray_3<-vegdist(matrix3, method='bray')
adonis2(dist_bray_3~as.factor(genus_3_anosim$IBD_type))

#Control-Quiescent Crohn
otu_4<-freq_abs_otu%>%select(c(1,control_quicrohn$`sample-id`))
otu_4_genus<-otu_4%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_4_genus<-otu_4_genus%>%group_by(Genus)%>%summarise_at(vars(MS12331:MS12352),list(sum))
otu_4_genust<-t(otu_4_genus)
columnsnames4<-otu_4_genust[1,]
otu_4_genust<-otu_4_genust[-1,]
firstcol4<-rownames(otu_4_genust)
rownames(otu_4_genust)<-1:12
otu_4_genust<-as.data.frame(otu_4_genust)
otu_4_genust$`sample-id`<-firstcol4
colnames(otu_4_genust)[1:349]<-columnsnames4
genus_4_anosim<-otu_4_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix4<-as.matrix(genus_4_anosim[,1:349])
matrix4<-matrix(as.numeric(matrix4),ncol = 349,byrow = FALSE)
anosim(matrix4,genus_4_anosim$IBD_type,distance = "bray")
dist_bray_4<-vegdist(matrix4, method='bray')
adonis2(dist_bray_4~as.factor(genus_4_anosim$IBD_type))

#Active UC-Quiescent UC
otu_5<-freq_abs_otu%>%select(c(1,acuc_quiuc$`sample-id`))
otu_5_genus<-otu_5%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_5_genus<-otu_5_genus%>%group_by(Genus)%>%summarise_at(vars(MS12330:MS12357),list(sum))
otu_5_genust<-t(otu_5_genus)
columnsnames5<-otu_5_genust[1,]
otu_5_genust<-otu_5_genust[-1,]
firstcol5<-rownames(otu_5_genust)
rownames(otu_5_genust)<-1:12
otu_5_genust<-as.data.frame(otu_5_genust)
otu_5_genust$`sample-id`<-firstcol5
colnames(otu_5_genust)[1:349]<-columnsnames5
genus_5_anosim<-otu_5_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC"|meta_data$IBD_type=="QUIESCENT UC",],by="sample-id")

matrix5<-as.matrix(genus_5_anosim[,1:349])
matrix5<-matrix(as.numeric(matrix5),ncol = 349,byrow = FALSE)
anosim(matrix5,genus_5_anosim$IBD_type,distance = "bray")
anosim(matrix4,genus_4_anosim$IBD_type,distance = "bray")
dist_bray_5<-vegdist(matrix5, method='bray')
adonis2(dist_bray_5~as.factor(genus_5_anosim$IBD_type))

#Active UC-Active Crohn
otu_6<-freq_abs_otu%>%select(c(1,acuc_accrohn$`sample-id`))
otu_6_genus<-otu_6%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_6_genus<-otu_6_genus%>%group_by(Genus)%>%summarise_at(vars(MS12330:MS12359),list(sum))
otu_6_genust<-t(otu_6_genus)
columnsnames6<-otu_6_genust[1,]
otu_6_genust<-otu_6_genust[-1,]
firstcol6<-rownames(otu_6_genust)
rownames(otu_6_genust)<-1:11
otu_6_genust<-as.data.frame(otu_6_genust)
otu_6_genust$`sample-id`<-firstcol6
colnames(otu_6_genust)[1:349]<-columnsnames6
genus_6_anosim<-otu_6_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC"|meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrix6<-as.matrix(genus_6_anosim[,1:349])
matrix6<-matrix(as.numeric(matrix6),ncol = 349,byrow = FALSE)
anosim(matrix6,genus_6_anosim$IBD_type,distance = "bray")
dist_bray_6<-vegdist(matrix6, method='bray')
adonis2(dist_bray_6~as.factor(genus_6_anosim$IBD_type))

#Active UC-Quiescent Crohn
otu_7<-freq_abs_otu%>%select(c(1,acuc_quicrohn$`sample-id`))
otu_7_genus<-otu_7%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_7_genus<-otu_7_genus%>%group_by(Genus)%>%summarise_at(vars(MS12330:MS12353),list(sum))
otu_7_genust<-t(otu_7_genus)
columnsnames7<-otu_7_genust[1,]
otu_7_genust<-otu_7_genust[-1,]
firstcol7<-rownames(otu_7_genust)
rownames(otu_7_genust)<-1:12
otu_7_genust<-as.data.frame(otu_7_genust)
otu_7_genust$`sample-id`<-firstcol7
colnames(otu_7_genust)[1:349]<-columnsnames7
genus_7_anosim<-otu_7_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix7<-as.matrix(genus_7_anosim[,1:349])
matrix7<-matrix(as.numeric(matrix7),ncol = 349,byrow = FALSE)
anosim(matrix7,genus_7_anosim$IBD_type,distance = "bray")
dist_bray_7<-vegdist(matrix7, method='bray')
adonis2(dist_bray_7~as.factor(genus_7_anosim$IBD_type))

#Quiescent UC-Active Crohn
otu_8<-freq_abs_otu%>%select(c(1,quiuc_accrohn$`sample-id`))
otu_8_genus<-otu_8%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_8_genus<-otu_8_genus%>%group_by(Genus)%>%summarise_at(vars(MS12333:MS12359),list(sum))
otu_8_genust<-t(otu_8_genus)
columnsnames8<-otu_8_genust[1,]
otu_8_genust<-otu_8_genust[-1,]
firstcol8<-rownames(otu_8_genust)
rownames(otu_8_genust)<-1:13
otu_8_genust<-as.data.frame(otu_8_genust)
otu_8_genust$`sample-id`<-firstcol8
colnames(otu_8_genust)[1:349]<-columnsnames8
genus_8_anosim<-otu_8_genust%>%left_join(meta_data[meta_data$IBD_type=="QUIESCENT UC"|meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrix8<-as.matrix(genus_8_anosim[,1:349])
matrix8<-matrix(as.numeric(matrix8),ncol = 349,byrow = FALSE)
anosim(matrix8,genus_8_anosim$IBD_type,distance = "bray")
dist_bray_8<-vegdist(matrix8, method='bray')
adonis2(dist_bray_8~as.factor(genus_8_anosim$IBD_type))

#Quiescent UC-Quiescent Crohn
otu_9<-freq_abs_otu%>%select(c(1,quiuc_quicrohn$`sample-id`))
otu_9_genus<-otu_9%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_9_genus<-otu_9_genus%>%group_by(Genus)%>%summarise_at(vars(MS12333:MS12357),list(sum))
otu_9_genust<-t(otu_9_genus)
columnsnames9<-otu_9_genust[1,]
otu_9_genust<-otu_9_genust[-1,]
firstcol9<-rownames(otu_9_genust)
rownames(otu_9_genust)<-1:14
otu_9_genust<-as.data.frame(otu_9_genust)
otu_9_genust$`sample-id`<-firstcol9
colnames(otu_9_genust)[1:349]<-columnsnames9
genus_9_anosim<-otu_9_genust%>%left_join(meta_data[meta_data$IBD_type=="QUIESCENT UC"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix9<-as.matrix(genus_9_anosim[,1:349])
matrix9<-matrix(as.numeric(matrix9),ncol = 349,byrow = FALSE)
anosim(matrix9,genus_9_anosim$IBD_type,distance = "bray")
dist_bray_9<-vegdist(matrix9, method='bray')
adonis2(dist_bray_9~as.factor(genus_9_anosim$IBD_type))

#Active Crohn-Quiescent Crohn
otu_10<-freq_abs_otu%>%select(c(1,accrohn_quicrohn$`sample-id`))
otu_10_genus<-otu_10%>%left_join(tax_mat,by="Feature ID")%>%select(-colnames(tax_mat)[-7])
otu_10_genus<-otu_10_genus%>%group_by(Genus)%>%summarise_at(vars(MS12336:MS12359),list(sum))
otu_10_genust<-t(otu_10_genus)
columnsnames10<-otu_10_genust[1,]
otu_10_genust<-otu_10_genust[-1,]
firstcol10<-rownames(otu_10_genust)
rownames(otu_10_genust)<-1:13
otu_10_genust<-as.data.frame(otu_10_genust)
otu_10_genust$`sample-id`<-firstcol10
colnames(otu_10_genust)[1:349]<-columnsnames10
genus_10_anosim<-otu_10_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE CROHN"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix10<-as.matrix(genus_10_anosim[,1:349])
matrix10<-matrix(as.numeric(matrix10),ncol = 349,byrow = FALSE)
anosim(matrix10,genus_10_anosim$IBD_type,distance = "bray")
dist_bray_10<-vegdist(matrix10, method='bray')
adonis2(dist_bray_10~as.factor(genus_10_anosim$IBD_type))