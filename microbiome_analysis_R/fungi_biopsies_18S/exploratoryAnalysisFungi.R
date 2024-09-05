#######################################################################################################
##Our aim is to create a stacked bar plot with the 15 more represented genus of fungi in our data##
##Also, at the end of this file both ANOSIM and PERMANOVA tests will be performed######################
#######################################################################################################
#First of all, the required libraries will be loaded: to read data, to process it and to plot it
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)

#First of all, the working directory must be changed to where all the data is, so it is not necessary to set the complete path when loading other objects
setwd("your_path_to_where_the_data_is")

#The otu dataset will be loaded, and changed the column name to the same as the taxonomy
uparse_otus <- read_excel("otu_mat.xlsx",col_types = c("text",rep("numeric",25)))
colnames(uparse_otus)[1]<-"Feature ID"
#The taxonomy dataset is loaded
taxonomy <- read_delim("tax_mat.csv",delim = ";", escape_double = FALSE, trim_ws = TRUE)
#The metadata dataset is loaded
meta_data <- read_excel("metadatafungi.xlsx")

#In the otu's dataset, the column with the otu ids will be changed to the Genus that correspond to each otu
#It will be done my joining the otu dataset with the taxonomy one by the column of otus (Feature ID). Then, all the columns except for the one corresponding to Genus, the number 7, will be deleted from the final dataset
otu_tax<-uparse_otus%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])

#With the next 5 lines of code the names of the samples for each IBD_type are obtained
control<-meta_data[meta_data[,2]=="CONTROL",1]
acuc<-meta_data[meta_data[,2]=="ACTIVE UC",1]
quiuc<-meta_data[meta_data[,2]=="QUIESCENT UC",1]
accrohn<-meta_data[meta_data[,2]=="ACTIVE CROHN",1]
quicrohn<-meta_data[meta_data[,2]=="QUIESCENT CROHN",1]

#A new variable for the mean of each group will be created by sub-setting the columns corresponding to the sample names that were obtained before
otu_tax$mean_control<-rowMeans(subset(otu_tax,select = control$`sample-id`))
otu_tax$mean_acuc<-rowMeans(subset(otu_tax,select = acuc$`sample-id`))
otu_tax$mean_quiuc<-rowMeans(subset(otu_tax,select = quiuc$`sample-id`))
otu_tax$mean_accrohn<-rowMeans(subset(otu_tax,select = accrohn$`sample-id`))
otu_tax$mean_quicrohn<-rowMeans(subset(otu_tax,select = quicrohn$`sample-id`))

#The columns beginning with GEI are deleted, as they are the ones corresponding to the sample names and they are no longer needed
otu_tax<-otu_tax%>%select(-starts_with("GEI"))

#Now, a new dataframe containing 3 columns (IBD_type, Genus and the mean) needs to be created
#To do so, the otu_tax dataset will be divided into 5 different datasets, one for each IBD_type level to create the column IBD_type
#In order to obtain a dataset for each of the levels, the first column and the one that corresponds to the mean of the group are selected. The result of that is then agroupated by Genus, 
#it is summarized (adding the ones that correspond to the same Genus, as different otus can get to the same Genus) and the observations that have a frequency of 0 are deleted
fungi_control<-otu_tax%>%select(Genus,mean_control)%>%group_by(Genus)%>%summarise_at(vars(mean_control),list(mean_value=sum))%>%filter(mean_value!=0)
fungi_acuc<-otu_tax%>%select(Genus,mean_acuc)%>%group_by(Genus)%>%summarise_at(vars(mean_acuc),list(mean_value=sum))%>%filter(mean_value!=0)
fungi_quiuc<-otu_tax%>%select(Genus,mean_quiuc)%>%group_by(Genus)%>%summarise_at(vars(mean_quiuc),list(mean_value=sum))%>%filter(mean_value!=0)
fungi_accrohn<-otu_tax%>%select(Genus,mean_accrohn)%>%group_by(Genus)%>%summarise_at(vars(mean_accrohn),list(mean_value=sum))%>%filter(mean_value!=0)
fungi_quicrohn<-otu_tax%>%select(Genus,mean_quicrohn)%>%group_by(Genus)%>%summarise_at(vars(mean_quicrohn),list(mean_value=sum))%>%filter(mean_value!=0)

#The new column containing the IBD_type is created. In R, if the length of the vector is not enough it will restart with the first element of the vector again, by giving just the text it puts the same to all observations
fungi_control$IBD_type<-"CONTROL"
fungi_acuc$IBD_type<-"ACTIVE UC"
fungi_quiuc$IBD_type<-"QUIESCENT UC"
fungi_accrohn$IBD_type<-"ACTIVE CROHN"
fungi_quicrohn$IBD_type<-"QUIESCENT CROHN"

#A dataframe to save the percentages of fungi represented and of NA within the 5 levels of IBD_type
percentages<-as.data.frame(matrix(rep(0,10),ncol = 5))
rownames(percentages)<-c("Represented","NA%")
colnames(percentages)<-c("CONTROL","ACTIVE UC","QUIESCENT UC","ACTIVE CROHN","QUIESCENT CROHN")

#It is now ordered by the mean frequency, from most to least and then keeping the top 15 without taking into account the NA's
#Warning, the top 15 that represents x% will be then standarize to 0-1 to represent the stacked barplot
#For control
fungi_control_top<-fungi_control%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentages["Represented","CONTROL"]<-sum(fungi_control_top$mean_value)
percentages["NA%","CONTROL"]<-as.numeric(fungi_control%>%filter(is.na(Genus))%>%select(mean_value))
fungi_control_top$mean_value<-fungi_control_top$mean_value/percentages["Represented","CONTROL"]
#For active uc
fungi_acuc_top<-fungi_acuc%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentages["Represented","ACTIVE UC"]<-sum(fungi_acuc_top$mean_value)
percentages["NA%","ACTIVE UC"]<-as.numeric(fungi_acuc%>%filter(is.na(Genus))%>%select(mean_value))
fungi_acuc_top$mean_value<-fungi_acuc_top$mean_value/percentages["Represented","ACTIVE UC"]
#For quiescent uc
fungi_quiuc_top<-fungi_quiuc%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentages["Represented","QUIESCENT UC"]<-sum(fungi_quiuc_top$mean_value)
percentages["NA%","QUIESCENT UC"]<-as.numeric(fungi_quiuc%>%filter(is.na(Genus))%>%select(mean_value))
fungi_quiuc_top$mean_value<-fungi_quiuc_top$mean_value/percentages["Represented","QUIESCENT UC"]
#For active crohn
fungi_accrohn_top<-fungi_accrohn%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentages["Represented","ACTIVE CROHN"]<-sum(fungi_accrohn_top$mean_value)
percentages["NA%","ACTIVE CROHN"]<-as.numeric(fungi_accrohn%>%filter(is.na(Genus))%>%select(mean_value))
fungi_accrohn_top$mean_value<-fungi_accrohn_top$mean_value/percentages["Represented","ACTIVE CROHN"]
#For quiescent crohn
fungi_quicrohn_top<-fungi_quicrohn%>%filter(!is.na(Genus))%>%arrange(desc(mean_value))%>%top_n(15,wt=mean_value)
percentages["Represented","QUIESCENT CROHN"]<-sum(fungi_quicrohn_top$mean_value)
percentages["NA%","QUIESCENT CROHN"]<-as.numeric(fungi_quicrohn%>%filter(is.na(Genus))%>%select(mean_value))
fungi_quicrohn_top$mean_value<-fungi_quicrohn_top$mean_value/percentages["Represented","QUIESCENT CROHN"]

#Now, the datasets for each level of IBD_type will be merged one bellow the other
fungi_represent<-union(fungi_control_top,fungi_acuc_top)
fungi_represent<-union(fungi_represent,fungi_quiuc_top)
fungi_represent<-union(fungi_represent,fungi_accrohn_top)
fungi_represent<-union(fungi_represent,fungi_quicrohn_top)

#The stacked barplot, 36 colors where choosen manually so that they are easily distinguishable 
color_choice<-c("#000000","#FFCC00","#0099FF","#FF00CC","#99CC00","#FF3300",
             "#000099","#FFCCFF","#666666","#FFFF00","#006633","#CCFFFF",
             "#330000","#9900FF","#99FF66","#6666FF","#CCCC33","#00CC99",
             "#FFFFCC","#996699","#CC0000","#66FFFF","#FF9933","#003366",
             "#F2F2F2","#663300","#0066CC","#FF6699","#E5C494","#6A3D9A",
             "#B3E2CD","#FC8D62","#A6CEE3","#4DAF4A","#8DA0CB","#F781BF")
x11()
ggplot(fungi_represent,aes(fill=Genus,y=mean_value,x=IBD_type))+
  geom_bar(position = "fill",stat="identity",color="black")+
  scale_fill_manual(values=color_choice)+ggtitle("Top 15 Genus of Fungi in Biopsies")+
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"))


#####################################################################
####              ANOSIM  &  PERMANOVA                            ###
#####################################################################
#To do the ANOSIM test, the library vegan needs to be loaded
library(vegan)
#For the PERMANOVA test, the library ggvegan is needed
library(ggvegan)

freq_abs_otu <- read_table("freq_abs_fungi.csv")
colnames(freq_abs_otu)[1]<-"Feature ID"

#This 10 lines of code will give the sample names for the different comparisons between the IBD_type levels. 
#Alongside this code some objects will be named using a number, it corresponds to the comparison that is being made, in the order as they can be seen in the following lines
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

#It is the same procedure for the 10 of them so there will only be an explanation on the first one
#Control-Active UC
#The otus column and the ones corresponding to control samples are selected
otu_1<-freq_abs_otu%>%select(c(1,control_acuc$`sample-id`))
#The result from before is joined with the taxonomy matrix by the column of otus and then the columns corresponding to the taxonomy matrix, except for the one of Genus (7) are deleted as they are no longer needed
#It is then grouped by Genus and the columns that refer to samples are summarized, as done in the first lines of this same code for the plots
otu_1_genus<-otu_1%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-1`:`GEICYL-HURH-4JAAR`),list(sum))
#The dataset is transposed. When doing this some of the structure is affected, that is why we save the column names into a vector, then delete the row corresponding to it and set them back as column names
#Also, the row names are preserved as they are a column itself, that will then be added as "sample-id" in line 154
#It is important to note that the colunm names only correspond to the first 37 columns as it is the number of Genus that are in this dataset, check the size of controlcolsnames to know how many if you try to replicate this code
#Moreover, when setting the the rownames to just the number they correspond to, be aware that 10 is because there are 10 samples between control and active uc, check the size of controlfirstcol
otu_1_genust<-t(otu_1_genus)
colsnames1<-otu_1_genust[1,]
otu_1_genust<-otu_1_genust[-1,]
firstcol1<-rownames(otu_1_genust)
rownames(otu_1_genust)<-1:10
otu_1_genust<-as.data.frame(otu_1_genust)
otu_1_genust$`sample-id`<-firstcol1
colnames(otu_1_genust)[1:37]<-colsnames1
genus_1_anosim<-otu_1_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="ACTIVE UC",],by="sample-id")

#Both ANOSIM and PERMANOVA need a matrix containing just the absolute frequencies, that is why all rows are kept and all the columns corresponding to the Genus
matrix1<-as.matrix(genus_1_anosim[,1:37])
matrix1<-matrix(as.numeric(matrix1),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix1,genus_1_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_1<-vegdist(matrix1, method='bray')
adonis2(dist_bray_1~as.factor(genus_1_anosim$IBD_type))

#Control-Quiescent UC
otu_2<-freq_abs_otu%>%select(c(1,control_quiuc$`sample-id`))
otu_2_genus<-otu_2%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])
otu_2_genus<-otu_2_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-11LMGM`:`GEICYL-HURH-3PGP`),list(sum))
otu_2_genust<-t(otu_2_genus)
colsnames2<-otu_2_genust[1,]
otu_2_genust<-otu_2_genust[-1,]
firstcol2<-rownames(otu_2_genust)
rownames(otu_2_genust)<-1:10
otu_2_genust<-as.data.frame(otu_2_genust)
otu_2_genust$`sample-id`<-firstcol2
colnames(otu_2_genust)[1:37]<-colsnames2
genus_2_anosim<-otu_2_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="QUIESCENT UC",],by="sample-id")

matrix2<-as.matrix(genus_2_anosim[,1:37])
matrix2<-matrix(as.numeric(matrix2),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix2,genus_2_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_2<-vegdist(matrix2, method='bray')
adonis2(dist_bray_2~as.factor(genus_2_anosim$IBD_type))

#Control-Active Crohn
otu_3<-freq_abs_otu%>%select(c(1,control_accrohn$`sample-id`))
otu_3_genus<-otu_3%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])
otu_3_genus<-otu_3_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-11LMGM`:`GEICYL-HURH-8JAMF`),list(sum))
otu_3_genust<-t(otu_3_genus)
colsnames3<-otu_3_genust[1,]
otu_3_genust<-otu_3_genust[-1,]
firstcol3<-rownames(otu_3_genust)
rownames(otu_3_genust)<-1:10
otu_3_genust<-as.data.frame(otu_3_genust)
otu_3_genust$`sample-id`<-firstcol3
colnames(otu_3_genust)[1:37]<-colsnames3
genus_3_anosim<-otu_3_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrix3<-as.matrix(genus_3_anosim[,1:37])
matrix3<-matrix(as.numeric(matrix3),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix3,genus_3_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_3<-vegdist(matrix3, method='bray')
adonis2(dist_bray_3~as.factor(genus_3_anosim$IBD_type))

#Control-Quiescent UC
otu_4<-freq_abs_otu%>%select(c(1,control_quicrohn$`sample-id`))
otu_4_genus<-otu_4%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])
otu_4_genus<-otu_4_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-11LMGM`:`GEICYL-HURH-9CAG`),list(sum))
otu_4_genust<-t(otu_4_genus)
colsnames4<-otu_4_genust[1,]
otu_4_genust<-otu_4_genust[-1,]
firstcol4<-rownames(otu_4_genust)
rownames(otu_4_genust)<-1:10
otu_4_genust<-as.data.frame(otu_4_genust)
otu_4_genust$`sample-id`<-firstcol4
colnames(otu_4_genust)[1:37]<-colsnames4
genus_4_anosim<-otu_4_genust%>%left_join(meta_data[meta_data$IBD_type=="CONTROL"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix4<-as.matrix(genus_4_anosim[,1:37])
matrix4<-matrix(as.numeric(matrix4),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix4,genus_4_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_4<-vegdist(matrix4, method='bray')
adonis2(dist_bray_4~as.factor(genus_4_anosim$IBD_type))

#Active UC-Quiescent UC
otu_5<-freq_abs_otu%>%select(c(1,acuc_quiuc$`sample-id`))
otu_5_genus<-otu_5%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])
otu_5_genus<-otu_5_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-1`:`GEICYL-HURH-4JAAR`),list(sum))
otu_5_genust<-t(otu_5_genus)
colsnames5<-otu_5_genust[1,]
otu_5_genust<-otu_5_genust[-1,]
firstcol5<-rownames(otu_5_genust)
rownames(otu_5_genust)<-1:10
otu_5_genust<-as.data.frame(otu_5_genust)
otu_5_genust$`sample-id`<-firstcol5
colnames(otu_5_genust)[1:37]<-colsnames5
genus_5_anosim<-otu_5_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC"|meta_data$IBD_type=="QUIESCENT UC",],by="sample-id")

matrix5<-as.matrix(genus_5_anosim[,1:37])
matrix5<-matrix(as.numeric(matrix5),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix5,genus_5_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_5<-vegdist(matrix5, method='bray')
adonis2(dist_bray_5~as.factor(genus_5_anosim$IBD_type))

#Active UC-Active Crohn
otu_6<-freq_abs_otu%>%select(c(1,acuc_accrohn$`sample-id`))
otu_6_genus<-otu_6%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])
otu_6_genus<-otu_6_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-1`:`GEICYL-HURH-8JAMF`),list(sum))
otu_6_genust<-t(otu_6_genus)
colsnames6<-otu_6_genust[1,]
otu_6_genust<-otu_6_genust[-1,]
firstcol6<-rownames(otu_6_genust)
rownames(otu_6_genust)<-1:10
otu_6_genust<-as.data.frame(otu_6_genust)
otu_6_genust$`sample-id`<-firstcol6
colnames(otu_6_genust)[1:37]<-colsnames6
genus_6_anosim<-otu_6_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC"|meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrix6<-as.matrix(genus_6_anosim[,1:37])
matrix6<-matrix(as.numeric(matrix6),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix6,genus_6_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_6<-vegdist(matrix6, method='bray')
adonis2(dist_bray_6~as.factor(genus_6_anosim$IBD_type))

#Active UC-Quiescent Crohn
otu_7<-freq_abs_otu%>%select(c(1,acuc_quicrohn$`sample-id`))
otu_7_genus<-otu_7%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])
otu_7_genus<-otu_7_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-1`:`GEICYL-HURH-9CAG`),list(sum))
otu_7_genust<-t(otu_7_genus)
colsnames7<-otu_7_genust[1,]
otu_7_genust<-otu_7_genust[-1,]
firstcol7<-rownames(otu_7_genust)
rownames(otu_7_genust)<-1:10
otu_7_genust<-as.data.frame(otu_7_genust)
otu_7_genust$`sample-id`<-firstcol7
colnames(otu_7_genust)[1:37]<-colsnames7
genus_7_anosim<-otu_7_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE UC"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix7<-as.matrix(genus_7_anosim[,1:37])
matrix7<-matrix(as.numeric(matrix7),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix7,genus_7_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_7<-vegdist(matrix7, method='bray')
adonis2(dist_bray_7~as.factor(genus_7_anosim$IBD_type))

#Quiescent UC-Active Crohn
otu_8<-freq_abs_otu%>%select(c(1,quiuc_accrohn$`sample-id`))
otu_8_genus<-otu_8%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])
otu_8_genus<-otu_8_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-18ASS`:`GEICYL-HURH-8JAMF`),list(sum))
otu_8_genust<-t(otu_8_genus)
colsnames8<-otu_8_genust[1,]
otu_8_genust<-otu_8_genust[-1,]
firstcol8<-rownames(otu_8_genust)
rownames(otu_8_genust)<-1:10
otu_8_genust<-as.data.frame(otu_8_genust)
otu_8_genust$`sample-id`<-firstcol8
colnames(otu_8_genust)[1:37]<-colsnames8
genus_8_anosim<-otu_8_genust%>%left_join(meta_data[meta_data$IBD_type=="QUIESCENT UC"|meta_data$IBD_type=="ACTIVE CROHN",],by="sample-id")

matrix8<-as.matrix(genus_8_anosim[,1:37])
matrix8<-matrix(as.numeric(matrix8),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix8,genus_8_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_8<-vegdist(matrix8, method='bray')
adonis2(dist_bray_8~as.factor(genus_8_anosim$IBD_type))

#Quiescent UC-Quiescent Crohn
otu_9<-freq_abs_otu%>%select(c(1,quiuc_quicrohn$`sample-id`))
otu_9_genus<-otu_9%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])
otu_9_genus<-otu_9_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-15RVL`:`GEICYL-HURH-9CAG`),list(sum))
otu_9_genust<-t(otu_9_genus)
colsnames9<-otu_9_genust[1,]
otu_9_genust<-otu_9_genust[-1,]
firstcol9<-rownames(otu_9_genust)
rownames(otu_9_genust)<-1:10
otu_9_genust<-as.data.frame(otu_9_genust)
otu_9_genust$`sample-id`<-firstcol9
colnames(otu_9_genust)[1:37]<-colsnames9
genus_9_anosim<-otu_9_genust%>%left_join(meta_data[meta_data$IBD_type=="QUIESCENT UC"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix9<-as.matrix(genus_9_anosim[,1:37])
matrix9<-matrix(as.numeric(matrix9),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix9,genus_9_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_9<-vegdist(matrix9, method='bray')
adonis2(dist_bray_9~as.factor(genus_9_anosim$IBD_type))

#Active Crohn-Quiescent Crohn
otu_10<-freq_abs_otu%>%select(c(1,accrohn_quicrohn$`sample-id`))
otu_10_genus<-otu_10%>%left_join(taxonomy,by="Feature ID")%>%select(-colnames(taxonomy)[-7])
otu_10_genus<-otu_10_genus%>%group_by(Genus)%>%summarise_at(vars(`GEICYL-HCU-15RVL`:`GEICYL-HURH-9CAG`),list(sum))
otu_10_genust<-t(otu_10_genus)
colsnames10<-otu_10_genust[1,]
otu_10_genust<-otu_10_genust[-1,]
firstcol10<-rownames(otu_10_genust)
rownames(otu_10_genust)<-1:10
otu_10_genust<-as.data.frame(otu_10_genust)
otu_10_genust$`sample-id`<-firstcol10
colnames(otu_10_genust)[1:37]<-colsnames10
genus_10_anosim<-otu_10_genust%>%left_join(meta_data[meta_data$IBD_type=="ACTIVE CROHN"|meta_data$IBD_type=="QUIESCENT CROHN",],by="sample-id")

matrix10<-as.matrix(genus_10_anosim[,1:37])
matrix10<-matrix(as.numeric(matrix10),ncol = 37,byrow = FALSE)
#Anosim
anosim(matrix10,genus_10_anosim$IBD_type,distance = "bray")
#Permanova
dist_bray_10<-vegdist(matrix10, method='bray')
adonis2(dist_bray_10~as.factor(genus_10_anosim$IBD_type))
