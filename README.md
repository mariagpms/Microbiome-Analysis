# Microbiome Analysis

In this repository you will find some additional resources to the following article: [Fungal-Bacterial Dysbiosis in IBD: Microbial Biomarkers of Disease Activity](https://doi.org/10.1002/mbo3.70088).

## Author
The code of this project was developed by [María García Pizarro](https://orcid.org/0009-0001-3647-1010) in colaboration with [Elisa Arribas Rodriguez](https://orcid.org/0000-0002-4470-075X) and Francisco Javier López Pascual.
We would also like to acknowledge the valuable insights provided by [Francisco Illanes Álvarez](https://orcid.org/0000-0002-0837-8265). 

## What will you find in this repository
There are two folders, one for fungi and another for bacteria. Moreover, there is a .md file containing the relevant commands to use qiime2 to obtain the files that are needed for the further R analysis. In addition, there is a file, [correlation.R](microbiome_analysis_R/correlation.R), that contains the script used to plot the heatmaps of the correlation between fungi and the most represented bacteri. <br />Type of patients: control, active uc, quiescent uc, active crohn, quiescent crohn.
### Fungi
The data corresponds to fungi found on biopsies of the different types of patients mentioned above.
The folder named as plots, has different plots that come as the result of running the R codes.<br />
Content of the files available:
- [exploratoryAnalysisFungi.R](./microbiome_analysis_R/fungi_biopsies_18S/exploratoryAnalysisFungi.R) contains the script we have used to plot the stacked bar plot of the 15 top Genus of Fungi found within the different types of patients. It also contains the script for running both ANOSIM and PERMANOVA tests on our data, checking the differences between all 5 types of patients with each other.
- [lefseAnnalysisFungi.R](./microbiome_analysis_R/fungi_biopsies_18S/lefseAnnalysisFungi.R) contains the script that has been used to create a phyloseq object, to plot the alpha diversity measures and to run the lefse algorithm.
## Bacteria
The data corresponds to bacteria found on both biopsies and stools of the different types of patients metioned above.
The folder named as plots, has different plots that come as the result of running the R codes.<br />
Content of the files available:
- [bacteria16s_exploratory.R](microbiome_analysis_R/bacteria_biopsies_stools_16S/bacteria16s_exploratory.R) contains the script that has been used to plot the stacked bar plots of the 15 top Genus of bacteria, one for biopsies and another for stools. It also contains the script for running both ANOSIM and PERMANOVA tests on our data, checking the differences between stools and biopsies and within each of those groups, checking the differences between the different types of patients one with each other.
- [bacteria16s_compare.R](microbiome_analysis_R/bacteria_biopsies_stools_16S/bacteria16s_exploratory.R) contains the script that has been used to create a phyloseq object, to plot the alpha diversity meassures of both stools and biopsies and to run the lefse algorithm.
