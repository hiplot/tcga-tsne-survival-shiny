#Main Control Panel
library(shiny)
library(survminer)
library(survival)
library(plotly)
library(ComplexHeatmap)
library(tsne)
library(Rtsne)
library(plotly)
library(shinydashboard)
library(dashboardthemes)
library(dplyr)
library(umap)
library(dbscan)
#library(org.Hs.eg.db)
#library(biomaRt)
#library(rsconnect)

source("helper_functions.R")


main_dir=getwd()
data_dir=paste(main_dir,"data/",sep="/")
clinical_data_file=paste(data_dir,"all_patients_clinical_data.txt",sep="")
clinical_data_file1=paste(data_dir,"all_patients_clinical_data2.txt",sep="")
tsne_table_dir=paste(data_dir,"tsne_tables",sep="")
cancer_mapping_file=paste(data_dir,"cancer_mapping.txt",sep="")
dendro_group_dir=paste(data_dir,"dendro_groups",sep="")
dir1=paste(data_dir,"dendro_annotations",sep="")
dir2=paste(data_dir,"heatmap_base",sep="")
dir3=paste(data_dir,"dtsne_annotations",sep="")
pal_file=paste(data_dir,"aeb_colors.rds",sep="")





