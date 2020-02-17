library(rsconnect)

#setwd("~/Documents/Programming/Prochownik/SeqRShiny/app_attempt2")
setwd("~/Documents/Programming/Prochownik/SeqRShiny/final_app_workspace/Survival_analysis_tsne_umap_TCGA")
source("main.R")
source("helper_functions.R")

############################ 

# dc_in=load_metadata(data_dir)
# 
# #CONTROLS
# use_tpm <- FALSE
# disease_filt <- "LAML"
# #sample_type_filt <- "Primary Tumor"
# 
# sample_type_filt <- "Primary Blood Derived Cancer - Peripheral Blood"
# 
# use_pathway="cell_cycle"
# 
# dc_in=filter_metadata(dc_in,disease_filt,sample_type_filt)
# 
# index=c(1:nrow(dc_in$metadata_filtered))
# custom_control=FALSE
# 
# dc_in$filtered_metadata_old=dc_in$metadata_filtered
# dc_in$metadata_filtered=dc_in$metadata_filtered[index,]
# 
# dc_in=load_filter_rename_gen_data2(dc_in,data_dir,use_tpm,use_pathway,custom_control)
# 
# 
# dc_in=merge_gen_meta(dc_in)
# dc_in=normalize_data(dc_in)
# dc_in=dim_reduce(dc_in)
# dc_in=cluster(dc_in)
# dc_in=plot(dc_in)

############################

#rsconnect::setAccountInfo(name='chpupsom19',token='04D9E130FCBEFD71BBF550C31D3CDF67',secret='6Y2jq88rh3tiUNXQ7yUMP6/H34N5Y56M4VYfUBYx')

shiny::runApp()
#deployApp(appTitle = "Survival_Analysis_tsne_umap_TCGA") 
