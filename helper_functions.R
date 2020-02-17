# The functions in this file were generated and tested in analysis.R
# That file can be run to learn how each function is used in the analysis


# Misc Functions ----------------------------------------------------------
h1=function(x){return(strsplit(x,".txt")[[1]])} #removes .txt from filename
h2=function(x,y){
  temp=rep(0,length(x))
  for (z in c(1:length(x))){
    temp[z]=paste(paste("P1_Cluster:",x[z],sep=""),paste("P2_Cluster:",y[z],sep=""),sep=" & ")
    #temp[z]=paste("P2_C:",y[z],sep="")
  }
  return(temp)
} 
h3=function(x,y){
  temp=rep(0,length(x))
  for (z in c(1:length(x))){
    temp[z]=paste(paste("P1_C:",x[z],sep=""),paste("P2_C:",y[z],sep=""),sep=" & ")
    #temp[z]=paste("P2_Cluster:",y[z],sep="")
  }
  return(temp)
} 
h4=function(x){paste("P2 Cluster",x,sep=" ")}
h5=function(x){paste("P1 Cluster",x,sep=" ")}
h6=function(x,y){
  temp=rep(0,length(x))
  for (z in c(1:length(x))){
    temp[z]=paste(paste("",x[z],sep=""),paste("",y[z],sep=""),sep=" vs ")
    #temp[z]=paste("P2_Cluster:",y[z],sep="")
  }
  return(temp)
} 
h10=function(x){return(paste(strsplit(x," ")[[1]],collapse="_"))}


# Filtering Data ----------------------------------------------------------

#### Merges 1 loaded tsne_table file and clincal data
merge_clinical=function(data,clinicalD){
  return(merge(data,clinicalD,by.x="Sample",by.y="samples",all.x=TRUE))
} 

#### Merges 2 TSNE tables and creates labels for survivial curves 
sequential_filter=function(data_1,data_2,Cluster,reverse=FALSE,Pathway_Names,Cluster2){
  if (reverse==TRUE){
    data1=data_2
    data2=data_1
    Pathway_Names=c(Pathway_Names[2],Pathway_Names[1])
  } else {
    data1=data_1
    data2=data_2
  }
  
  if (colnames(data1)[1]=="sample" || colnames(data2)=="sample"){
    colnames(data1)[1]="Sample"
    colnames(data2)[1]="Sample"
  }
  
  temp=data2[,c("Sample","Cluster")]
  temp2=merge(data1,temp,by.x="Sample",by.y="Sample")
  Merged=h2(temp2$Cluster.x,temp2$Cluster.y)
  Merged2=h6(temp2$Cluster.x,temp2$Cluster.y)
  temp3=cbind(temp2,Merged)
  temp4=cbind(temp3,Merged2)
  
  
  if("all"%in%Cluster){final=temp4
  } else {final=temp4[temp4$Cluster.x%in%Cluster,]}
  
  if("all"%in%Cluster2){final=final
  } else {final=final[final$Cluster.y%in%Cluster2,]}
  
  return(list(final,Pathway_Names))
  
} 
sequential_filter_dendro=function(data1,dendro_groups,Cluster,Group2){
  if (colnames(data1)[1]=="sample"){ colnames(data1)[1]="Sample"}
  
  data_1=data1
  data_1[,1]=sapply(data_1[,1], function(x){return(substr(as.vector(x),1,15)[[1]])})
  temp2=merge(data_1,dendro_groups,by.x="Sample",by.y="Sample")
  Merged=h2(as.vector(temp2$Cluster),as.vector(temp2$Group))
  Merged2=h6(as.vector(temp2$Cluster),as.vector(temp2$Group))
  temp3=cbind(temp2,Merged)
  temp4=cbind(temp3,Merged2)
  
  if("all"%in%Cluster){final=temp4
  } else {final=temp4[temp4$Cluster%in%Cluster,]}
  
  if("all"%in%Group2){final=final
  } else {final=final[final$Group%in%Group2,]}
  
  Pathway_Names=c("NA","NA")
  
  return(list(final,Pathway_Names))
}

#### Filter patients using phenotype input
choose_pheno_3D=function(table,phenotypes,thresholds,comparisons,clinical_data_types){
  print("Filtering Patients - ")
  print(phenotypes)
  print("Thresholds:")
  print(thresholds)
  print("Comparisons:")
  print(comparisons)
  filter_index=c(1:nrow(table))
  table=as.matrix(table)
  i=1
  for (pheno in phenotypes){
    compar=comparisons[[i]]
    th=thresholds[[i]]
    column=as.vector(table[,pheno])
    col_type=as.vector(unlist(unname(clinical_data_types[pheno])))
    print("A")
    
    if(col_type=="integer"){column=as.numeric(column)}
    print("B")
    if(compar=="<"){filter_index=intersect(which(column<th),filter_index)}
    else if (compar=="none"){filter_index}
    else if(compar=="="){filter_index=intersect(which(column == th),filter_index)}
    else if(compar=="<="){filter_index=intersect(which(column <= th),filter_index)}
    else if(compar==">"){filter_index=intersect(which(column > th),filter_index)}
    else if(compar==">="){filter_index=intersect(which(column >= th),filter_index)}
    else if(compar=="include"){filter_index=intersect(which(column %in% th),filter_index)}
    else if(compar=="exclude"){
      tmp=intersect(which(column %in% th),filter_index)
      filter_index=setdiff(filter_index,tmp)}
    else if(compar=="between"){
      th1=th[1]
      th2=th[2]
      if(!is.numeric(th) || !is.numeric(th2)){
        stop(paste("Threshold",as.character(i),"is invalid",sep=" "))}
      filter_index=intersect(union(which(column>=th1),which(column<=th2)),filter_index)
    }
    i=i+1
  }
  return(filter_index)
}

#use to control phenotype selection (ex. age multiply by 365 to match)
modify_filter_options=function(phenotypes,thresholds){
  if ("Age" %in% phenotypes){
    i=which(phenotypes == "Age")
    thresholds[[i]]=thresholds[[i]]*365
  }
  return(thresholds)
}


# Creating Survival Curves - Fits (survival) ----------------------------------------
#### Creates survival "FIT" for sequential analysis
my_sfit=function(data){
  
  # if ("all" %in% Cluster){colnames(data)[which(colnames(data)=="Merged")]="C" #Changes labels used for curves
  # } else {colnames(data)[which(colnames(data)=="Merged2")]="C"}
  colnames(data)[which(colnames(data)=="Merged2")]="C"
  fit<-surv_fit(Surv(OS.time,OS)~ C, data=data)
  return(list(fit,data))
} 

#### Creates original survival "FIT" (NO SEQUENTIAL)
my_sfit_P1=function(data){
  fit<-surv_fit(Surv(OS.time,OS)~ Cluster, data=data)
  return(list(fit,data))
} 

#### returns the pval from a generated "FIT"
my_spval=function(fit,data){return(round(surv_pvalue(fit,data=data)$pval,3))} 


# Plots -------------------------------------------------------------------
#### Creates survivial plot for 1 pathway
my_splot_create1=function(fit,data,Pathway_Names,Cluster,pal){
  if("all" %in% Cluster){ C=Cluster
  } else {C=paste(Cluster,collapse =",")}
  pal=pal[1:length(unique(data$C))]
  #title=paste(c("Pathway1:",Pathway_Names[1],"Pathway2:",Pathway_Names[2],"Cluster:",C),collapse =" ")
  return(ggsurvplot(fit, data=data,pval=TRUE,title='',font.tile=1,font.x=8,font.y=8,palette = pal))
} 

#### Creates survivial plot for Sequential Analysis
my_splot_P1_create1=function(fit,data,Pathway,pal){
  pal=pal[1:length(unique(data$Cluster))]
  return(ggsurvplot(fit, data=data,pval=TRUE,title='',font.tile=1,font.x=8,font.y=8,palette = pal))
}

#Uses survivial plot objects in splot and arranges them into grid
my_splot_generate=function(splots,Cancer,pdf=FALSE,filename){
  nplot=length(splots)
  nrow_plot=ceiling(nplot/2)
  if(nplot==1){
    ncol_plot=1
  } else {ncol_plot=2}
  res <- arrange_ggsurvplots(splots, print = TRUE,ncol = ncol_plot, nrow = nrow_plot,title='')
  if(pdf==TRUE){ggsave(filename, res,width=10,height=10)}
  return(res)
} 

# Tables ------------------------------------------------------------------

#Counts number of patients represented in each curve of created survival curve 
count_table=function(table){
  temp=table[,c("Sample","Cluster.x","Cluster.y")]
  u_c.x=sort(unique(temp$Cluster.x))#unique clusters pathway 1
  u_c.y=sort(unique(temp$Cluster.y))#unique cluster pathway 2
  out_table=matrix(-1,length(u_c.x),length(u_c.y))
  row_count=1
  for (i in u_c.x){
    col_count=1
    temp2=temp[which(temp$Cluster.x==i),]
    for (j in u_c.y){
      out_table[row_count,col_count]=length(which(temp2$Cluster.y==j))
      col_count = col_count+1
    }
    row_count=row_count+1
  }
  rownames(out_table)=u_c.x
  colnames(out_table)=u_c.y
  return(out_table)
}

#Counts number of patients in dendro groups
count_table_dendro=function(table){
  temp=table[,c("Sample","Cluster","Group")]
  u_c.x=sort(unique(temp$Cluster))#unique clusters pathway 1
  u_c.y=sort(unique(temp$Group))#unique cluster pathway 2
  out_table=matrix(-1,length(u_c.x),length(u_c.y))
  row_count=1
  for (i in u_c.x){
    col_count=1
    temp2=temp[which(temp$Cluster==i),]
    for (j in u_c.y){
      out_table[row_count,col_count]=length(which(temp2$Group==j))
      col_count = col_count+1
    }
    row_count=row_count+1
  }
  rownames(out_table)=u_c.x
  colnames(out_table)=u_c.y
  return(out_table)
}

# change_pheno_names=function(phenotypes,type){
#     phenotypes[which(phenotypes == "age_at_diagnosis.diagnoses")]=="Age"
#     phenotypes[which(phenotypes == "name.tissue_souce_site")]=="Institution"
# }


# tSNE/UMAP Custom Analysis -----------------------------------------------
# Filtering Data ----------------------------------------------------------
merge_metadata=function()

# Helper_Functions_JM -----------------------------------------------------
# takes a vector and checks for duplicates, INCLUDING the first instance
all_dup <- function (value){
  duplicated(value) | duplicated(value, fromLast = TRUE)
} 




#collapses contingency table to be a 2x2  yes/no contingenccy for one category specified by row and column
contingency_collapse <- function(table,r,c){
  out_mat <- matrix(NA,2,2)
  out_mat[1,1]  <- table[r,c]
  out_mat[1,2]  <- sum(table[r,])-table[r,c]
  out_mat[2,1]  <- sum(table[,c])-table[r,c]
  out_mat[2,2]  <- sum(table) -(out_mat[1,1] + out_mat[1,2] +out_mat[2,1])
  return(out_mat)
}


#exponentiates base two and decrements by 1
exp_dec <- function(x) {
  (2^x)-1
}


load_gen_data_subset <- function(data_dir,use_tpm){
  if(!exists("expression_data")){
    print("data not loaded")
    if(use_tpm){ 
      expression_data <<- readRDS(paste0(data_dir,"data2/merged/joined_tpm_mat_genefilter.rds"))
      tpm_loaded <<- TRUE
    } else {
      expression_data <<- readRDS(paste0(data_dir,"data2/merged/joined_fpkm_mat_genefilter.rds"))
      tpm_loaded <<- FALSE 
    }
    
  }  else  {   
    print("data are loaded")
    if(use_tpm){
      if(tpm_loaded){
      }else{
        expression_data <<- readRDS(paste0(data_dir,"data2/merged/joined_tpm_mat_genefilter.rds"))
        tpm_loaded <<- TRUE
      }
    } else {
      if(tpm_loaded){
        expression_data <<- readRDS(paste0(data_dir,"data2/merged/joined_fpkm_mat_genefilter.rds"))
        tpm_loaded <<- FALSE
      }else{
      }
      
    }  
  }
  return(expression_data)
  
}

#Loads genetic  data we want only if it is not already loaded
#is a function but used like a script
load_gen_data <- function(){
  if(!exists("expression_data")){
    print("data not loaded")
    
    if(use_tpm){ 
      expression_data <<- readRDS(file.path(base_dir,data_dir,tpm_file))
      tpm_loaded <<- TRUE  #asignment operator traces back to first parent environment with 'tpm_loaded' until base environment
    } else {
      expression_data <<- readRDS(file.path(base_dir,data_dir,fpkm_file))
      tpm_loaded <<- FALSE 
    } 
    
    
  } else {   
    print("data are loaded")
    if(use_tpm){
      if(tpm_loaded){
      }else{
        expression_data <<- readRDS(file.path(base_dir,data_dir,tpm_file))
        tpm_loaded <<- TRUE
      }
    } else {
      if(tpm_loaded){
        expression_data <<- readRDS(file.path(base_dir,data_dir,fpkm_file))
        tpm_loaded <<- FALSE
      }else{
      }
      
    }  
  }
  return(expression_data)
}

#gets the ggplot color hues but excludes the inital red hue, useful for using it for noise points
ggcolor_hue_minus <- function(n) {
  hues = seq(15, 375, length = n + 1)
  full <- hcl(h = hues, l = 65, c = 100)[1:n]
  minus_one <-full[2:numel(full)]
  return(minus_one)
}



# Main_Functions_JM -------------------------------------------------------

load_metadata <- function(data_dir)  { #dc_in stands for data container
  #loading the metadata about each sample
  all_metadata <-  readRDS(paste(data_dir,"data2/merged/joined_meta.rds",sep="")) #recall that a single primary can have multiple mets
  surv_data <- readRDS(paste(data_dir,"data2/merged/gdc_survival.rds",sep=""))
  surv_data <- surv_data[,c("sample","OS","OS.time")]
  all_metadata <- inner_join(all_metadata,surv_data,by = "sample")
  dc_out <- list()
  
  
  dc_out$all_metadata <- all_metadata
  return(dc_out)
}

filter_metadata <- function(dc_in,disease_filt,sample_type_filt){
  filter_dataset=TRUE
  filter_disease=TRUE
  filter_sample_type=TRUE
  dataset_filt <- "TCGA"  
  #sample_type_filt <- "Primary Tumor"
  print("TESTB")
  print(dim(dc_in$all_metadata))
  #sample_type_filt <- "Primary Blood Derived Cancer - Peripheral Blood"
  metadata_filtered <- dplyr::filter(dc_in$all_metadata,!is.na("disease"))
  if(filter_dataset){ metadata_filtered <- dplyr::filter(metadata_filtered,dataset %in% dataset_filt)  }
  if(filter_disease){  metadata_filtered <- dplyr::filter(metadata_filtered,disease %in% disease_filt)}
  if(filter_sample_type){  metadata_filtered <- dplyr::filter(metadata_filtered,sample_type %in% sample_type_filt)}
  
  dc_out <- dc_in
  
  dc_out$metadata_filtered <- metadata_filtered
  return(dc_out)
}

load_filter_rename_gen_data2 <- function(dc_in,data_dir,use_tpm,use_pathway,custom_control){
  expression_data <- load_gen_data_subset(data_dir,use_tpm)
  tx=read.csv(paste(data_dir,"data2/merged/merged_pathway_lists.txt",sep=""),sep="\t")
  #tx=read.csv(list.files(file.path(base_dir,data_dir_1,pathway_lists_dir),full.names = TRUE,pattern = ".txt"),sep="\t")
  
  if(custom_control==TRUE){
    tx1=use_pathway
    print(head(tx1))
    index=which(colnames(expression_data) %in% tx1)
    wanted_fpkm <- expression_data[,index]
  } else {
    tx1=as.vector(tx[which(tx$Pathway==use_pathway),"Gene"])
    print(head(tx1))
    wanted_fpkm <- expression_data[,tx1]
  }
  
  print(dim(wanted_fpkm))
  
  print(head(dc_in$metadata_filtered$sample))
  
  wanted_fpkm <- wanted_fpkm[dc_in$metadata_filtered$sample,]
  wanted_fpkm <- wanted_fpkm[, colSums(is.na(wanted_fpkm)) == 0]
  #recall that a single primary can have multiple mets
  #matching up rownames of metadata_filtered and expression
  
  #making the expression data into a tibble that can do it right 
  fpkm_df <- dplyr::as_tibble(wanted_fpkm)
  fpkm_df$sample <- rownames(wanted_fpkm)
  dc_out <- dc_in
  dc_out$fpkm_df <- fpkm_df
  dc_out$expression_data <- expression_data
  return(dc_out)
  
}

merge_gen_meta <- function(dc_in){
  
  metadata_end <- dim(dc_in$metadata)[2]
  meta_fpkm_joined <- merge(dc_in$metadata,dc_in$fpkm_df,by = "sample")
  
  # where are genetic data stored
  expression_indices <- (metadata_end+1):NCOL(meta_fpkm_joined)
  #deleting samples with no gen_data
  has_gen <-  apply(meta_fpkm_joined[,expression_indices], 1, function(x) !any(is.na(x)))
  meta_fpkm_joined <- meta_fpkm_joined[has_gen,]
  
  
  dc_out <- dc_in
  dc_out$expression_indices <- expression_indices
  dc_out$meta_fpkm_joined <- meta_fpkm_joined
  
  return(dc_out)  
}

normalize_data <- function(dc_in){
  standardize <-  TRUE  
  spherize <- TRUE   
  dr_matrix <- data.matrix(dc_in$meta_fpkm_joined[,dc_in$expression_indices])
  
  if(standardize){
    dr_matrix <- sweep(dr_matrix,1,rowSums(dr_matrix),"/")
    dr_matrix[is.na(dr_matrix)] <- 0
  }
  
  
  #centering mean and projecting onto unit hypersphere
  if(spherize){
    centroid <- colMeans(dr_matrix)
    dr_matrix <- sweep(dr_matrix,2,centroid,"-")
    radius <- sqrt(rowSums(dr_matrix^2))
    dr_matrix <- sweep(dr_matrix,1,radius,"/")
  }
  
  
  rownames(dr_matrix) <- dc_in$meta_fpkm_joined$sample
  
  dc_out <- dc_in
  dc_out$dr_matrix <- dr_matrix
  return(dc_out)
  
}

dim_reduce <- function(dc_in,use_umap,custom.config,tsne_perplexity,tsne_max_iter){
  use_umap_wrapper <- FALSE
  
  dr_matrix <- dc_in$dr_matrix
  
  if (use_umap){
    
    umap_method <-    ifelse(use_umap_wrapper,  "umap-learn",  "naive")
    #embedding <- umap(dr_matrix,config = custom.config, method = umap_method,verbose = TRUE)
    embedding <- try(umap(dr_matrix,config = custom.config, method = umap_method,verbose = TRUE),silent=TRUE)
    if (length(embedding)==1){
      if(strsplit(embedding[1],":")[[1]][1]=="Error ") {return("FAILED")}
    } 
    layout_df <- data.frame(embedding$layout)
    
    
    #joining umap layout with meta_tpm joined.  important to not reorder meta_tpm_joined before running umap
    #dr stands for 'dimension reductin'a
    
    
  } else{
    
    use_tsne_wrapper=FALSE
    tsne_dims <- 3 #how many dims should we use
    
    #basic parameters to change for tsne: note that this is a mixture of the parameters from both the tsne and Rtsne package
    #tsne_perplexity <- 10
    #tsne_max_iter <- 500
    
    tsne_eta <- 200
    
    #slightly more advanced parameters to change for tsne
    tsne_stop_lying_iter <- 250L
    tsne_mom_switch_iter <- 250L
    tsne_momentum <- .5
    tsne_final_momentum <- .8
    tsne_exaggeration_factor <- 12
    tsne_theta <- .5    # set to zero for pure t-sne, positive <1 for bhtsne
    tsne_min_cost <- 0   #only tsne package
    
    #things having to do with normalization that we shouldn't set because we have a normalization function
    tsne_pca <- FALSE 
    tsne_pca_center <- TRUE
    tsne_pca_scale <- TRUE
    tsne_normalize <- TRUE
    tsne_initial_dims <- NA #should we reduce the dimension with PCA first? not done here because we want to do this in the normalization function
    tsne_partial_pca <- FALSE
    tsne_whiten <- FALSE   #only tsne package.  should be done in normalization function. 
    tsne_Y_init <- NULL
    
    #things we shouldn't change for this app
    tsne_check_duplicates <- FALSE
    tsne_verbose <- TRUE
    tsne_is_distance <- FALSE  #can pre-compute a distance matrix
    tsne_index <- NA            # for nearest neighbors function (unused)
    tsne_distance <- NA         #for  nearest neighbors function (unused)
    tsne_epoch <- 100            #only tsne package; how long an epoch is 
    tsne_epoch_callback <- NULL  #only tsne package; a function to excecute every epoch
    tsne_num_threads <- 1
    
    
    if (use_tsne_wrapper){
      
      #could include the Rtsne_neighbors here
      
      embedding <-   Rtsne(dr_matrix, 
                           dims = tsne_dims,  
                           initial_dims = dim(dr_matrix)[2],
                           perplexity = tsne_perplexity, 
                           theta = tsne_theta, 
                           check_duplicates = tsne_check_duplicates,
                           pca = tsne_pca, 
                           partial_pca = tsne_partial_pca,
                           max_iter = tsne_max_iter,
                           verbose = tsne_verbose, 
                           is_distance = tsne_is_distance,
                           Y_init = tsne_Y_init,
                           pca_center = tsne_pca_center, 
                           pca_scale = tsne_pca_scale,
                           normalize = tsne_normalize,
                           stop_lying_iter = ifelse(is.null(tsne_Y_init), 
                                                    tsne_stop_lying_iter,0L), 
                           mom_switch_iter = ifelse(is.null(tsne_Y_init), tsne_mom_switch_iter, 0L),
                           momentum = tsne_momentum, 
                           final_momentum = tsne_final_momentum, 
                           eta = tsne_eta,
                           exaggeration_factor = tsne_exaggeration_factor, 
                           num_threads = tsne_num_threads)
      rownames(embedding$Y) <- rownames(dr_matrix)
      layout_df <- data.frame(embedding$Y)
      
    } else {
      embedding <- tsne::tsne(dr_matrix, 
                              initial_config = tsne_Y_init, 
                              k = tsne_dims, 
                              initial_dims = dim(dr_matrix)[2], 
                              perplexity = tsne_perplexity,
                              max_iter = tsne_max_iter, 
                              min_cost = tsne_min_cost, 
                              epoch_callback = tsne_epoch_callback, 
                              whiten = tsne_whiten,
                              epoch = tsne_epoch)
      rownames(embedding) <- rownames(dr_matrix)
      layout_df <- data.frame(embedding)
      
      
    }
    
    
    
    
    
    
  }
  
  joined_dr <- bind_cols(dc_in$meta_fpkm_joined,layout_df) 
  dc_out <- dc_in
  dc_out$embedding <- embedding
  dc_out$joined_dr <- joined_dr
  dc_out$layout_df <- layout_df
  return(dc_out)
}

cluster <- function(dc_in){ 
  jitter_cluster <- FALSE  #should we add noise to the data before clustering?
  jitter_plot <- FALSE     #should we add the clustering noise to the output plot?
  jitter_factor <- NULL   #parameter for jitter function
  jitter_amount <- .005    #parameter for jitter function
  hdbscan_min_pts <- 20    
  
  to_clust <- cbind(dc_in$joined_dr$X1,dc_in$joined_dr$X2,dc_in$joined_dr$X3)
  colnames(to_clust) <- c("X1","X2","X3")
  if(jitter_cluster){
    to_clust %<>% jitter(amount = jitter_amount,factor = jitter_factor)
  } else{ }
  
  dc_in$joined_dr$Cluster <- factor(hdbscan(to_clust,minPts = hdbscan_min_pts)$cluster) 
  
  
  dc_out <- dc_in
  return(dc_out)
}


