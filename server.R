source("main.R")

# shinyServer -------------------------------------------------------------

server <- function(input, output) {
  
  # Load Data and Modify Input Tables -------------------------------------------------------------------
  
  
  #################### pre-run t-SNE Groups Analysis 1 ####################
  #### Load Data
  load_clinical_data<- reactive({ #Loading Clinical Data (without phenotype variable type)
    load_clinical_data=read.csv(clinical_data_file,sep="\t")})
  
  load_clinical_data1<- reactive({ #Loading Clinical Data (with phenotype variable type)
    load_clinical_data=read.csv(clinical_data_file1,sep="\t")})
  
  load_tsne_table1<-reactive ({ # Load TSNE Table 1
    req(input$pathway1)
    clinicalD=load_clinical_data()
    cancer_name=paste(strsplit(input$cancer," ")[[1]],collapse="_")
    data1=read.csv(paste(c(tsne_table_dir,cancer_name,paste(input$pathway1,".txt",sep="")),collapse = "/"))
    load_tsne_table1=merge_clinical(data1,clinicalD)
    
  })
  
  load_tsne_table2<-reactive ({ # Load TSNE Table 2
    req(input$pathway2)
    clinicalD=load_clinical_data()
    cancer_name=paste(strsplit(input$cancer," ")[[1]],collapse="_")
    data1=read.csv(paste(c(tsne_table_dir,cancer_name,paste(input$pathway2,".txt",sep="")),collapse = "/"))
    load_tsne_table2=merge_clinical(data1,clinicalD)
  })
  
  #merges tsnetables
  load_tsne_tables<-reactive ({ # Merge 
    req(input$clusters)
    req(input$clusters2)
    clinicalD=load_clinical_data()
    table1=load_tsne_table1()
    table2=load_tsne_table2()
    
    index=filter_patients()
    data1=table1[index,]
    data2=table2[index,]
    
    Pathway_Names=c(input$pathway1,input$pathway2)
    Cluster=input$clusters
    Cluster2=input$clusters2
    if ("all" %in% Cluster){Clusters1=sort(unique(data1$Cluster))
    } else{Clusters1=as.integer(Cluster)}
    load_tsne_tables=sequential_filter(data1,data2,Clusters1,reverse = FALSE,Pathway_Names,Cluster2)
  })
  
  display_data<- reactive ({ #Not Used but can show Merged Data used in Analysis
    req(input$pathway2)
    data=load_tsne_tables()[[1]]
    n=ncol(data)
    display_data=data[,-c(n-1,n)]
  })
  
  load_mapping<-reactive({ #Mapping Data for heirarchical clustering link
    load_mapping=read.csv(cancer_mapping_file,sep="\t",header=TRUE)})
  
  load_dendrogram_groups<-reactive ({ 
    req(input$cancer)
    mapping=load_mapping()
    cancer=input$cancer
    cancer_name=paste(strsplit(input$cancer," ")[[1]],collapse="_")
    cancer_index=which(mapping$Full_Name==cancer_name)
    symbol=as.vector(mapping[cancer_index,"Symbol"])
    dendg=readRDS(paste(c(dendro_group_dir,"/",symbol,"_dendro_groups.RDS"),collapse=""))
    dendro_groups=cbind(names(dendg),unname(dendg))
    colnames(dendro_groups)=c("Sample","Group")
    load_dendogram_groups=dendro_groups
    
  })
  
  load_tsne_tables_dendro <-reactive ({ #merge t-sne data for pathway 1 with dendrogram grouping for survival curves
    req(input$clusters_heatmap)
    req(input$group2)
    data1=load_tsne_table1()
    dendro_groups=load_dendrogram_groups()
    Cluster=input$clusters_heatmap
    Group2=input$group2
    if ("all" %in% Cluster){Clusters1=sort(unique(data1$Cluster))
    } else{Clusters1=as.integer(Cluster)}
    load_tsne_tables_dendro=sequential_filter_dendro(data1,dendro_groups,Cluster,Group2)
  })
  
  #### Count tables from input data 
  counttable<- reactive({ #Generated table of counts under survivial plots in Tab 2 
    req(input$pathway2)
    req(input$clusters)
    data=load_tsne_tables()[[1]]
    counttable1=count_table(data)
    rownames(counttable1)=sapply(rownames(counttable1),function(x){paste("Pathway 1: Cluster",as.character(x))})
    colnames(counttable1)=sapply(colnames(counttable1),function(x){paste("Pathway 2: Cluster",as.character(x))})
    counttable=counttable1
  })
  
  #### Count Table for Dendro Groupos
  counttable_dendro<- reactive({ #Generated table of counts under survivial plots in Tab 2 
    req(input$group2)
    req(input$clusters_heatmap)
    table=load_tsne_tables_dendro()[[1]]
    index=filter_patients()
    data=table[index,]
    counttable_dendro1=count_table_dendro(data)
    rownames(counttable_dendro1)=sapply(rownames(counttable_dendro1),function(x){paste("Pathway 1: Cluster",as.character(x))})
    colnames(counttable_dendro1)=sapply(colnames(counttable_dendro1),function(x){paste("Pathway 2: Cluster",as.character(x))})
    counttable_dendro = counttable_dendro1
  })
  
  #### Count Table for Pathway 1
  count_table_scatter1<-reactive({
    req(input$cancer)
    req(input$pathway1)
    data <- load_tsne_table1()
    index=filter_patients()
    validate(need(index>0,"No patients with selected filter options"))
    if(!is.null(data) && !is.null(index) ){
      data=data[index,]
      
      table=table(data$Cluster)
      m=matrix(0,1,length(table))
      m[1,]=unname(table)
      rownames(m)=c("Number of Patients")
      colnames(m)=sapply(c(1:length(table)),function(x){return(paste("Cluster",x,sep=" "))})
      count_table_scatter1=m
    } else {count_table_scatter1=NULL}
    
  })
  
  #### Count Table for Pathway 2
  count_table_scatter2<-reactive({
    req(input$cancer)
    req(input$pathway1)
    req(input$pathway2)
    data <- load_tsne_table2()
    index=filter_patients()
    validate(need(index>0,"No patients with selected filter options"))
    if(!is.null(data) && !is.null(index)){
      data=data[index,]
      table=table(data$Cluster)
      m=matrix(0,1,length(table))
      m[1,]=unname(table)
      rownames(m)=c("Number of Patients")
      colnames(m)=sapply(c(1:length(table)),function(x){return(paste("Cluster",x,sep=" "))})
      count_table_scatter2=m
    } else {count_table_scatter2=NULL}
  })
  
  #Load Color palette for plots 
  load_pal <-reactive ({
    load_pal=readRDS(pal_file)
    load_pal=unname(load_pal)
  })
  
  #################### tSNE/UMAP Analysis 2 ####################
  
  load_dendrogram_groups_custom <-reactive ({ 
    req(input$custom_cancer)
    mapping=load_mapping()
    cancer=input$custom_cancer
    cancer_name=paste(strsplit(input$custom_cancer," ")[[1]],collapse="_")
    cancer_index=which(mapping$Full_Name==cancer_name)
    symbol=as.vector(mapping[cancer_index,"Symbol"])
    dendg=readRDS(paste(c(dendro_group_dir,"/",symbol,"_dendro_groups.RDS"),collapse=""))
    dendro_groups=cbind(names(dendg),unname(dendg))
    colnames(dendro_groups)=c("Sample","Group")
    load_dendogram_groups=dendro_groups
    
  })
  
  #Gene/Pathway List
  load_gene_pathways <- reactive({
    load_gene_pathways=read.csv(paste(data_dir,"data2/merged/merged_pathway_lists.txt",sep=""),sep="\t")
  })
  
  load_phenotype_map_custom <- reactive({
    phenotype_var_type=read.csv(paste0(data_dir,"data2/merged/phenotype_mapping.txt"),sep="\t")
    clinical_data_types=as.vector(phenotype_var_type$var_type)
    names(clinical_data_types)=as.vector(phenotype_var_type$Phenotype)
    return(clinical_data_types)
  })
  
  #Pathway 1
  
  #Step 1
  load_custom_analysis_cancer_data <- reactive({
    req(input$custom_cancer)
    cancer_map=read.csv(paste0(data_dir,"data2/cancer_name_map.txt"),sep="\t")
    disease_filt=as.vector(cancer_map[which(cancer_map$Full == h10(input$custom_cancer)),"Abrv"])
    #if(input$use_tpm==TRUE){use_tpm=TRUE} else {use_tpm=FALSE}
    
    if(disease_filt=="LAML"){sample_type_filt="Primary Blood Derived Cancer - Peripheral Blood"} else {sample_type_filt="Primary Tumor"}
    
    print("LOADING DATA")
    dc_in=load_metadata(data_dir)
    dc_in=filter_metadata(dc_in,disease_filt,sample_type_filt)
    return(dc_in)
  })
  
  #Step 2
  filter_custom_p1 <- eventReactive(input$run_custom,{
    req(input$custom_pathway1)
    req(input$use_tpm)
    dc_in=load_custom_analysis_cancer_data()
    data=dc_in$metadata_filtered
    #input$filter_control=TRUE
    phenotypes=c(input$custom_phenotypes)
    if(is.null(data)){
      filter_index=NULL
    } else {
      print(input$filter_control_custom)
      if(input$filter_control_custom!="No Filter" && !is.null(phenotypes) && !is.null(input$custom_threshold1) && !is.null(input$custom_compare1)){
        num_phenotypes=length(phenotypes)
        if(num_phenotypes==2){
          req(input$threshold2)
          req(input$compare2)
        }
        thresholds=list()
        comparisons=list()
        for (p in c(1:num_phenotypes)){
          if(p==1){
            thresholds[[1]]=input$custom_threshold1
            comparisons[[1]]=input$custom_compare1
          }
          if(p==2){
            print(input$threshold2)
            thresholds[[2]]=input$custom_threshold2
            comparisons[[2]]=input$custom_compare2
          }
        }
        thresholds=modify_filter_options(phenotypes,thresholds) 
        clinical_data_types=load_phenotype_map_custom ()
        filter_index=choose_pheno_3D(data,phenotypes,thresholds,comparisons,clinical_data_types)
      } else {
        print("A")
        filter_index=c(1:nrow(data))
      }
      print(length(filter_index))
      return(filter_index)
    }
  })
  
  #Step 3
  custom_analysis_pathway1 <- eventReactive(input$run_custom,{
    req(input$custom_pathway1)
    req(input$use_tpm)
    index=filter_custom_p1()
    #print(length(index))
    if(input$use_tpm==TRUE){use_tpm=TRUE} else {use_tpm=FALSE}
    use_pathway=input$custom_pathway1
    dc_in=load_custom_analysis_cancer_data()
    
    
    if(use_pathway=="Custom Gene Set"){
      #req(input$multigene_p1)
      custom_control=TRUE
      use_pathway=input$multigene_p1
      use_pathway=use_pathway
    } else { custom_control=FALSE}

    dc_in$filtered_metadata_old=dc_in$metadata_filtered
    dc_in$metadata_filtered=dc_in$metadata_filtered[index,]
    dc_in=load_filter_rename_gen_data2(dc_in,data_dir,use_tpm,use_pathway,custom_control)
    dc_in=merge_gen_meta(dc_in)
    dc_in=normalize_data(dc_in)
    
    
    print("Running UMAP")
    custom.config <- umap.defaults 
    #custom.config$n_epochs <- 100    #how much should we run umap for
    #custom.config$min_dist <- 10^(-66)    #minimum distance for umap
    #custom.config$n_neighbors <- 30      #number of neighbors for umap
    custom.config$n_components <- 3     #shouldn't be changed by user: how many dimensions we are doing umap into
    
    if(input$use_UMAP=="UMAP"){
      use_umap=TRUE
      #print(input$n_neighbors)
      custom.config$n_epochs=as.numeric(input$n_epochs)
      custom.config$min_dist=10^(as.numeric(input$min_dist))
      custom.config$n_neighbors=as.numeric(input$n_neighbors)
      tsne_perplexity=NULL
      tsne_max_iter=NULL
    } else {
        use_umap=FALSE
        tsne_perplexity=as.numeric(input$perplexity)
        tsne_max_iter=as.numeric(input$max_iter)
    }
    
    dc_out=dim_reduce(dc_in,use_umap,custom.config,tsne_perplexity,tsne_max_iter)
    
    
    print("TEST")
    if (length(dc_out) == 1){
      if (dc_out=="FAILED") {return("FAILED")
      } else {
        dc_in=cluster(dc_in)
        print("DONE")
     }
    }
    #temp=dc_in$joined_dr
    dc_in=cluster(dc_out)
    print("DONE")
    
    #temp=dc_in$joined_dr
    print(colnames(dc_in$joined_dr))
    custom_analysis_pathway2=dc_in$joined_dr
    
  })
  
  custom_analysis_pathway2 <- eventReactive(input$run_custom,{
    req(input$custom_pathway2)
    req(input$use_tpm)
    index=filter_custom_p1()
    print(length(index))
    if(input$use_tpm==TRUE){use_tpm=TRUE} else {use_tpm=FALSE}
    use_pathway=input$custom_pathway2
    dc_in=load_custom_analysis_cancer_data()
    
    if(use_pathway=="Custom Gene Set"){
      #req(input$multigene_p1)
      custom_control=TRUE
      use_pathway=input$multigene_p2
      use_pathway=use_pathway
    } else {custom_control=FALSE}

    dc_in$filtered_metadata_old=dc_in$metadata_filtered
    dc_in$metadata_filtered=dc_in$metadata_filtered[index,]
    dc_in=load_filter_rename_gen_data2(dc_in,data_dir,use_tpm,use_pathway,custom_control)
    dc_in=merge_gen_meta(dc_in)
    dc_in=normalize_data(dc_in)
    
    print("Running UMAP")
    custom.config <- umap.defaults 
    #custom.config$n_epochs <- 100    #how much should we run umap for
    #custom.config$min_dist <- 10^(-66)    #minimum distance for umap
    #custom.config$n_neighbors <- 30      #number of neighbors for umap
    custom.config$n_components <- 3     #shouldn't be changed by user: how many dimensions we are doing umap into
    
    if(input$use_UMAP=="UMAP"){
      use_umap=TRUE
      print(input$n_neighbors)
      custom.config$n_epochs=as.numeric(input$n_epochs)
      custom.config$min_dist=10^(as.numeric(input$min_dist))
      custom.config$n_neighbors=as.numeric(input$n_neighbors)
      tsne_perplexity=NULL
      tsne_max_iter=NULL
    } else {
      use_umap=FALSE
      tsne_perplexity=as.numeric(input$perplexity)
      tsne_max_iter=as.numeric(input$max_iter)
    }

    
    dc_out=dim_reduce(dc_in,use_umap,custom.config,tsne_perplexity,tsne_max_iter)
    
    
    print("TEST")
    if (length(dc_out) == 1){
      if (dc_out=="FAILED") {return("FAILED")
      } else {
        dc_in=cluster(dc_in)
        print("DONE")
      }
    }
    #temp=dc_in$joined_dr
    dc_in=cluster(dc_out)
    print("DONE")
    
    #temp=dc_in$joined_dr
    print(colnames(dc_in$joined_dr))
    custom_analysis_pathway2=dc_in$joined_dr
    
  })
  
  #Num Patients Tables for Scatter Plots
  count_table_custom_scatter1<-reactive({
    req(input$custom_cancer)
    req(input$custom_pathway1)
    req(input$use_tpm)
    data=custom_analysis_pathway1()
    validate(
      need(data!="FAILED", "")
    )
    #data1=data[,c(1,c(ncol(data)-3:ncol(data)))]
    data1=data[,c("sample","X1","X2","X3","Cluster")]
    colnames(data1)=c("Sample","X","Y","Z","Cluster")
    table=table(data1$Cluster)
    m=matrix(0,1,length(table))
    m[1,]=unname(table)
    rownames(m)=c("Number of Patients")
    colnames(m)=sapply(c(1:length(table)),function(x){return(paste("Cluster",x,sep=" "))})
    count_table_custom_scatter1=m
  })
  
  count_table_custom_scatter2<-reactive({
    req(input$custom_cancer)
    req(input$custom_pathway2)
    req(input$use_tpm)
    data=custom_analysis_pathway2()
    validate(
      need(data!="FAILED", "")
    )
    #data1=data[,c(1,c(ncol(data)-3:ncol(data)))]
    data1=data[,c("sample","X1","X2","X3","Cluster")]
    colnames(data1)=c("Sample","X","Y","Z","Cluster")
    table=table(data1$Cluster)
    m=matrix(0,1,length(table))
    m[1,]=unname(table)
    rownames(m)=c("Number of Patients")
    colnames(m)=sapply(c(1:length(table)),function(x){return(paste("Cluster",x,sep=" "))})
    count_table_custom_scatter1=m
  })
  
  custom_warning<- eventReactive(input$run_custom,{
      custom_warning="Please wait, analysis may take a few seconds to load "
  })
  
  
  
  # Analysis ----------------------------------------------------------------
  
    #################### pre-run t-SNE Groups Analysis 1 ####################
  
  
  ######## Individual Pathway Survival Curves 
  fit_P1<- reactive ({
    data=load_tsne_table1()
    index=filter_patients()
    table=data[index,]
    fit_P1=my_sfit_P1(table)
  })
  
  fit_P2<- reactive ({
    data=load_tsne_table2()
    index=filter_patients()
    table=data[index,]
    fit_P2=my_sfit_P1(table)
  })
  
  ######## Sequential Analysis Survivial Fit
  fit<- reactive ({
    table=load_tsne_tables()[[1]]
    Clusters=input$clusters
    # if ("all" %in% Cluster){Clusters1="all"
    # } else{Clusters1=as.integer(Cluster)}
    fit=my_sfit(table)
  })
  
  ######## Dendrogram Analysis Survivial Fit
  fit_dendro<-reactive({
    data=load_tsne_tables_dendro()[[1]]
    index=filter_patients()
    table=data[index,]
    fit_dendro=my_sfit(table)
  })
  
  ######## Choose patients by phenotype
  filter_patients<-reactive({
    req(input$cancer)
    req(input$pathway1)
    data=load_tsne_table1()
    if(is.null(data)){
      filter_index=NULL
    } else {
      if(input$filter_control==TRUE && !is.null(input$phenotype) && !is.null(input$threshold1) && !is.null(input$compare1)){
        phenotypes=c(input$phenotype)
        num_phenotypes=length(phenotypes)
        if(num_phenotypes==2){
          req(input$threshold2)
          req(input$compare2)
        }
        thresholds=list()
        comparisons=list()
        for (p in c(1:num_phenotypes)){
          if(p==1){
            thresholds[[1]]=input$threshold1
            comparisons[[1]]=input$compare1
          }
          if(p==2){
            print(input$threshold2)
            thresholds[[2]]=input$threshold2
            comparisons[[2]]=input$compare2
          }
        }
        thresholds=modify_filter_options(phenotypes,thresholds)
        clinicalD=load_clinical_data1()
        clinical_data_types=as.vector(clinicalD[1,])
        filter_index=choose_pheno_3D(data,phenotypes,thresholds,comparisons,clinical_data_types)
      } else {
        filter_index=c(1:nrow(data))
      }
      print(length(filter_index))
      return(filter_index)
    }
  })
  
  
  
    #################### tSNE/UMAP Analysis 2 ####################
  fit_P1_custom <- reactive ({
    req(input$custom_cancer)
    req(input$custom_pathway1)
    req(input$use_tpm)
    table=custom_analysis_pathway1()
    validate(
      need(table!="FAILED", "")
    )
    fit_P1_custom=my_sfit_P1(table)
  })
  
  fit_P2_custom <- reactive ({
    req(input$custom_cancer)
    req(input$custom_pathway2)
    req(input$use_tpm)
    table=custom_analysis_pathway2()
    validate(
      need(table!="FAILED", "")
    )
    fit_P1_custom=my_sfit_P1(table)
  })
  
  merge_custom_pathways<- reactive({
    req(input$custom_cancer)
    req(input$custom_pathway1)
    req(input$custom_pathway2)
    req(input$use_tpm)
    data1=custom_analysis_pathway1()
    data2=custom_analysis_pathway2()
    Pathway_Names=c(input$custom_pathway1,input$custom_pathway2)
    Cluster=input$clusters_custom
    Cluster2=input$clusters2_custom
    if ("all" %in% Cluster){Clusters1=sort(unique(data1$Cluster))
    } else{Clusters1=as.integer(Cluster)}
    merge_data=sequential_filter(data1,data2,Clusters1,reverse = FALSE,Pathway_Names,Cluster2)
    return(merge_data)
  })
  
  fit_custom <- reactive({
    merge_data=merge_custom_pathways()
    table=merge_data[[1]]
    fit=my_sfit(table)
  })
  
  filter_patients_custom_p1<-reactive({
    req(input$cancer)
    req(input$custom_pathway1)
    req(input$use_tpm)
    data=custom_analysis_pathway1()
    
    #TEST INPUT
    input=list()
    input$filter_control=TRUE
    input$phenotype="OS"
    input$threshold1=0
    input$compare1="="
    if(is.null(data)){
      filter_index=NULL
    } else {
      if(input$filter_control==TRUE && !is.null(input$phenotype) && !is.null(input$threshold1) && !is.null(input$compare1)){
        phenotypes=c(input$phenotype)
        num_phenotypes=length(phenotypes)
        if(num_phenotypes==2){
          req(input$threshold2)
          req(input$compare2)
        }
        thresholds=list()
        comparisons=list()
        for (p in c(1:num_phenotypes)){
          if(p==1){
            thresholds[[1]]=input$threshold1
            comparisons[[1]]=input$compare1
          }
          if(p==2){
            print(input$threshold2)
            thresholds[[2]]=input$threshold2
            comparisons[[2]]=input$compare2
          }
        }
        thresholds=modify_filter_options(phenotypes,thresholds) 
        phenotype_var_type=read.csv(paste0(data_dir,"data2/merged/phenotype_mapping.txt"),sep="\t")
        clinical_data_types=as.vector(phenotype_var_type$var_type)
        names(clinical_data_types)=as.vector(phenotype_var_type$Phenotype)
        filter_index=choose_pheno_3D(data,phenotypes,thresholds,comparisons,clinical_data_types)
      } else {
        filter_index=c(1:nrow(data))
      }
      print(length(filter_index))
      return(filter_index)
    }
  })
  
  counttable_custom<- reactive({ #Generated table of counts under survivial plots in Tab 2 
    req(input$custom_pathway2)
    req(input$clusters_custom)
    req(input$clusters2_custom)
    data=merge_custom_pathways()[[1]]
    counttable1=count_table(data)
    rownames(counttable1)=sapply(rownames(counttable1),function(x){paste("Pathway 1: Cluster",as.character(x))})
    colnames(counttable1)=sapply(colnames(counttable1),function(x){paste("Pathway 2: Cluster",as.character(x))})
    counttable_custom=counttable1
  })
  
  load_tsne_tables_dendro_custom <-reactive ({ #merge t-sne data for pathway 1 with dendrogram grouping for survival curves
    req(input$custom_clusters_heatmap)
    req(input$custom_group2)
    data1=custom_analysis_pathway1()
    dendro_groups=load_dendrogram_groups_custom()
    Cluster=input$custom_clusters_heatmap
    Group2=input$custom_group2
    if ("all" %in% Cluster){Clusters1=sort(unique(data1$Cluster))
    } else{Clusters1=as.integer(Cluster)}
    load_tsne_tables_dendro=sequential_filter_dendro(data1,dendro_groups,Cluster,Group2)
  })
  
  fit_custom_dendro<-reactive({
    req(input$custom_pathway1)
    table=load_tsne_tables_dendro_custom()[[1]]
    fit_dendro=my_sfit(table)
  })
  
  custom_counttable_dendro<- reactive({ #Generated table of counts under survivial plots in Tab 2 
    req(input$custom_group2)
    req(input$custom_clusters_heatmap)
    data=load_tsne_tables_dendro_custom()[[1]]
    counttable_dendro1=count_table_dendro(data)
    rownames(counttable_dendro1)=sapply(rownames(counttable_dendro1),function(x){paste("Pathway 1: Cluster",as.character(x))})
    colnames(counttable_dendro1)=sapply(colnames(counttable_dendro1),function(x){paste("Pathway 2: Cluster",as.character(x))})
    counttable_dendro = counttable_dendro1
  })
  
  
  
  
  
  # Plots -------------------------------------------------------------------
  
  #################### pre-run t-SNE Groups Analysis ####################
  
  ######## Survival Plots 
  survivalplot_P1_gen<-reactive({
    validate(
      need(input$cancer != "", "")
    )
    validate(
      need(input$pathway1 != "", "")
    )
    req(input$cancer)
    req(input$pathway1)
    f=fit_P1()
    pal=load_pal()
    survivalplot_P1_gen1=my_splot_P1_create1(f[[1]],f[[2]],input$pathway1,pal)
    splots=list(survivalplot_P1_gen1)
    survivalplot_P1_gen=my_splot_generate(splots,input$cancer,pdf=FALSE,filename='test_P1.pdf')
    
  })
  
  survivalplot_dendo_groups_gen<-reactive({
    validate(
      need(input$cancer != "", "Please select a Cancer")
    )
    validate(
      need(input$cancer != "Acute Myeloid Leukemia Marrow", "Acute Myeloid Leukemia marrow does not have heiarchical clustering data for analysis")
    )
    validate(
      need(input$pathway1 != "", "Please select 1st Pathway")
    )
    validate(
      need(input$clusters_heatmap != "", "Please select clusters in 1st pathway to analyze")
    )
    validate(
      need(input$group2!= "", "Please select groups from dendrogram in heirarchical clustering to analyze")
    )
    
    req(input$cancer)
    req(input$pathway1)
    req(input$clusters_heatmap)
    req(input$group2)
    data=load_tsne_tables_dendro()[[1]]
    index=filter_patients()
    table=data[index,]
    
    f=fit_dendro()
    Cluster=input$clusters_heatmap
    if ("all" %in% Cluster){Clusters1="all"
    } else{Clusters1=as.integer(Cluster)}
    Group2=input$group2
    if ("all" %in% Group2){Group2a="all"
    } else{Group2a=as.integer(Group2)}
    pal=load_pal()
    survivalplot_gen1=my_splot_create1(f[[1]],f[[2]],table[[2]],Clusters1,pal)
    splots=list(survivalplot_gen1)
    survivalplot_gen=my_splot_generate(splots,input$cancer,FALSE,filename=NULL)
  })
  
  survivalplot_P2_gen<-reactive({
    req(input$cancer)
    req(input$pathway1)
    req(input$pathway2)
    f=fit_P2()
    pal=load_pal()
    survivalplot_P2_gen1=my_splot_P1_create1(f[[1]],f[[2]],input$pathway2,pal)
    splots=list(survivalplot_P2_gen1)
    survivalplot_P2_gen=my_splot_generate(splots,input$cancer,pdf=FALSE,filename='test_P2.pdf')
    
  })
  
  survivalplot_gen<-reactive({
    validate(
      need(input$cancer != "", "Please select a Cancer")
    )
    validate(
      need(input$pathway1 != "", "Please select 1st Pathway")
    )
    validate(
      need(input$pathway2 != "", "Please select 2nd Pathway")
    )
    validate(
      need(input$clusters != "", "Please select clusters in 1st pathway to analyze")
    )
    validate(
      need(input$clusters2 != "", "Please select clusters in 2nd pathway to analyze")
    )
    req(input$cancer)
    req(input$pathway1)
    req(input$pathway2)
    req(input$clusters)
    table=load_tsne_tables()
    f=fit()
    Cluster=input$clusters
    if ("all" %in% Cluster){Clusters1="all"
    } else{Clusters1=as.integer(Cluster)}
    Cluster2=input$clusters2
    if ("all" %in% Cluster2){Clusters2="all"
    } else{Clusters2=as.integer(Cluster2)}
    pal=load_pal()
    survivalplot_gen1=my_splot_create1(f[[1]],f[[2]],table[[2]],Clusters1,pal)
    splots=list(survivalplot_gen1)
    survivalplot_gen=my_splot_generate(splots,input$cancer,FALSE,filename=NULL)
  })
  
  ######## 3D Scatter Plots 
  scatterplot_gen2_p1<-reactive({
    req(input$cancer)
    req(input$pathway1)
    
    data <- load_tsne_table1()
    index=filter_patients()
    table=data[index,]
    
    pal=load_pal()
    colors=pal[1:length(unique(table$Cluster))]
    
    
    table$Cluster=as.character(table$Cluster)
    scatterplot_gen2<-plot_ly(table,x=~X,y=~Y,z=~Z,color=~Cluster,type='scatter3d',mode="markers",opacity=.3,colors=colors,marker=list(size=7))

  })
  
  scatterplot_gen2_p2<-reactive({
    req(input$cancer)
    req(input$pathway1)
    
    data <- load_tsne_table2()
    index=filter_patients()
    table=data[index,]
    
    pal=load_pal()
    colors=pal[1:length(unique(table$Cluster))]
    
    
    table$Cluster=as.character(table$Cluster)
    scatterplot_gen2<-plot_ly(table,x=~X,y=~Y,z=~Z,type='scatter3d',mode="markers",color=~Cluster,opacity=.3,colors=colors,marker=list(size=7))
  })
  
  ######## Heatmaps
  heatmap_gen<-eventReactive(input$create_heatmap,{
    validate(
      need(input$cancer != "", "Please select a Cancer")
    )
    validate(
      need(input$cancer != "Acute Myeloid Leukemia Marrow", "Acute Myeloid Leukemia marrow does not have heiarchical clustering data for analysis")
    )
    validate(
      need(input$pathway1 != "", "Please select 1st Pathway")
    )
    req(input$cancer)
    req(input$pathway1)
    mapping=load_mapping()
    cancer_name=paste(strsplit(input$cancer," ")[[1]],collapse="_")
    cancer_index=which(mapping$Full_Name==cancer_name)
    cancer=as.vector(mapping[cancer_index,"Symbol"])
    pathway=input$pathway1
    base_dir=data_dir
    hm <- readRDS(file.path(base_dir,"heatmap_base",paste0(cancer,"_heatmap.RDS")))
    base_ano <- readRDS(file.path(base_dir,"dendro_annotations",paste0(cancer,"_annotation.RDS")))
    tsne_ano <- readRDS(file.path(base_dir,"tsne_annotations",cancer,paste0(pathway,"#tsneAnnotation.RDS")))
    heatmap_gen <- draw(tsne_ano%v% base_ano %v% hm)
  })
  
  #################### tSNE/UMAP Analysis 2 ####################
  custom_scatterplot_gen2_p1<-eventReactive(input$run_custom,{
    req(input$custom_cancer)
    req(input$custom_pathway1)
    req(input$use_tpm)
    data=custom_analysis_pathway1()
    validate(
      need(data!="FAILED", "Phenotypes choosen not sufficent for analysisT")
    )
    data1=data[,c("sample","X1","X2","X3","Cluster")]
    colnames(data1)=c("Sample","X","Y","Z","Cluster")
    table=data1
    pal=load_pal()
    colors=pal[1:length(unique(table$Cluster))]
    #table$Cluster=as.character(table$Cluster)
    #print(table$Cluster)
    plot_ly(table,x=~X,y=~Y,z=~Z,color=~Cluster,type='scatter3d',mode="markers",opacity=.3,colors=colors,marker=list(size=7))
  })
  
  custom_scatterplot_gen2_p2<-eventReactive(input$run_custom,{
    req(input$custom_cancer)
    req(input$custom_pathway2)
    req(input$use_tpm)
    data=custom_analysis_pathway2()
    validate(
      need(data!="FAILED", "Phenotypes choosen not sufficent for analysis")
    )
    data1=data[,c("sample","X1","X2","X3","Cluster")]
    colnames(data1)=c("Sample","X","Y","Z","Cluster")
    table=data1
    pal=load_pal()
    colors=pal[1:length(unique(table$Cluster))]
    #table$Cluster=as.character(table$Cluster)
    #print(table$Cluster)
    plot_ly(table,x=~X,y=~Y,z=~Z,color=~Cluster,type='scatter3d',mode="markers",opacity=.3,colors=colors,marker=list(size=7))
  })
  
  custom_survivalplot_P1_gen<-eventReactive(input$run_custom,{
    req(input$custom_cancer)
    req(input$custom_pathway1)
    req(input$use_tpm)
    f=fit_P1_custom()
    pal=load_pal()
    survivalplot_P1_gen1=my_splot_P1_create1(f[[1]],f[[2]],input$pathway1,pal)
    splots=list(survivalplot_P1_gen1)
    survivalplot_P1_gen=my_splot_generate(splots,input$cancer,pdf=FALSE,filename='test_P1.pdf')
  })
  
  custom_survivalplot_P2_gen<-eventReactive(input$run_custom,{
    req(input$custom_cancer)
    req(input$custom_pathway2)
    req(input$use_tpm)
    f=fit_P2_custom()
    pal=load_pal()
    survivalplot_P2_gen1=my_splot_P1_create1(f[[1]],f[[2]],input$pathway1,pal)
    splots=list(survivalplot_P2_gen1)
    survivalplot_P2_gen=my_splot_generate(splots,input$cancer,pdf=FALSE,filename='test_P1.pdf')
  })
  
  custom_survivalplot_gen<-reactive({
    req(input$custom_cancer)
    req(input$custom_pathway1)
    validate(
      need(input$custom_pathway2 != "", "Please select 2nd Pathway for sequential analysis")
    )
    req(input$custom_pathway2)
    validate(
      need(input$clusters_custom != "", "Please select clusters in 1st pathway to analyze")
    )
    validate(
      need(input$clusters2_custom != "", "Please select clusters in 2nd pathway to analyze")
    )
    req(input$clusters_custom)
    req(input$clusters2_custom)
    table=merge_custom_pathways()
    f=fit_custom()
    Cluster=input$clusters_custom
    if ("all" %in% Cluster){Clusters1="all"
    } else{Clusters1=as.integer(Cluster)}
    Cluster2=input$clusters2_custom
    if ("all" %in% Cluster2){Clusters2="all"
    } else{Clusters2=as.integer(Cluster2)}
    pal=load_pal()
    survivalplot_gen1=my_splot_create1(f[[1]],f[[2]],table[[2]],Clusters1,pal)
    splots=list(survivalplot_gen1)
    survivalplot_gen=my_splot_generate(splots,input$cancer,FALSE,filename=NULL)
  })
  
  custom_survivalplot_dendo_groups_gen<-reactive({
    validate(
      need(input$custom_cancer != "", "Please select a Cancer")
    )
    validate(
      need(input$custom_cancer != "Acute Myeloid Leukemia Marrow", "Acute Myeloid Leukemia marrow does not have heiarchical clustering data for analysis")
    )
    validate(
      need(input$custom_pathway1 != "", "Please select 1st Pathway")
    )
    validate(
      need(input$custom_clusters_heatmap != "", "Please select clusters in 1st pathway to analyze")
    )
    validate(
      need(input$custom_group2!= "", "Please select groups from dendrogram in heirarchical clustering to analyze")
    )
    
    req(input$custom_cancer)
    req(input$custom_pathway1)
    req(input$custom_clusters_heatmap)
    req(input$custom_group2)
    print("TEST2")
    table=load_tsne_tables_dendro_custom()[[1]]
    
    f=fit_custom_dendro()
    Cluster=input$custom_clusters_heatmap
    if ("all" %in% Cluster){Clusters1="all"
    } else{Clusters1=as.integer(Cluster)}
    Group2=input$group2
    if ("all" %in% Group2){Group2a="all"
    } else{Group2a=as.integer(Group2)}
    pal=load_pal()
    survivalplot_gen1=my_splot_create1(f[[1]],f[[2]],table[[2]],Clusters1,pal)
    splots=list(survivalplot_gen1)
    survivalplot_gen=my_splot_generate(splots,input$cancer,FALSE,filename=NULL)
  })
  
  custom_heatmap_gen<-eventReactive(input$custom_create_heatmap,{
    validate(
      need(input$custom_cancer != "", "Please select a Cancer")
    )
    validate(
      need(input$custom_cancer != "Acute Myeloid Leukemia Marrow", "Acute Myeloid Leukemia marrow does not have heiarchical clustering data for analysis")
    )
    validate(
      need(input$custom_pathway1 != "", "Please select 1st Pathway")
    )
    req(input$custom_cancer)
    req(input$custom_pathway1)
    mapping=load_mapping()
    cancer_name=paste(strsplit(input$custom_cancer," ")[[1]],collapse="_")
    cancer_index=which(mapping$Full_Name==cancer_name)
    cancer=as.vector(mapping[cancer_index,"Symbol"])
    pathway=input$custom_pathway1
    base_dir=data_dir
    hm <- readRDS(file.path(base_dir,"heatmap_base",paste0(cancer,"_heatmap.RDS")))
    base_ano <- readRDS(file.path(base_dir,"dendro_annotations",paste0(cancer,"_annotation.RDS")))
    #tsne_ano <- readRDS(file.path(base_dir,"tsne_annotations",cancer,paste0(pathway,"#tsneAnnotation.RDS")))
    custom_heatmap_gen <- draw(base_ano %v% hm)
  })
  
  
  
  # User Input Buttons ------------------------------------------------------
  
  #################### pre-run t-SNE Groups Analysis 1 ####################

  #### Step 0: Select Cancer and Gene Sets (Sidebar)
  output$cancerSelector <- renderUI({
    Cancers=list.dirs(tsne_table_dir,full.names = FALSE,recursive = TRUE) #Get list of cancers from tsne_table_dir
    Cancers=unname(sapply(Cancers,function(x){return(paste(strsplit(x,"_")[[1]],collapse=" "))}))
    Cancers=c("",Cancers)
    selectizeInput("cancer",label="Choose a Cancer",selected=Cancers[1],multiple=FALSE,choices=Cancers) #Cancers[1] vs NA
  })
  output$pathway1Selector <- renderUI({
    req(input$cancer)
    cancer_name=paste(strsplit(input$cancer," ")[[1]],collapse="_")
    Pathways=c("",as.vector(sapply(list.files(paste(tsne_table_dir,cancer_name,sep="/"),full.names = FALSE),h1)))
    selectizeInput("pathway1",label="Choose 1st Pathway",selected=Pathways[1],multiple=FALSE,choices=Pathways) #NA Pathways[1]
  })
  output$pathway2Selector <- renderUI({
    req(input$pathway1)
    cancer_name=paste(strsplit(input$cancer," ")[[1]],collapse="_")
    Pathways2=c("",as.vector(sapply(list.files(paste(tsne_table_dir,cancer_name,sep="/"),full.names = FALSE),h1)))
    pathways2a=Pathways2[which(Pathways2!=input$pathway1)]
    selectizeInput("pathway2",label="Choose 2nd Pathway",selected=pathways2a[1],multiple=FALSE,choices=pathways2a)  #NA pathways2a[1]
  })
  
  #### Step 2a: Select Clusters in Gene Sets (Tab 3)
  output$clusterSelector <- renderUI({
    req(input$pathway2)
    tsne_table1=load_tsne_table1()
    Clusters=sort(unique(tsne_table1$Cluster))
    Clusters1=c("all",Clusters)
    selectizeInput("clusters",label="Choose Clusters in Pathway 1",multiple=TRUE,selected=NA,choices=Clusters1)
  })
  output$clusterSelector2 <- renderUI({
    req(input$clusters)
    tsne_table1=load_tsne_table2()
    Clusters=sort(unique(tsne_table1$Cluster))
    Clusters1=c("all",Clusters)
    selectizeInput("clusters2",label="Choose Clusters in Pathway 2",multiple=TRUE,selected="all",choices=Clusters1)
  })
  
  #### Step 2b: Select Clusters in Dendrogram Groups (Tab 4)
  output$clusterSelector_heatmap <- renderUI({
    req(input$pathway1)
    tsne_table1=load_tsne_table1()
    Clusters=sort(unique(tsne_table1$Cluster))
    Clusters1=c("all",Clusters)
    selectizeInput("clusters_heatmap",label="Choose Clusters in Pathway 1",multiple=TRUE,selected=NA,choices=Clusters1)
  })
  output$groupSelector2 <- renderUI({
    req(input$clusters_heatmap)
    req(input$cancer!="Acute Myeloid Leukemia Marrow")
    
    table=load_dendrogram_groups()
    Groups=sort(unique(table[,'Group']))
    Groups1=c("all",Groups)
    selectizeInput("group2",label="Choose Dendogram Groups",multiple=TRUE,selected=NA,choices=Groups1)
  })
  output$go_to_web <- renderUI({
    validate(
      need(input$cancer != "", "Please select a Cancer")
    )
    #fluidRow(column(4,conditionalPanel(condition = "input.cancer != ''",uiOutput("go_to_web"))))
    actionButton("go_to_web", "View on Next-Generation Clustered Heat Map (NG-CHM) Viewer")
  })
  
  output$create_heatmap_button<- renderUI({
    req(input$pathway1)
    actionButton("create_heatmap", "Create Heatmap")
  })
  
  #### Step 1: Filter patients in original analysis (Tab 1)
  
  ######## Checkbox
  output$filter_option <- renderUI({
    req(input$pathway1)
    checkboxInput("filter_control",label="Filter (Y/N)",value=FALSE,width=NULL)
  })
  
  ######## Choose Variables 
  output$filter_pheno_selector <- renderUI({
    req(input$pathway1)
    table=load_tsne_table1()
    variables=c("",colnames(table)[c(9,11)])
    
    selectizeInput("phenotype",label="Phenotype(s)",multiple=TRUE,selected=variables[1],choices=variables)
  })

  ######## Dynamic Selectors for threshold and comparision type (3 sets - 3 possible variables compared at a time)
  output$filter_comptype_selector1 <- renderUI({
    req(input$phenotype != "")
    req(length(input$phenotype) >= 1)
    i=1
    if(length(input$phenotype)>1){phenotype=input$phenotype[i]} else { phenotype=input$phenotype}
    CD1=load_clinical_data1()
    clinical_data_types=as.vector(CD1[1,])
    type=clinical_data_types[which(colnames(CD1)==phenotype)] 
    if (type == "integer"){
      options=c("none","<","<=",">",">=")
    }
    else if (type =="string"){
      options=c("none","include","exclude")
    }
    selectizeInput("compare1",label=paste(phenotype," - Option",sep=" "),multiple=FALSE,selected=options[1],choices=options)
  })
  output$filter_comptype_selector2 <- renderUI({
    req(input$phenotype != "")
    req(length(input$phenotype) >= 2)
    phenotype=input$phenotype[2]
    CD1=load_clinical_data1()
    clinical_data_types=as.vector(CD1[1,])
    type=clinical_data_types[which(colnames(CD1)==phenotype)] 
    if (type == "integer"){
      options=c("none","<","<=",">",">=")
    }
    else if (type =="string"){
      options=c("none","include","exclude")
    }
    selectizeInput("compare2",label=paste(phenotype," - Option",sep=" "),multiple=FALSE,selected=options[1],choices=options)
  })
  output$filter_pheno_value1 <- renderUI({
    req(input$phenotype != "")
    req(length(input$phenotype) >= 1)
    i=1
    if(length(input$phenotype)>1){phenotype=input$phenotype[i]} else { phenotype=input$phenotype}
    CD=load_tsne_table1()
    CD1=load_clinical_data1()
    clinical_data_types=as.vector(CD1[1,])
    index=which(colnames(CD)==phenotype)
    column=CD[,index]
    type=clinical_data_types[which(colnames(CD1)==phenotype)] 
    if(type=="integer"){
      column=as.numeric(column)
      range_max=max(column)
      range_min=min(column)
      if(phenotype=="Age"){
        range_max=round(range_max/365)+1
        range_min=round(range_min/365)-1
      }
      sliderInput("threshold1",label=paste(phenotype,"- Threshold",sep=" "),min=range_min,max=range_max,value="")
    }
    else if(type=="string"){
      options=c(" ",as.vector(unique(column)))
      selectizeInput("threshold1",label=paste(phenotype,"- Threshold",sep=" "),multiple=TRUE,selected=options[1],choices=options)
    }
  })
  output$filter_pheno_value2 <- renderUI({
    req(length(input$phenotype) >= 2)
    req(input$phenotype[1] != "")
    phenotype=input$phenotype[2]
    CD=load_tsne_table1()
    CD1=load_clinical_data1()
    clinical_data_types=as.vector(CD1[1,])
    index=which(colnames(CD)==phenotype)
    column=CD[,index]
    type=clinical_data_types[which(colnames(CD1)==phenotype)] 
    if(type=="integer"){
      column=as.numeric(column)
      range_max=max(column)
      range_min=min(column)
      if(phenotype=="Age"){
        range_max=round(range_max/365)+1
        range_min=round(range_min/365)-1
      }
      sliderInput("threshold2",label=paste(phenotype,"- Threshold",sep=" "),min=range_min,max=range_max,value="")
    }
    else if(type=="string"){
      options=c(" ",as.vector(unique(column)))
      selectizeInput("threshold2",label=paste(phenotype,"- Threshold",sep=" "),multiple=TRUE,selected=options[1],choices=options)
    }
  })
  
  #################### tSNE/UMAP Analysis 2 ####################
  output$custom_cancerSelector <- renderUI({
    Cancers=list.dirs(tsne_table_dir,full.names = FALSE,recursive = TRUE) #Get list of cancers from tsne_table_dir
    Cancers=Cancers[-2]
    print(head(Cancers))
    Cancers=unname(sapply(Cancers,function(x){return(paste(strsplit(x,"_")[[1]],collapse=" "))}))
    Cancers=c("",Cancers)
    selectizeInput("custom_cancer",label="Choose a Cancer",selected=Cancers[1],multiple=FALSE,choices=Cancers) #Cancers[1] vs NA
  })
  
  output$custom_pathway1_selector<- renderUI({
    req(input$custom_cancer)
    gene_pathways=load_gene_pathways()
    Pathways=c("","Custom Gene Set",as.vector(gene_pathways$Pathway))
    selectizeInput("custom_pathway1",label="Pathway 1",selected=Pathways[1],multiple=FALSE,choices=Pathways) 
  })
  
  output$custom_pathway2_selector<- renderUI({
    req(input$custom_pathway1)
    gene_pathways=load_gene_pathways()
    Pathways=as.vector(gene_pathways$Pathway)
    Pathways2=Pathways[which(Pathways!=input$custom_pathway1)]
    Pathways3=c("","Custom Gene Set",Pathways2)
    selectizeInput("custom_pathway2",label="Pathway 2",selected=Pathways3[1],multiple=FALSE,choices=Pathways3) 
  })
  
  output$custom_pathway1_multigene_selector <- renderUI({
    req(input$custom_pathway1)
    req(input$custom_pathway1=="Custom Gene Set")
    gene_pathways=load_gene_pathways()
    Genes=c("",as.vector(gene_pathways$Gene))
    selectizeInput("multigene_p1",label="Select Custom Gene Set as Pathway 1",selected=Genes[1],multiple=TRUE,choices=Genes)
  })
  
  output$custom_pathway2_multigene_selector <- renderUI({
    req(input$custom_pathway2)
    req(input$custom_pathway2=="Custom Gene Set")
    gene_pathways=load_gene_pathways()
    Genes=c("",as.vector(gene_pathways$Gene))
    selectizeInput("multigene_p2",label="Select Custom Gene Set as Pathway 2",selected=Genes[1],multiple=TRUE,choices=Genes)
  })
  
  output$run_custom<- renderUI({
    req(input$custom_pathway1)
    actionButton("run_custom", "Run Analysis!")
  })
  
  output$exp_data_type<- renderUI({
    selectInput("use_tpm",label="Choose expression data type",choices=c("","TPM","FPKM"),selected="FPKM")
  })
  
  output$run_algo<- renderUI({
    selectInput("use_UMAP",label="Choose algorithm to use",choices=c("","UMAP","t-SNE"),selected="UMAP")
  })
  
  output$custom_create_heatmap_button<- renderUI({
    req(input$custom_pathway1)
    actionButton("custom_create_heatmap", "Create Heatmap")
  })
  
  #### Step 1: Filter patients in custom analysis before tSNE/UMAP (Tab 4)
  
  output$filter_option_custom <- renderUI({
    req(input$custom_cancer)
    options=c("No Filter","Filter prior to clustering","Filter patients after clustering")
    selectizeInput("filter_control_custom",label="Choose Filtering Options",multiple=FALSE,selected=options[1],choices=options)
  })
  
  output$filter_pheno_selector_custom <- renderUI({
    req(input$custom_cancer)
    dc_in=load_custom_analysis_cancer_data()
    data=dc_in$metadata_filtered
    variables=c("",colnames(data)[c(8,9)])
    #print(variables)
    selectizeInput("custom_phenotypes",label="Phenotype(s)",multiple=TRUE,selected=variables[1],choices=variables)
  })
  
  output$filter_comptype_custom_selector1 <- renderUI({
    req(input$custom_phenotypes != "")
    req(length(input$custom_phenotypes) >= 1)
    i=1
    if(length(input$custom_phenotypes)>1){phenotype=input$custom_phenotypes[i]} else { phenotype=input$custom_phenotypes}
    
    CD=load_custom_analysis_cancer_data()
    CD1=CD$metadata_filtered
    
    clinical_data_types=load_phenotype_map_custom() 
    
    type=clinical_data_types[which(colnames(CD1)==phenotype)] 
    if (type == "integer"){
      options=c("none","<","<=",">",">=")
    }
    else if (type =="string"){
      options=c("none","include","exclude")
    }
    selectizeInput("custom_compare1",label=paste(phenotype," - Option",sep=" "),multiple=FALSE,selected=options[1],choices=options)
  })
  
  output$filter_custom_pheno_value1 <- renderUI({
    req(input$custom_phenotypes != "")
    req(length(input$custom_phenotypes) >= 1)
    i=1
    if(length(input$custom_phenotypes)>1){phenotype=input$custom_phenotypes[i]} else { phenotype=input$custom_phenotypes}
    
    CD=load_custom_analysis_cancer_data()
    CD1=as.matrix(CD$metadata_filtered)
    clinical_data_types=load_phenotype_map_custom() 
    
    index=which(colnames(CD1)==phenotype)
    column=CD1[,index]
    type=clinical_data_types[index] 
    if(type=="integer"){
      column=as.numeric(column)
      range_max=max(column)
      range_min=min(column)
      if(phenotype=="Age"){
        range_max=round(range_max/365)+1
        range_min=round(range_min/365)-1
      }
      sliderInput("custom_threshold1",label=paste(phenotype,"- Threshold",sep=" "),min=range_min,max=range_max,value="")
    }
    else if(type=="string"){
      options=c(" ",as.vector(unique(column)))
      selectizeInput("custom_threshold1",label=paste(phenotype,"- Threshold",sep=" "),multiple=TRUE,selected=options[1],choices=options)
    }
  })
  
  ##### Seqeuntial Analysis 
  
  output$custom_group_clusterSelector <- renderUI({
    req(input$custom_pathway2)
    table=custom_analysis_pathway1()
    validate(
      need(table!="FAILED", "")
    )
    Clusters=sort(unique(table$Cluster))
    Clusters1=c("all",Clusters)
    selectizeInput("clusters_custom",label="Choose Clusters in Pathway 1",multiple=TRUE,selected=NA,choices=Clusters1)
  })
  
  output$custom_group_clusterSelector2 <- renderUI({
    req(input$clusters_custom)
    table=custom_analysis_pathway2()
    validate(
      need(table!="FAILED", "")
    )
    Clusters=sort(unique(table$Cluster))
    Clusters1=c("all",Clusters)
    selectizeInput("clusters2_custom",label="Choose Clusters in Pathway 2",multiple=TRUE,selected="all",choices=Clusters1)
  })
  
  ##### Dendrogram Analysis 
  
  
  output$custom_clusterSelector_heatmap <- renderUI({
    req(input$custom_pathway1)
    table=custom_analysis_pathway1()
    Clusters=sort(unique(table$Cluster))
    Clusters1=c("all",Clusters)
    selectizeInput("custom_clusters_heatmap",label="Choose Clusters in Pathway 1",multiple=TRUE,selected=NA,choices=Clusters1)
  })
  
  output$custom_groupSelector2 <- renderUI({
    req(input$custom_clusters_heatmap)
    req(input$custom_cancer!="Acute Myeloid Leukemia Marrow")
    
    table=load_dendrogram_groups_custom()
    Groups=sort(unique(table[,'Group']))
    Groups1=c("all",Groups)
    selectizeInput("custom_group2",label="Choose Dendogram Groups",multiple=TRUE,selected=NA,choices=Groups1)
  })
  
  ###### Algo Options
  
  #t-SNE
  
  output$tsne_input_max_iter <- renderUI({
    sliderInput("max_iter",label="Iterations",min=0,max=5000,step=10,value=500)
  })
  
  output$tsne_input_perplexity <- renderUI({
    sliderInput("perplexity",label="Perplexity",min=0,max=100,step=1,value=10)
  })

  #UMAP
  
  output$umap_input_n_epochs <- renderUI({
    textInput("n_epochs",label="# Epochs",value=500)
  })
  
  output$umap_input_min_dist<- renderUI({
    sliderInput("min_dist",label="Minimum Distance (10^-(choosen value))",min=-500,max=1,step=1,value=-66)
  })
  
  output$umap_n_neighbors <- renderUI({
    textInput("n_neighbors",label="# Neighbors",value=15)
  })
  
  #Spherize (used during normalization)
  
  
  #hdb scan
  
  
  
  
  # Output Text (Analysis) -------------------------------------------------------------
  
  #################### pre-run t-SNE Groups Analysis 1 ####################
  
  output$text1 <- renderText(input$pathway1)
  #output$text2 <- renderText("Number of Patients in Each Cluster")
  
  # Tab 3
  #"Number of patients in each cluster per pathways"
  output$text2 <- renderText({
    req(input$clusters2);
    paste(c("<h5>","<b>Number of patients in each cluster per pathways<b>","</h5>"),collapse="")})
  #"Sequential analysis with" 
  output$survivalplotmid<- renderText({
    req(input$clusters2);
    "sequential analysis with"})
  output$survivalplotheader <- renderText({
    req(input$clusters2);
    cluster=paste(input$clusters,collapse=" , ")
    paste(c("<h5>",paste(c("t-SNE clusters in ",input$pathway1,": <b>",cluster,'<b>'),collapse=" "),"</h5>"),collapse="")})
  output$survivalplotheader2 <- renderText({
    req(input$clusters2);
    cluster=paste(input$clusters2,collapse=" , ")
    paste(c("<h5>",paste(c("t-SNE clusters in",input$pathway2,": <b>",cluster,'<b>'),collapse=" "),"</h5>"),collapse="")})
  
  # Tab 4 
  #"Number of patients in each cluster per pathway"
  output$text3 <- renderText({
    req(input$group2);
    paste(c("<h5>","<b>Number of patients in each t-SNE cluster and dendrogram groups<b>","</h5>"),collapse="")})
  output$heatmaptitle<- renderText({
    req(input$pathway1);
    paste(c("<h4><b>",input$cancer,"<b><h4>"),collapse="")})  
  output$heatmaptitle2<- renderText({
    req(input$pathway1);
    paste(c("<h4><b>",input$pathway1,"<b><h4>"),collapse="")})  
  output$dendrosurvivalplotmid<- renderText({
    req(input$group2);
    "sequential analysis with"})
  output$dendrosurvivalplotheader <- renderText({
    req(input$group2);
    cluster=paste(input$clusters,collapse=" , ")
    paste(c("<h5>",paste(c("t-SNE clusters in ",input$pathway1,": <b>",cluster,'<b>'),collapse=" "),"</h5>"),collapse="")})
  output$dendrosurvivalplotheader2 <- renderText({
    req(input$group2);
    cluster=paste(input$group2,collapse=" , ")
    paste(c("<h5>",paste(c("Dendrogram Groups: <b>",cluster,'<b>'),collapse=" "),"</h5>"),collapse="")})
  
  # Tab 2 - "Choose Phenotype to Filter" 
  output$text_filter <- renderText({
    req(input$pathway1);
    paste(c("<h5>","<b>Choose Phenotype(s) to Filter<b>","</h5>"),collapse="")})
  output$pathway1_title_3d <- renderText({
    req(input$pathway1);
    paste(c("<h5>","<b>Pathway  1 : <b>",input$pathway1,"</h5>"),collapse="")})
  output$pathway2_title_3d <- renderText({
    req(input$pathway2);
    paste(c("<h5>","<b>Pathway 2 : <b>",input$pathway2,"</h5>"),collapse="")})
  
  
  # Tab 1 - Cancer Name
  output$cancer_name_about <- renderText({
    req(input$cancer);
    paste(c("<h4>",input$cancer,"</h4>"),collapse="")})
  output$cancer_name <- renderText({
    req(input$pathway1);
    paste(c("<h4><b>",input$cancer,"<b></h4>"),collapse="")})
  
  #Misc
  output$pdfview <- renderUI({
    tags$iframe(style="height:600px; width:100%", src="Test.pdf")
  })
  
  output$appdes <- renderText({'<h6>TEST<h6>/'})
  
  #################### tSNE/UMAP Analysis 2 ####################
  output$custom_test <-renderText({
    req(input$custom_pathway1)
    req(input$custom_pathway2)
    custom_test()
  })
  
  output$text_filter_custom <- renderText({
    req(input$pathway1);
    paste(c("<h5>","<b>Choose Phenotype(s) to Filter<b>","</h5>"),collapse="")})
  
  output$custom_cancer_name <- renderText({
    req(input$custom_pathway1);
    req(input$run_custom)
    paste(c("<h4><b>",input$custom_cancer,"<b></h4>"),collapse="")})
  
  output$custom_pathway1_title_3d <- renderText({
    req(input$custom_pathway1);
    req(input$run_custom)
    paste(c("<h5>","<b>Pathway  1 : <b>",input$custom_pathway1,"</h5>"),collapse="")})
  
  output$custom_pathway2_title_3d <- renderText({
    req(input$custom_pathway2);
    req(input$run_custom)
    paste(c("<h5>","<b>Pathway 2 : <b>",input$custom_pathway2,"</h5>"),collapse="")})

  output$text_run_custom_warning <- renderText({
    req(input$custom_pathway1)
    paste(c("<h5>","Analysis may take a few seconds to run","</h5>"),collapse="")})
  
  
  # Output Tables ------------------------------------------------------------
  
  #################### pre-run t-SNE Groups Analysis 1 ####################
  output$table <- DT::renderDataTable({ # Not currently used 
    #table <- display_data()
    table <- load_tsne_tables_dendro()[[1]]
    DT::datatable(table,filter="top",rownames=FALSE)
  })
  
  output$counttable <- DT::renderDataTable({
    table <- counttable()
    print(table)
    DT::datatable(table,rownames=TRUE)
  })
  
  #output$counttable <- renderTable({counttable()},bordered = TRUE,align="c",digits=0,rownames = TRUE)

  output$counttable_scatter_p1 <- renderTable({count_table_scatter1()},bordered = TRUE,align="c",digits=0,rownames = TRUE)
  output$counttable_scatter_p2 <- renderTable({count_table_scatter2()},bordered = TRUE,align="c",digits=0,rownames = TRUE)
  
  output$counttable_dendro <- DT::renderDataTable({
    table <- counttable_dendro()
    DT::datatable(table,rownames=TRUE)
  })
  
  
  #################### tSNE/UMAP Analysis 2 ####################
  output$counttable_custom_scatter_p1 <- renderTable({count_table_custom_scatter1()},bordered = TRUE,align="c",digits=0,rownames = TRUE)
  output$counttable_custom_scatter_p2 <- renderTable({count_table_custom_scatter2()},bordered = TRUE,align="c",digits=0,rownames = TRUE)
  
  output$custom_counttable <- DT::renderDataTable({
    table <- counttable_custom()
    print(table)
    DT::datatable(table,rownames=TRUE)
  })
  
  output$custom_counttable_dendro <- DT::renderDataTable({
    table <- custom_counttable_dendro()
    DT::datatable(table,rownames=TRUE)
  })
  
  # Output Plots ------------------------------------------------------------
  
  
  #################### pre-run t-SNE Groups Analysis 1 ####################
  
  #### Tab 2
  output$survivalplot_P1_out <- renderPlot(survivalplot_P1_gen())
  output$survivalplot_P2_out <- renderPlot(survivalplot_P2_gen())
  output$scatterplot_out <- renderPlotly({scatterplot_gen2_p1()})
  output$scatterplot_out2 <- renderPlotly({scatterplot_gen2_p2()})
  
  #### Tab 3  
  output$survivalplot_out <- renderPlot(survivalplot_gen())
  
  #### Tab 4
  output$survivalplot_dendro_out <- renderPlot(survivalplot_dendo_groups_gen())
  output$heatmap <- renderPlot(heatmap_gen())
  
  #################### tSNE/UMAP Analysis 2 ####################
  output$custom_scatterplot_P1_out <- renderPlotly(custom_scatterplot_gen2_p1())
  
  output$custom_survivalplot_P1_out <- renderPlot(custom_survivalplot_P1_gen())
  
  output$custom_scatterplot_P2_out <- renderPlotly(custom_scatterplot_gen2_p2())
  
  output$custom_survivalplot_P2_out <- renderPlot(custom_survivalplot_P2_gen())
  
  output$custom_survivalplot_out <- renderPlot(custom_survivalplot_gen())
  
  output$custom_survivalplot_dendro_out <- renderPlot(custom_survivalplot_dendo_groups_gen())
  
  output$custom_heatmap <- renderPlot(custom_heatmap_gen())
  
  
  
}#/server