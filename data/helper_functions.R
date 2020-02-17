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

#Merges 1 loaded tsne_table file and clincal data
merge_clinical=function(data,clinicalD){
  return(merge(data,clinicalD,by.x="Sample",by.y="samples",all.x=TRUE))
} 

#Merges 2 TSNE tables and creates labels for survivial curves 
sequential_filter=function(data_1,data_2,Cluster,reverse=FALSE,Pathway_Names,Cluster2){
  if (reverse==TRUE){
    data1=data_2
    data2=data_1
    Pathway_Names=c(Pathway_Names[2],Pathway_Names[1])
  } else {
    data1=data_1
    data2=data_2
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


# Survival Curve - Fits (survival) ----------------------------------------
#Creates survival "FIT" for sequential analysis
my_sfit=function(data){
  
  # if ("all" %in% Cluster){colnames(data)[which(colnames(data)=="Merged")]="C" #Changes labels used for curves
  # } else {colnames(data)[which(colnames(data)=="Merged2")]="C"}
  colnames(data)[which(colnames(data)=="Merged2")]="C"
  fit<-surv_fit(Surv(OS.time,OS)~ C, data=data)
  return(list(fit,data))
} 

#Creates original survival "FIT" (NO SEQUENTIAL)
my_sfit_P1=function(data){
  fit<-surv_fit(Surv(OS.time,OS)~ Cluster, data=data)
  return(list(fit,data))
} 

#returns the pval from a generated "FIT"
my_spval=function(fit,data){return(round(surv_pvalue(fit,data=data)$pval,3))} 


# Plots -------------------------------------------------------------------
my_splot_create1=function(fit,data,Pathway_Names,Cluster,pal){
  if("all" %in% Cluster){ C=Cluster
  } else {C=paste(Cluster,collapse =",")}
  pal=pal[1:length(unique(data$C))]
  #title=paste(c("Pathway1:",Pathway_Names[1],"Pathway2:",Pathway_Names[2],"Cluster:",C),collapse =" ")
  return(ggsurvplot(fit, data=data,pval=TRUE,title='',font.tile=1,font.x=8,font.y=8,palette = pal))
} 

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

#Counts number of patients represented in ach curve of created survival curve 
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


