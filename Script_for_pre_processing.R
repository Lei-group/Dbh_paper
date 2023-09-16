###This script is for R-based pre-processing and preliminary analysis of single-cell RNA seq R objects used in the paper.
#Data was processed as in methods and is loaded into R first as raw count matrices.
#We then process these count matrices, and we include our preliminary analysis to demonstrate how we first engaged and interacted with the data.
#We then process this data into our final processed R object of the filtered and clustered scRNAseq data.
#We include representative pseudo-code of the clustering processes - as these were done in an interactive and prolonged manner and not recorded in the script.
#However, we provide the final clustered and labeled objects. 
#We then subset this dataset to produce the cardiomyocyte lineage (CM lineage) and then the Dbh positive CM lineage (Dbh p).
#We then combine parts of the whole scRNAseq dataset, with CM-specific clusters to produce a dataset for intergration with ST data by BGI.
#We present each figure independently in 'chronological order'.
#We include a number of 'save points' in the script as these large files can caught excessive instability in R and result in crashes.
#This script is particularly resource-intensive and was run on a c. 512 GB RAM machine. We apologise for any lack of optimisation.
#If you are interested in just replicating a figure: we have provided our final processed R objects and the scripts needed to use them, in the additional Figures.R script.
#This script demonstrates representative, but not exhaustive, code of our pre-processing workflow and generates a number of figures present in the paper and supplementary figures, and also includes a number of figures not included in the paper but left here for the interest of discerning readers.

#### Read in libraries and set-up global variables

#### Read in count matrices

#### QC and filtering

#### Preliminary exploration of data

#### Final whole scRNAseq dataset construction

#### Cardiomyocyte lineage dataset construction

#### Dbhp cardiomyocyte dataset construction

#### Combined dataset construction for integration with Spatial data


####################READ IN LIBRARIES AND SET UP GLOBAL VARIABLES
library(reticulate)
use_python("#####/anaconda3/") #this should be wherever you have your installation of anaconda3
library(tidyverse)
library(ggplot2)
library(stringr)
library(dplyr) 
library(BiocGenerics)
library(tidyr)
library(R.utils)
library(Matrix)
library(Seurat)
library(MASS)
library(sctransform)
library(base)
library(tibble)
library(sctransform)
library(phateR)
library(ComplexHeatmap)
library(viridis)
library(plotly)
library(htmlwidgets)
library(ComplexHeatmap)
library(stringr)
library(viridis)
library(RColorBrewer)
library(circlize)
#library(robustbase)
library(Seurat)
library(utils)
library(slingshot)
library(rgl)
library(tradeSeq)
library(BiocParallel)
library(hablar)
library(scone)
library(parglm)
library(devtools)
library(kBET)
library(processx)
library(bayNorm)
library(devtools)
library(harmony)
library(devtools)
library(multtest)
library(metap)
library(Rmagic)
library(forcats)
library(see)
library(ggh4x)
pyphate <- reticulate::import("phate")

#Set global options for multiprocessing
options(SnowParam=SnowParam(workers=25, type=c('SOCK'),tasks=0L,progressbar=TRUE,RNGseed=804L,exportglobals = TRUE,log=TRUE))
register(SnowParam(workers=25, type=c('SOCK'),tasks=0L,progressbar=TRUE,RNGseed=804L,exportglobals = TRUE,log=TRUE))
set.seed(804L)


RAW_FILE_LOC='###' #DIR WHERE RAW FILES DOWNLOADED

OUTPUT_FILE_LOC='###' #DIR WHERE YOU WANT OUTPUT FILES TO GO - SUBDIRS ARE NOT CREATED AUTOMATICALLY BY THIS SCRIPT

INT_FILE_LOC='###' #WHERE YOU WANT INTERMEDIATE FILES TO GO E.G. PROCESSING LARGE DATA.



#I make regular 'save points', as R doesn't always like working with objs this big
#And so can be prone to crashing - please try and avoid reading in your file just after saving it, or this will just waste your time as these take a while to read/write


#### Read in count matrices

#read in the raw datasets - warning: these are very large and this will not be possile on a standard laptop etc.
stage=c('E8_5','E10_5','E12_5','E14_5','E16_5','Neo_')
inlist <- list()
for(i in 1:6){
  batchlist <- list()
  for(j in 1:2){
    batchlist[[j]] = data.table::data.table::fread(paste0(RAW_FILE_LOC,"/",stage[i],"_",j,"_Spliced.csv"),header=TRUE,check.names=FALSE,logical01=FALSE)  #this should be wherever you have saved the composite stage-wise sections
    batchlist[[j]] = t(batchlist[[j]])
    colnames(batchlist[[j]]) = make.unique(batchlist[[j]]['V1',])
    batchlist[[j]]=batchlist[[j]][,intersect(colnames(batchlist[[j]]),batchlist[[j]]['V1',])][2:nrow(batchlist[[j]]),]
    batchlist[[j]] = as.matrix(batchlist[[j]])
    batchlist[[j]] = as.sparse(batchlist[[j]])
  }
  inlist[[i]] = batchlist
}
#bind together per stage
outlist=list()
seurout=list()
for(k in 1:6){
  if(k=6){
    rownames(inlist[[k]][[1]]) = paste0("Neo_1",str_sub(rownames(inlist[[k]][[1]]),5L,30L))
    rownames(inlist[[k]][[2]]) = paste0("Neo_2",str_sub(rownames(inlist[[k]][[2]]),5L,30L))
  }
  outlist[[k]] = rbind2(inlist[[k]][[1]],inlist[[k]][[2]])
  
  
  #normalise
  outlist[[k]] <- t(outlist[[k]])
  spliced <- CreateSeuratObject(counts = outlist[[k]], project = "spliced", meta.data=NULL)
  spliced@meta.data$nUMI <- colSums(outlist[[k]])
  spliced@meta.data$nGene <-  colSums(outlist[[k]]>0)
  spliced@meta.data$batch <- str_sub(colnames(outlist[[k]]),-18L,-18L)
  spliced@meta.data$log10GenesPerUMI <- log10(spliced@meta.data$nGene) / log10(spliced@meta.data$nUMI)
  spliced@meta.data$mtUMI <- Matrix::colSums(outlist[[k]][grep(pattern="^mt",x=rownames(outlist[[k]]),ignore.case=FALSE,value=TRUE),], na.rm=TRUE)
  spliced@meta.data$mitoRatio <- spliced@meta.data$mtUMI/spliced@meta.data$nUMI
  seurout[[k]] = spliced
  pdf(paste0(OUTPUT_FILE_LOC,"/QC/QC_",stage[k],".pdf")) #wherever you want your output to be stored. This is preliminary QC based on entire dataset.
  print(VlnPlot(spliced, features = c("nUMI"),split.by = 'batch'))
  print(VlnPlot(spliced, features = c("nGene"),split.by = 'batch'))
  print(VlnPlot(spliced, features = c("mitoRatio"),split.by = 'batch'))
  print(VlnPlot(spliced, features = c("mtUMI"),split.by = 'batch'))
  dev.off()
}

#outlist is raw combined by stage
#seurout is combined by stage and then into seurat objects
saveRDS(seurout,paste0(INT_FILE_LOC,"/stages_seur_prefilt.csv")) # INT_FILE_LOCATION IS WHEREVER YOU WOULD LIKE TO STORE INTERMEDIATE FILES
saveRDS(outlist,paste0(INT_FILE_LOC,"/stages_tab_prefilt.csv"))


####QC AND FILTERING WITH DYNAMIC MT RANGE ACROSS DEVELOPMENT
mtbar = c(0.05,0.05,0.075,0.10,0.15,0.20)
outlistz = list()
for(l in 1:6){
  
  metadata_cell_ids <- data.frame(row.names=colnames(outlist[[l]]), cells=colnames(outlist[[l]]), stringsAsFactors=FALSE)
  metadata_cell_ids$nUMI <- colSums(outlist[[l]])
  metadata_cell_ids$nGene <- colSums(outlist[[l]]>0)
  metadata_cell_ids$batch <-  str_sub(colnames(outlist[[l]]),-18L,-18L)
  metadata_cell_ids$log10GenesPerUMI <-  log10(metadata_cell_ids$nGene)/log10(metadata_cell_ids$nUMI)
  metadata_cell_ids$mtUMI <- Matrix::colSums(outlist[[l]][grep(pattern="^mt",x=rownames(outlist[[l]]),ignore.case=FALSE,value=TRUE),], na.rm=TRUE)
  metadata_cell_ids$mtUMI[is.na(metadata_cell_ids$mtUMI)] <- 0
  metadata_cell_ids$mitoRatio <- metadata_cell_ids$mtUMI/metadata_cell_ids$nUMI
  filtered_data <- metadata_cell_ids %>% dplyr::filter(nUMI > 300, nGene > 270, log10GenesPerUMI >0.8, mitoRatio < mtbar[[l]]) %>% pull(cells)
  outlistz[[l]] = outlist[[l]][,intersect(filtered_data,colnames(outlist[[l]]))]
}
saveRDS(outlistz, paste0(INT_FILE_LOC,"/stages_tab_filt.csv")) #this is now a list of filtered seurat objects

seuroutz =list()
for(m in 1:6){
  spliced <- CreateSeuratObject(counts = outlistz[[m]], project = "spliced", meta.data=NULL)
  spliced@meta.data$nUMI <- colSums(outlistz[[m]])
  spliced@meta.data$nGene <-  colSums(outlistz[[m]]>0)
  spliced@meta.data$batch <- str_sub(colnames(outlistz[[m]]),-18L,-18L)
  spliced@meta.data$log10GenesPerUMI <- log10(spliced@meta.data$nGene) / log10(spliced@meta.data$nUMI)
  spliced@meta.data$mtUMI <- Matrix::colSums(outlistz[[m]][grep(pattern="^mt",x=rownames(outlistz[[m]]),ignore.case=FALSE,value=TRUE),], na.rm=TRUE)
  spliced@meta.data$mitoRatio <- spliced@meta.data$mtUMI/spliced@meta.data$nUMI
  seuroutz[[m]] = spliced
  pdf(paste0(OUTPUT_FILE_LOC,"/QC/QC_",stage[m],"_postfilter.pdf")) #These are the outputs of the post-filtered dataset QC
  print(VlnPlot(spliced, features = c("nUMI"),split.by = 'batch'))
  print(VlnPlot(spliced, features = c("nGene"),split.by = 'batch'))
  print(VlnPlot(spliced, features = c("mitoRatio"),split.by = 'batch'))
  print(VlnPlot(spliced, features = c("mtUMI"),split.by = 'batch'))
  dev.off()
}
saveRDS(seuroutz,paste0(INT_FILE_LOC,"/stages_seur_postfilt.csv")) #These are the filtered seurat objects by stage


# OPTIONAL ASSAY TO PREDICT FNR AT EACH STAGE (Adapted from SCONE - https://www.bioconductor.org/packages/devel/bioc/vignettes/scone/inst/doc/sconeTutorial.html#drop-out-characteristics) 
data("housekeeping_revised")
hk = tolower(housekeeping_revised$V1)
# Define a color scheme
cc <- c(brewer.pal(9, "Set1"))
for(n in 1:length(hk)){
  hk[[n]] = paste0(toupper(str_sub(hk[[n]],0L,1L)),str_sub(hk[[n]],2L,nchar(hk[[n]])))
}
for(o in 1:6){
  hz = intersect(rownames(outlistz[[o]]),hk)  
  mu_obs=rowMeans(log10(outlistz[[o]][hz,]+1))
  drop_outs = outlistz[[o]][hz,]==0
  batch=factor(str_sub(colnames(outlistz[[o]]),-18L,-18L))
  ref.glms=list()
  for(si in 1:dim(drop_outs)[2]){
    fit=glm(cbind(drop_outs[,si],1 - drop_outs[,si])~mu_obs,family=binomial(logit))
    ref.glms[[si]] = fit$coefficients
  }
  par(mfrow=c(1,2))
  
  # Plot Failure Curves and Calculate AUC
  pdf(paste0(OUTPUT_FILE_LOC,"/QC/",stage[o],"False_negative_rates.pdf"))
  plot(NULL, main = "False Negative Rate Curves",
       ylim = c(0,1),xlim = c(0,6), 
       ylab = "Failure Probability", xlab = "Mean log10 Expression")
  x = (0:60)/10
  AUC = NULL
  for(si in 1:ncol(seuroutz[[o]]@assays$RNA)){
    y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
    AUC[si] = sum(y)/10
    lines(x, 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1),
          type = 'l', lwd = 2, col = cc[batch][si])
  }
  ord = order(AUC)[order(batch[order(AUC)])]
  
  barplot(AUC[ord], col=cc[batch][ord], border=cc[batch][ord], main="FNR AUC")
  legend("topright", legend=levels(batch), fill=cc, cex=0.4)
  dev.off()
}


#PRELIMINARY CLUSTERING AND OPTIONAL EXAMINATION OF BATCH EFFECT

outlistz = readRDS(paste0(INT_FILE_LOC,"/stages_tab_filt.csv"))
seuroutz=readRDS(paste0(INT_FILE_LOC,"/stages_seur_postfilt.csv"))
kbetoprim = list()
phateoutprim = list()
pmbcoutprim =list()
for(p in 1:6){
  #cluster in UMAP
  pbmc = FindVariableFeatures(seuroutz[[p]],nfeatures=5000)
  pbmc = ScaleData(pbmc,features=VariableFeatures(pbmc))
  pbmc = RunPCA(pbmc,features=VariableFeatures(pbmc))
  pbmc = FindNeighbors(pbmc,dims=1:16)
  pbmc = FindClusters(pbmc,resolution=0.5)
  pbmc = RunUMAP(pbmc,dims=1:16)
  #CLUSTERING USING SIMPLE NORMALISATION
  saveRDS(pbmc,paste0(INT_FILE_LOC,"/pbmc_",stage[p],"_postfilt.csv"))
  pdf(paste0(OUTPUT_FILE_LOC,"/QC/",stage[p],"Initial_clustering_pre_SCT.pdf"))
  print(DimPlot(pbmc,reduction='umap',group.by='batch'))
  print(DimPlot(pbmc,reduction='umap',group.by='seurat_clusters'))
  dev.off()
  pmbcoutprim[[p]] = pbmc
  #PHATE DIM REDUC AND CLUSTERING ON OPERATOR - SAVED AS INTERACTIVE WIDGETS FOR 2D INTERROGATION
  
  phateoutprim[[p]]=phate(t(outlistz[[p]]),ndim=20,n.jobs=-10,seed=804L)  
  phatez = cluster_phate(phateoutprim[[p]],k=length(pbmc@meta.data$seurat_clusters))
  saveRDS(phateoutprim[[p]],paste0(OUTPUT_FILE_LOC,"/QC/phate_",stage[p],"_postfilt.csv"))
  fig1=plot_ly(as.data.frame(phateoutprim[[p]]), mode="markers", x = ~PHATE1, y = ~PHATE2, color =pbmc@meta.data$batch, colors = viridis(n=2),
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1),
               text = ~paste('<br>PHATE1',PHATE1,'<br>PHATE2', PHATE2))
  fig2=plot_ly(as.data.frame(phateoutprim[[p]]), mode="markers", x = ~PHATE1, y = ~PHATE3, color =pbmc@meta.data$batch, colors = viridis(n=2),
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1),
               text = ~paste('<br>PHATE1',PHATE1,'<br>PHATE2', PHATE3))
  fig3=plot_ly(as.data.frame(phateoutprim[[p]]), mode="markers", x = ~PHATE2, y = ~PHATE3, color =pbmc@meta.data$batch, colors = viridis(n=2),
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1),
               text = ~paste('<br>PHATE1',PHATE2,'<br>PHATE2', PHATE3))
  htmlwidgets::saveWidget(as_widget(fig1), paste0(OUTPUT_FILE_LOC,"/QC/",stage[p],"Initial_phate_clustering_pre_SCT_1_2.html"))
  htmlwidgets::saveWidget(as_widget(fig2), paste0(OUTPUT_FILE_LOC,"/QC/",stage[p],"Initial_phate_clustering_pre_SCT_1_3.html"))
  htmlwidgets::saveWidget(as_widget(fig3), paste0(OUTPUT_FILE_LOC,"/QC/",stage[p],"Initial_phate_clustering_pre_SCT_2_3.html"))
  
  #OPTIONAL kBET TO QUANTIFY BATCH EFFECT IN QUANTITATIVE MANNER (https://github.com/theislab/kBET)
  #data: a matrix (rows: samples, columns: features (genes))
  #batch: vector or factor with batch label of each cell 
  #clusters: vector or factor with cluster label of each cell 
  subset_size <- 0.1 #subsample to 10% of the data
  subset_id <- sample.int(n = length(pbmc@meta.data$batch), size = floor(subset_size * length(pbmc@meta.data$batch)), replace=FALSE)
  batch.estimate <- kBET(t(outlistz[[p]])[subset_id,], pbmc@meta.data$batch[subset_id])
  kbetoprim[[p]] = batch.estimate
  saveRDS(batch.estimate,paste0(INT_FILE_LOC,"/kbet_",stage[p],"_postfilt.csv"))
  pdf(paste0(OUTPUT_FILE_LOC,"/QC/",stage[p],"Initial_kbet_pre_SCT.pdf"))
  plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                    each=length(batch.estimate$stats$kBET.observed)), 
                          data =  c(batch.estimate$stats$kBET.observed,
                                    batch.estimate$stats$kBET.expected))
  try(print(g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
              labs(x='Test', y='Rejection rate',title='kBET test results') +
              theme_bw() +  
              scale_y_continuous(limits=c(0,1))))
  dev.off()
  pca.data <- prcomp(t(outlistz[[p]])[subset_id], center=TRUE) #compute PCA representation of the data
  ###pca kbet next, check, also check PHATE
  batch.silhouette <- batch_sil(pca.data, pbmc@meta.data$batch[subset_id])
  batch.pca <- pcRegression(pca.data, pbmc@meta.data$batch[subset_id])
  saveRDS(batch.silhouette,paste0(INT_FILE_LOC,"/kbet_sil_",stage[q],"_preSCT.csv"))
  saveRDS(batch.pca,paste0(INT_FILE_LOC,"/kbet_pca_",stage[q],"_preSCT.csv"))
  
}
saveRDS(pmbcoutprim,paste0(INT_FILE_LOC,"/Seurat_clustering_prim.csv"))
saveRDS(phateoutprim,paste0(INT_FILE_LOC,"/Phate_prim.csv"))
saveRDS(kbetoprim,paste0(INT_FILE_LOC,"/kBET_prim.csv"))

#so can see how kBET is for all of them now
sctout=list()
phateout=list()
kbeto=list()

#Now use SCT to normalise them with a batch regression - NB, this was not used for final object, but forms part of our early exploration of the dataset.
for(q in 1:6){
  pbmc = readRDS(paste0(INT_FILE_LOC,"/pbmc_",stage[q],"_postfilt.csv"))
  sctout[[q]] <- SCTransform(pbmc, vars.to.regress = c("batch"), verbose = TRUE, min_cells=2, return.only.var.genes=FALSE,do.scale=TRUE)
  saveRDS(sctout[[q]],paste0(INT_FILE_LOC,"/SCT_",stage[q],"_filt.csv"))
  pbmc = RunPCA(sctout[[q]],assay='SCT')
  #print(ElbowPlot(pbmc))
  pbmc = FindNeighbors(pbmc,dims=1:16,assay='SCT')
  pbmc = FindClusters(pbmc,resolution=0.5,assay='SCT')
  pbmc = RunUMAP(pbmc,dims=1:16,assay='SCT')
  pdf(paste0(OUTPUT_FILE_LOC,"/QC/",stage[q],"_clustering_post_SCT.pdf"))
  print(DimPlot(pbmc,reduction='umap',group.by='batch'))
  print(DimPlot(pbmc,reduction='umap',group.by='seurat_clusters'))
  dev.off()
  #phate
  phateout[[q]]=phate(t(sctout[[q]]@assays$SCT@counts),ndim=50,n.jobs=10,seed=804L)  
  phatez = cluster_phate(phateout[[q]],k=length(pbmc@meta.data$seurat_clusters))
  saveRDS(phateout[[q]],paste0(INT_FILE_LOC,"/phate_",stage[q],"_postSCT.csv"))
  fig1=plot_ly(as.data.frame(phateout[[q]]), mode="markers", x = ~PHATE1, y = ~PHATE2, color =pbmc@meta.data$batch, colors = viridis(n=2),
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1),
               text = ~paste('<br>PHATE1',PHATE1,'<br>PHATE2', PHATE2))
  fig2=plot_ly(as.data.frame(phateout[[q]]), mode="markers", x = ~PHATE1, y = ~PHATE3, color =pbmc@meta.data$batch, colors = viridis(n=2),
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1),
               text = ~paste('<br>PHATE1',PHATE1,'<br>PHATE2', PHATE3))
  fig3=plot_ly(as.data.frame(phateout[[q]]), mode="markers", x = ~PHATE2, y = ~PHATE3, color =pbmc@meta.data$batch, colors = viridis(n=2),
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1),
               text = ~paste('<br>PHATE1',PHATE2,'<br>PHATE2', PHATE3))
  htmlwidgets::saveWidget(as_widget(fig1), paste0(OUPUT_FILE_LOC,"/QC/",stage[q],"_phate_clustering_post_SCT_1_2.html"))
  htmlwidgets::saveWidget(as_widget(fig2), paste0(OUPUT_FILE_LOC,"/QC/",stage[q],"_phate_clustering_post_SCT_1_3.html"))
  htmlwidgets::saveWidget(as_widget(fig3), paste0(OUPUT_FILE_LOC,"/QC/",stage[q],"_phate_clustering_post_SCT_2_3.html"))
  
  #kbet to re-examine what effect of SCT batch regression is
  subset_size <- 0.1 #subsample to 10% of the data
  subset_id <- sample.int(n = length(pbmc@meta.data$batch), size = floor(subset_size * length(pbmc@meta.data$batch)), replace=FALSE)
  batch.estimate <- kBET(t(as.matrix(pbmc@assays$SCT@counts))[subset_id,], pbmc@meta.data$batch[subset_id])
  saveRDS(batch.estimate,paste0(INT_FILE_LOC,"/kbet_",stage[q],"_postSCT.csv"))
  pdf(paste0(OUPUT_FILE_LOC,"/QC/",stage[q],"_kbet_post_SCT.pdf"))
  plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                    each=length(batch.estimate$stats$kBET.observed)), 
                          data =  c(batch.estimate$stats$kBET.observed,
                                    batch.estimate$stats$kBET.expected))
  try(print(g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
              labs(x='Test', y='Rejection rate',title='kBET test results') +
              theme_bw() +  
              scale_y_continuous(limits=c(0,1))))
  dev.off()
}

#ONLY THE NEONATAL STAGE DISPLAYED SIGNIFICANT BATCH EFFECT
#SO WE EXPLORED THE USE OF HARMONY TO CLUSTER THE NEONATAL STAGE BUT LATER DID NOT USE THIS AS OUR DIMENSIONALITY REDUCTION FOR DOWNSTREAM ANALYSIS.
#As detailed in the methods of the main paper, we acknowledge the batch effect present in this dataset, but we agree with the literature that methods to handle this
#are not yet able to separate out technical and biological noise, and have some a priori concerns around the use of high assumption techniques to distinguish
#when they have been shown to introduce increaed risk of type I failures downstream
#However, we present all our work here to demonstrate our early interaction with the dataset to establish its quality and quirks.

sctoutumi = list()
for(t in 1:6){
  pbmc = readRDS(paste0(INT_FILE_LOC,"/pbmc_",stage[t],"_postfilt.csv"))
  if(t!=6){
    sctoutumi[[t]] <- SCTransform(pbmc, vars.to.regress = c("nUMI"), verbose = TRUE, min_cells=2, return.only.var.genes=FALSE,do.scale=TRUE)
  }
  if(t==6){
    sctoutumi[[t]] <- SCTransform(pbmc, vars.to.regress = c("nUMI"), verbose = TRUE, min_cells=2, return.only.var.genes=FALSE,do.scale=FALSE)
  }
  saveRDS(sctoutumi[[t]],paste0(INT_FILE_LOC,"/SCTumi_",stage[t],".csv"))
  if(t==6){
    batchy = as.data.frame(sctoutumi[[t]]@meta.data$batch)
    colnames(batchy)=c('batch')
    harmony_embed=RunHarmony(sctoutumi[[t]],verbose=TRUE,plot_convergence=TRUE,"batch")
    saveRDS(harmony_embed,paste0(INT_FILE_LOC,"/Harmony_",stage[t],".csv"))
    options(repr.plot.height = 5, repr.plot.width = 12)
    p1 <- DimPlot(object = harmony_embed, reduction = "harmony", pt.size = .1, group.by = "batch")
    p2 <- VlnPlot(object = harmony_embed, features = "harmony_1", group.by = "batch", pt.size = .1)
    
    pdf(paste0(OUPUT_FILE_LOC,"/QC/Neonatal_harmony.pdf"))
    print(p1)
    print(p2)
    dev.off()
  }
  
}

####Harmony can't be used for kBEt as gives PCA out, but is minimally 'invasive' and only affects the dim reduction not the underlying counts.
harmony_embed=readRDS(paste0(INT_FILE_LOC,"/Harmony_",stage[t],".csv"))
clustout = list()
for(u in 1:6){
  if(u !=6){
    pbmc=readRDS(paste0(INT_FILE_LOC,"/SCTumi_",stage[u],".csv"))
    DefaultAssay(pbmc) = 'SCT'
    pbmc=RunPCA(pbmc,dims=1:35)
    pbmc = RunUMAP(pbmc,reduction="pca",dims=1:30)
    pbmc=FindNeighbors(pbmc,reduction='umap',dims=1:2)
    pbmc = FindClusters(pbmc,resolution=0.5)
    saveRDS(pbmc,paste0(INT_FILE_LOC,"/",stage[u],"res0.5_clustered.csv"))
    clustout[[u]] = pbmc
  }
  if(u==6){
    pbmc=readRDS(paste0(INT_FILE_LOC,"/Harmony_",stage[u],".csv"))
    DefaultAssay(pbmc) = 'SCT'
    pbmc = RunUMAP(pbmc,reduction="harmony",dims=1:30)
    pbmc=FindNeighbors(pbmc,reduction='harmony')
    pbmc = FindClusters(pbmc,graph.name = 'RNA_snn',resolution=0.5)
    saveRDS(pbmc,paste0(INT_FILE_LOC,"/",stage[u],"res0.5_clustered.csv"))
    clustout[[u]] = pbmc
  }
  pdf(paste0(OUTPUT_FILE_LOC,"/Pre_lim_clust/",stage[u],"res0.5_clustered.pdf"))
  print(DimPlot(pbmc,group.by='seurat_clusters'))
  dev.off()
}

for(s in 1:6){
  pbmc=readRDS(paste0(INT_FILE_LOC,"/",stage[s],"res0.5_clustered.csv"))
  pbmc.markers=FindAllMarkers(pbmc,only.pos=TRUE, min.pct=0.1,verbose=TRUE,assay='SCT',logfc.threshold=0.1)
  saveRDS(pbmc.markers,paste0(INT_FILE_LOC,"/",stage[s],"res0.5_clustermarkers.csv"))
  print(stage[s])
  options(dplyr.print_max = 1e9)
  pbmc.markers.group = pbmc.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
  pdf(paste0(OUTPUT_FILE_LOC,"/Pre_lim_clust/",stage[s],"res0.5_clustered_top.pdf"))
  print(DimPlot(pbmc,group.by='seurat_clusters'))
  for(i in pbmc.markers.group$gene){
    print(FeaturePlot(pbmc,features=c(i)))  
  }
  dev.off()
  
}

######
#The following represents some preliminary clustering on the basis of the above
clusterlists = list()
clusterlists[[3]]=c("eVM","eVM","eVM","eVM","CM","EC","CM","FB","eVM","EP","eVM","FB","eAM","EC","FB","EC","CM","SM","eVM","eAM","EC","EC","CCS","SM","CCS","SM","EC","MC","Platelets","Platelets")
clusterlists[[4]]=c("VM","VM","VM","EC","FB","VM","FB","VM","VM","VM","VM","FB","EC","VM","VM","AM","CCS","VM","VM","AM","EC","EP","MC","CCS","EC","FB","SM/pericyte")
clusterlists[[5]]=c("VM","VM","VM","VM","VM","VM","FB","EC","VM-CCS","VM","FB","EC","AM-CCS","AM","VM","VM","VM","VM","FB","FB","VM","VM","VM","EC","MC","EP","AM","EC","FB","EC")
clusterlists[[6]]=c("AM","AM","AM","AM","FB","VM","VM","AM-CCS","EC","VM","FB","AM-CCS","VM-CCS","SM","Testis","VM","MC","EP")


for(k in 1:6){
  pbmc=readRDS(paste0(INT_FILE_LOC,"/",stage[k],"res0.5_clustered.csv"))
  new.cluster.ids <- clusterlists[[k]]
  names(new.cluster.ids)=levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.cluster.ids)
  png(paste0(OUTPUT_FILE_LOC,"/Pre_lim_clustt/",stage[s],"res0.5_named_types.png"))
  DimPlot(pbmc,split.by='seurat_clusters')
  DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  dev.off()
  saveRDS(pbmc,paste0(INT_FILE_LOC,"/",stage[k],"res0.5_named_types.csv"))
}
pbmcout=list()
for(i in 1:6){
  pbmcout[[6]] = readRDS(paste0(INT_FILE_LOC,"/",stage[i],"res0.5_clustered.csv"))
}


#IDENTIFYING A PARTICULAR SUBCLUSTER IN NEONATAL STAGE THAT CONTAINED MARKERS OF CELL DAMAGE - AND SUBCLUSTERING AT HIGHER 'RESOLUTION'
output = pbmcout[[6]]
output=RunPCA(output,dims=1:100)
DimHeatmap(output, dims = 21:40, cells = 500, balanced = TRUE)
output = RunUMAP(output,reduction="pca",dims=1:28)

output=FindNeighbors(output,reduction='pca',dims=1:23)
output = FindClusters(output,resolution=0.9)

#####
Harmonyneo=RunHarmony(pbmcout[[6]],verbose=TRUE,plot_convergence=TRUE,"batch",assay.use='SCT')

DimHeatmap(Harmonyneo, dims = 1:5, cells = 500, balanced = TRUE,reduction='harmony')
Harmonyneo = RunUMAP(Harmonyneo,reduction="harmony",dims=1:28)
Harmonyneo=FindNeighbors(Harmonyneo,reduction='harmony')
Harmonyneo = FindClusters(Harmonyneo,graph.name = 'RNA_snn',resolution=1.8)


saveRDS(Harmonyneo,paste0(INT_FILE_LOC,"/output_harmony_res1_8.csv"))

DimPlot(Harmonyneo,reduction = "umap",group.by='seurat_clusters', label = TRUE, pt.size = 0.5,raster=TRUE) + NoLegend()

output.markers <- FindAllMarkers(Harmonyneo, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.10)
saveRDS(output.markers,paste0(INT_FILE_LOC,"/harmony_MARKERS_Neo.csv"))
View(output.markers[which(output.markers$avg_log2FC>0&output.markers$p_val_adj<0.05),])



DimPlot(output,group.by='seurat_clusters',reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(output,reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(Harmonyneo,features=c('Prm1','Prm2','Meig1','Dnajb3','Gsg1'))


###########################saving stage by stage
E85=c("Neuronal","Neuronal","Mixed cells replicating","Mesoderm","Cardiogenic Mesoderm","Neuronal","Neuronal","Mesoderm","Endoderm","Endoderm","Neuromesoderm","NC","Neuromesoderm","Neuronal","Neuronal","EC","Cardiogenic Mesoderm","CM","Endoderm","EC","Neuronal","CM","Endoderm","Notochord","Erythrocytes","Endoderm","Stem cell")
E105=c("Mixed","Mixed","Neural","Mesoderm","Endoderm-derived tissue","Neural","Mixed","Neural","Neural","Cardiogenic mesoderm","Mesoderm","Cardiogenic mesoderm","Cardiogenic mesoderm","Mesoderm","Mesoderm","Neural","EC","Neural","Neural","Neural Crest","EC","Spinal cord","Cardiomyocytes","Neural","Erythrocytes","Neuromesoderm","Mesoderm","Skeletal muscle","Immune cells","Neural","Developing gut","Erythrocytes")
E125=c("Early VM","Early VM","Cardiomyocyte","EC cells","FB","Cardiomyocyte","Smooth muscle","Early VM","Early VM","FB","Atrioventricular Node","EP","AM","AM","Mixed mesoderm-derived","EC cell","Cardiomyocyte","EC cells","Epicardial","EC/Immune cells","Platelets")
E145=c("VM","Early VM","VM","FB","VM","FB","VM","Early VM","EC","Smooth muscle","AM","Early AM","Atrioventricular canal","EC","Trabecular","EC","EP","Immune cells","Neural Crest")
E165=c("VM","VM","VM","VM","VM","FB","VM","EC","FB","AM","VM","Trabecular VM","AM","VM","EC","Pericyte","FB","Nodal","Immune cells","EP","VM","FB","EC")
Neo1.8=c("AM","AM","AM","AM","AM","AM","AM","AM","AM","AM","AM","AM","FB","VM","AM","AM","AM","VM","AM","EC","FB","AM","FB","SM_Pericytes","AM","VM","FB","AM","VM","AM/VM","EP","Immune cells","Germ cell","EC","Germ cell","FB","EC")
clustersnamed=list(E85,E105,E125,E145,E165,Neo1.8)
######################################

#The above is left in for completeness and transparency, such that one can see our attempt to quantify batch effect
#And the impact of clustering on a stage-wise basis using Harmony-adjusted PCA dimension reductions on a batch-effect affected stage.
#The above clusters are preliminary and are not used downstream, nor are any other attempts to reduce mitigate batch effect.
#As below, we observed that when we collated the whole dataset together, the batch effect appears to have little impact when viewed at a whole-dataset level




#### Final whole processed and filtered scRNAseq dataset construction

######################AT THIS STAGE WE PRODUCE A NEW OBJECT COMBINING ALL THE STAGES AND USE THIS FOR DOWNSTREAM CLUSTERING AND ANALYSIS
#warning - this creates a large object.
pbmc.big = merge(pbmcout[[1]],y=c(pbmcout[[2]],pbmcout[[3]],pbmcout[[4]],pbmcout[[5]],pbmcout[[6]]),add.cell.ids=c("E8.5","E10.5","E12.5","E14.5","E16.5","Neo"),project="scRNASeqxST")
saveRDS(pbmc.big,paste0(INT_FILE_LOC,"/merged_res0.5_named_types.csv"))

pbmc.big=SCTransform(pbmc.big, assay='RNA',new.assay.name='SCT',ncells=16000,variable.features.n=9000,return.only.var.genes=FALSE,seed.use=804L, vars.to.regress = c('nCount_RNA'), verbose = TRUE, min_cells=3,do.scale=TRUE)
pbmc=RunPCA(pbmc.big,dims=1:100)
DimHeatmap(pbmc, dims = 21:40, cells = 500, balanced = TRUE)
#take 36 dims - have a look at the QC heatmaps and this is where signals really degrade biologically
pbmc = RunUMAP(pbmc,reduction="pca",dims=1:36,n.components=3L)
pbmc=FindNeighbors(pbmc,reduction='pca',dims=1:36)
pbmc = FindClusters(pbmc,resolution=2)

#For clustering label assignment please see code segment just below.

DimPlot(pbmc,reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#Following cluster classification, the following object was generated.
#This file is the dataset used for the whole filtered scRNA datatset analysis downstream e.g. for Supplementary Figures 2 & 3 etc.
saveRDS(pbmc,paste0(INT_FILE_LOC,"/merged_all_res2_clustered_dim36.csv"))


pbmc.markers=FindAllMarkers(pbmc,only.pos=FALSE, min.pct=0.15,verbose=TRUE,assay='SCT',logfc.threshold=0.1)
saveRDS(pbmc.markers,paste0(INT_FILE_LOC,"/merged_all_markers_res2_clustered_dim36.csv"))

#plot 3D UMAP to examine clusterings
figclust=plot_ly(as.data.frame(pbmc@reductions$umap@cell.embeddings), type="scatter3d", mode="markers", x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color =pbmc@meta.data$updated7 , colors = viridis((length(unique(pbmc@meta.data$updated7)))),
                 marker = list(symbol = 'circle', sizemode = 'diameter'),alpha=1, sizes = c(1),
                 text = ~paste(pbmc@meta.data$updated7,'<br>UMAP_1',UMAP_1,'<br>UMAP_2', UMAP_2, '<br>UMAP_3',UMAP_3))
htmlwidgets::saveWidget(as_widget(figclust), paste0(OUTPUT_FILE_LOC,"/plotlys/merged_all_res2_dim36_3d_newclusters.html"))

########## The following assignment of labels to clusters, and any subclustering - is representative, not exhaustive, as specific code was re-used for different subclustering needs
#The sequential clustering and reclustering are stored in the final objects' meta.data from updated - > updated8. Most sub-clustering took place at the stage of update6/7.
#Primarily this involved a range of different visualisation processes looking at gene expression across clusters.
#Such as FeaturePlot(), ViolinPlot(), and some use of ComplexHeatmap to examine different clusters.

#Clusters were subclustered as appropriate and necessary using the below as a worked example of sub_clustering and clustering.
#Similar code was used for clustering, entered into the console but this was a several month process of iterative clustering and subclustering.
#Thus, for presentability, we show this representative example with the already clustered pbmc object loaded in here: 

############################ SUBCLUSTERING HERE INVOLVED CHANGING A RANGE OF CLUSTERS BEING INVESTIGATED AND AND ITERATIVE PROCESS.
#EARLY CARDIOMYOCYTES IS REPRESENTATIVE HERE AS MANY CLUSTERS WENT THROUGH THIS PROCESS OF SUBCLUSTERING

#Here we include a range of genesets that were used to distinguish between cell types:

Pijuang=unique(c('Dnmt3b','Pou5f1','Epcam','Utf1','Nanog','Eomes','Nkx1-2','Ifitm3','Dnd1','Dppa3','Hhex','Otx2','Mixl1','Gsc','Lhx1','Fgf8','Pax7','Foxa2','T','Chrd','Noto','Cer1','Apela','Mesp1','Lefty2','Mesp2','Osr1','Cdx2','Hes7','Hoxb1','Cdx1','Gbx2','Tbx1','Meox1','Tcf15','Aldh1a2','Dll1','Tbx6','Isl1','Tcf21','Smarcd3','Acta2','Tagln','Myl7','Nkx2-5','Myl4','Tnnt2','Hoxa10','Hoxa11','Tbx4','Bmp4','Ahnak','Pmp22','Krt18','Krt8','Igf2','Etv2','Kdr','Anxa5','Pecam1','Cdh5','Lmo2','Runx1','Cited4','Gata1','Hbb-bh1','Gypa','Hba-a1','Hba-a2','Hes3','Epha5','Cdx4','Hoxb9','Sox2','Irx3','Hesx1','Six3','Sox9','Pax3','Tfap2a','Foxd3','Dlx2','Sox10','Wnt1','En1','Pax2','Pax6','Foxg1','Trp63','Grhl2','Grhl3','Dkk1','Krt19','Amot','Spink1','Emb','Cystm1','Apoe','Apoa2','Ttr','Tfap2c','Ascl2','Elf5','Sparc','Plat','Lamb1'))
mesodermgenes=unique(c('Rhox5','Fabp3','Tfap2a','Plm2','Otx2','Lhx5','Hes3','Hoxaas3','Wnt3a','Nkx1-2','T','Tbx6','Gal','Dll1','Aldh1a2','Pcdh19','Meox1','Mesp1','Cdx4','Foxf1','Tcf15','Foxc2','Foxc1','Tbx1','Bik','Mef2c','Irx5','Jag1','Alx1','Cer1','Ptgfr','AV02606B','Sfrp5','Rgs5','Tbx5','Tdo2','Acta2','Pitx1','Hoxd9','Hoxc8','1110017D15Rik','Noto','Gpx2','Foxa1','Sox17','Kdr','Fll1','Flt1','Klf1','Cited4','Gata1','Hs3st1','Apoa2','Eif4a1','Eif5a','Bmp2','Bmp4','Wnt2','Cyp26a1','Fgf8','Fgf10','Dkk1','Raldh2','Pax3','Pax2','Sox2','Sox9','Sox10','Onecut','Stmn1','Hoxb9','Epha5','Ptn','Tubb5','En1','Noto','Pyy','Cldn8','Alcam','Epcam','Ttr','Prrx2','Pgf','Hoxa11','Ntrk2','Robo1','Slit2','Hells','Stmn2'))
hfgenes=unique(c('Nebl','Acta2','Myocd','Pmp22','Tnnt2','Myl4','Krt18','Krt8','Mest','Hand1','H19','Car4','Myl7','Rgs5','Nkx2-5','Cdkn1c','Csrp2','Id2','Peg10','Mef2c','Cfc1','Tpm1','Dkk1','Id3','Apoe','Alx1','Id1','Wls','Fgf8','Sfrp1','Cyp26a1','Bik','Tbx1','Fst','Irx5','Irx3','RP23-220F20.2','Cer1','Ppp1r1a','Socs2','Tcf15','Pcdh19','Meox1','Pou5f1','Pdlim4','Trp53i11','Rbp1','Nnat','Pcdh8','Ripply2','Dll1','Raldh2','Rspo3','Meis2','Lef1','Tnnt2','Pmp22','Hand1','Tbx5','Tbx1','Tbx18','Tbx20','Shox2','Rhoa','Isl1','Bmp10','Hcn4','Nkx2-5','Hand2','Mef2c','Mef2a','Fgf8','Fgf10','Pdpn','Dkk1','Bik','Id2','Jag1','Wnt2','Meox1','Ripply2','Gata6','Gata4','Gata5','Tbx3','Shox2','Hcn4','Pitx2','Pitx1','Evx2','Nodal','Lefty1','Foxa2','Foxh1','Gsc','Hoxb1','Tbx6','Raldh2','Osr1','Foxf1','Hoxb6','Fgf4','Fgf8','Shh'))
ccsgenes=unique(c('Myl2','Myl7','Myh6','Myh7','Nppa','Nppb','Smpx','Bmp10','Hey2','Hopx','Tbx5','Tbx3','Notch3','Jag1','Isl1','Shox2','Smoc2','Tbx18','Ephb3','Tbx2','Bmp2','Bmp4','Rspo3','Gata5','Gata4','Ntm','Etv1','Sema3a','Sema3c','Kcnj3','Kcnj5','Hcn4','Hcn2','Scn5a','Scn10a','Cacna2d2','Cacna1g','Cacna1h','Cacna1l','Gja1','Gja5','Gjc1','Gjd3','Gjb6','Igfbp5','Cpne5','Pdgfra','Prrx1','Smarca4','Myoz2','Cntn2','Kcne1','Hcn2','Hcn1','Kcnh2','Wnt5a','Gja3','Prox1','Bmpr1a','Ncam-1','Alcam','Scn3b','Kcnn4','Kcnd2','Kcnk2','Kcne1','Cacnb1','Nkx2-5','Dkk3','Rgs6','Id3','Bves','Bex4','Irx3','Irx2','Irx5','Dbh','Slc22a1'))
mesengenes=unique(c('Ecscr','Eng','Emcn','Tie1','Icam2','Cdh5','Tek','Pecam1','Cav1','Kdr','Cdh11','Ntrk2','Nfatc1','Vwf','Sox9','Ptn','Col1a1','Col1a2','Dcn','Ogn','Col5a1','Fbn1','Eln','Postn','Col3a1','Fbln5','Fbln2','Ncam1','Fbln1','Vim','Fn1','Fstl1','Gsn','Sparc','Mmp2','Cald1','Acta2','Cnn1','Myh11','Sfrp2','Cpz','Pdgfrl','Tagln','Tagln2','Kcnj8','Neb','Tnnt2','Myh3','Tpm2','Myod1','Myog','Myf5','Ryr1','Cthrc1','Pdgfra','Postn','Sox9','Fbln2','Apoe','Lef1','Csf1r','Cd163','Csar1','Fcgr1','Cd2','Cd69','Aldh1a2','Upk1b','Upk3b','Kcne1l'))

dbhcmz <- list("Mesoderm", "Early CMs", "Early VMs", "Early AMs", "VMs", "AMs")
stage <- c("E8_5","E10_5","E12_5","E14_5","E16_5","Neo")
Early <- list('Mesp1','Kdr','Hand2','Maoa','Maob','Comt','Th','Pah','Ddc')
CMs <- list('Pitx2','Tbx20','Cryab','Calc1','Foxh1','Hand1')
CMs2 <- list('Nppa','Slc18a2','Slc6a3','Pnmt','Slc22a1','Gata4','Fgf10','Fgf8','Mef2c','Foxc1','Foxc2','Six2','Tbx1','Tbx5','Osr1','Wnt2','Sfrp5')
Sarcoplasmic <- list('Ryr2','Ryr3', "Casq2" , 'Pln')
Myofibrillar_1 <-list('Ctnna3', 'Tnnt2' ,'Tnnc1','Tnni3')
Myofibrillar_2 <- list('Nebl','Myocd','Ttn','Nexn')
Myofibrillar_3 <-list( 'Myl7','Myh6', 'Myl4','Myl2')
Sarcolemmal <- list('Hcn4','Cacna1c','Actc1','Smpx')
MC1 <- list('Cthrc1', 'Pdgfra', 'Postn', 'Sox9')
MC2 <- list( 'Fbln2', 'Apoe', 'Lef1', 'Tagln')
EP <-  list('Aldh1a2', 'Upk3b', 'Kcne1l', 'Tmem255a')
EC1 <- list('Pecam1', 'Emcn', 'Tie1', 'Kdr')
EC2 <- list('Cdh5', 'Ecscr','Ptn','Pla2g7')
Blood <- list('Runx1', 'Hbb-y', 'Hba-x', 'Runx2')
NC1 <- list('Snai1', 'Foxd3', 'Snai2', 'Sox10')
NC2 <- list('Sox9', 'Twist1', 'Dbh', 'Pnmt')
FB1 <- list('Col1a2','Col3a1','Col1a1','Pde1a')
FB2 <- list('Col5a1','Thbs1','Dcn','Tcf21')
SMC1 <- list('Neb', 'Myh3', 'Cald1', 'Tpm2')
SMC2 <- list('Myod1','Myog','Myf5','Ryr1')
Atrial <- list('Nr2f2','Cav1','Vsnl1','Stard10','Phkg1','Crabp2','Slc24a5','Myl7','Myl4','Myh6','Sln','Kcnj3','Kcnj5','Epha4')
Ventricular <- list('Myl2','Pln','Mpped2','Sall4','Irx4','Hey2','Gja1','Myl3','Nppb','Scn5a','Myh7')
RA <- list('Shox2','Adm','Tigd2','Myl7')
LA <- list('Pitx2','Klhl41','Prss23', 'Myl7', 'Meis2')
AVC <- list('Rspo3','Tbx3','Bmp2','Smpx','Tbx2','Id2','Gata6','Ncx1','Gjb6','Gjc1')
OFT <- list('Itm2a','Gli3','Rspo3','Cxcl12')
PTO <-list('Homer2','Myl2','Tnnt2','Tnnc1')
DTO <- list('Bmp4','Col1a2','Isl1','Tnnt2')
LV <- list('Gja1','Malat1','Lyrm4','Hand1','Tbx5')
RV <- list('Gata4','Mcam','H19','Hspb7','Mfapb2')
Septum <- list('Tgfbr3','Pcsk6','Hs6st2','Tnnt2','Irx1','Irx2')
Trabec <- list('Bmp10', 'Myl2', 'Tnnt2', 'Myh6')
PKJ <- list('Irx3','Tbx18','Hcn4','Casq2','Gja5', 'Sema3a','Cntn2','Nkx2-5','Id','Tbx18','Hopx')

Depol_Repol <- c("Scn5a","Kcnd2","Kcnd3","Kcna4","Kcna7","Kcnc4","Kcna5","Kcnc1","Kcnh2","Kcnq1","Kcnj2","Kcnj12","Kcnj11","Kcnj3","Kcnj5","Kcnk1","Kcnk6","Kcnk3","Kcnk4","Hcn2","Hcn4")
Conduction <- c("Gja1","Gjc1","Gja5",'Gjb6',"Gja4",'Sema3a','Cntn2','Nkx2-5','Id','Tbx18','Hopx')
Calcium_handling <- c("Ryr1","Ryr2","Ryr3","Itpr3","Itpr2","Itpr1","Atp2b1","Atp2b2","Atp2b3","Atp2b4","Atp2a1","Atp2a2","Atp2a3","Mcu","Slc8a1","Slc8a2","Slc8a3","Calr","Cas1q","Casq2","Trpv1","Cacna1s","Cacna1c","Cacna1d","Cacna1f","Cacna1a","Cacna1b","Cacna1e","Cacna1g","Cacna1h","Cacna1l","Camkd2","Camkd1")
Contraction <- c("Ttn","Acta1","Actc1","Tpm1",'Nebl','Myocd','Nexn','Ctnna3', 'Tnnt2' ,'Tnnc1','Tnni3',"Myl7","Myl2","Myl1","Myl3","Myl4","Myh6","Myh7","Mybpc3")
Secretion <- c('Nr2f2','Nr2f1',"Chga","Chgb","Chgc","Scg3","Scg4","Cpe","Pc1","Pc2","Plat","Atp6v0a1","Atp6v0a2","Atp6v0a4","Atp6v0b","Atp6v0c","Atp6v0d1","Atp6v0d2","Atp6v0e1","Atp6v0e2","Slc18a1","Slc18a2","Npy","Adm","Adm1","Adm2","Adcyap1","Nppa","Nppb","Nppc")
Catecholamines <- c("Th","Ddc","Dbh","Pnmt","Slc22a1",'Slc18a2','Slc18a1','Comt','Maoa','Maob')
Neural_crest <- c("Sox10","Foxd3",'Sox9', 'Twist1','Snai1','Snai2')
Developmental <- c('Isl1','Bmp10','Bmp4','Pitx2','Shox2','Hey2','Hey1','Gata4','Gata6','Tbx1','Tbx2','Tbx3','Tbx5','Tbx18','Tbx20','Fgf10','Fgf8','Mef2c','Foxc1','Foxc2','Six2','Osr1','Wnt2','Sfrp5','Irx1','Irx2','Irx3','Irx4','Nkx2_5')
Stem_cell <- c('Ly6a','Kdr','Oct4','Sox2','Gata4','Otx2','Snai1','Foxa2','Pdx1','Vegfr2','Sox17')
Neuronal <- c('Th','Chat','Npy','Tac1','Uchl1')
Wnt <- c('Dvl1','Wnt1','Wnt2','Fhl2','Dkk2','Dkk3','Vangl2','Rock2')
Metabolic <- c('Tnfaip3','Tieg','Tieg2','Ptgfrn')
Mesodermal <- c('Fgf19','Nodal','Lefty1','Calca','T')
Progen = c('Snai1','Dll1','Tacr3','Mixl1','Dkk1')


genetypez <- list(Depol_Repol,Conduction,Calcium_handling,Contraction,Secretion,Catecholamines,Neural_crest,Developmental,Stem_cell,Neuronal,Wnt,Metabolic,Mesodermal,Progen)
Types <- list(CMs, CMs2, Sarcoplasmic, Myofibrillar_1, Myofibrillar_2, Myofibrillar_3, Sarcolemmal,NC1, NC2, Atrial, Ventricular, RA, LA, AVC, OFT, PTO, DTO, LV, Septum, Trabec, PKJ, Early, MC1, MC2, EP, EC1, EC2, Blood)
listo <- unlist(Types)
genezz <- unlist(genetypez)
goi <- unique(listo,genezz)

#########Identify subclusters
sub_39=FindSubCluster(pbmc,cluster='Early Cardiomyocytes',graph.name='SCT_snn',resolution=0.5,subcluster.name='neweCM')
sub_39@meta.data$updated = as.character(sub_39@meta.data$neweCM)
sub_39@meta.data$transfer=as.factor(sub_39@meta.data$updated)
#Demonstrates that this separated Early-CMs into subclusters from 0 -> 7
unique(sub_39@meta.data$transfer)

#pbmc@meta.data$updated6=sub_39@meta.data$transfer
sub_sub=subset(sub_39,subset=updated7=='Early Cardiomyocytes')

#THIS WAS DONE ON A SUBCLUSTER-WISE BASIS TO ASSIGN THE LABELS
#sub_39@meta.data$updated[which(sub_39@meta.data$updated=='Early Cardiomyocytes_7')]='Early Cardiomyocytes'

#DIFFERENT GENES WERE SUBSTITUTED IN HERE FOR EXAMINING DIFFERENT SUBCLLUSTERS
VlnPlot(sub_sub,group.by='neweCM',assay='SCT',features=c(intersect(unique(c('Kdr','Col1a1','Col1a2','Dcn','Ogn','Eln','Fbn','Cnn1','Acta2','Actc1','Ttn','Tbx5','Tbx20','Bmp10','Myh6','Myh7','Myl2','Myl7','Col1a1','Col1a2','Dcn','Myl7','Tnni3','Actc1','Ryr2','Nebl','Ttn','Hey2','Nr2f2','Myh7','Myh6')
),rownames(sub_sub[['SCT']]@counts))),cols=viridis(150,option="B"),slot='counts',stack=TRUE,flip=TRUE)+NoLegend()

heatclust=list()
n=0
for(i in unique(sub_sub@meta.data$neweCM)){
  n=n+1
  cells=subset(sub_sub,subset=neweCM==i)
  print(i)
  #ha=HeatmapAnnotation(Clust=cells@meta.data$cellidents,col=list(Clust=colorRamp2(c(1,20),c('maroon2','turquoise2'))),show_legend = FALSE,show_annotation_name = FALSE)
  ha=HeatmapAnnotation(Stage=cells@meta.data$orig.ident,col=list(Stage=c('E8'='blue4','E10'='cadetblue','E12'='blue3','E14'='royalblue1','E16'='brown2','Neo'='deeppink')),show_legend = FALSE,show_annotation_name = FALSE)
  heatclust[[as.numeric(n)]]=Heatmap(cells@assays$SCT@scale.data[c(intersect(unique(c(mesengenes)),rownames(cells@assays$SCT@scale.data))),],column_title=paste0('Cluster: ',i),column_title_rot = 90,top_annotation = ha,cluster_rows = FALSE,show_column_names = FALSE,show_heatmap_legend = FALSE,heatmap_legend_param = list(at=c(-3,0,1,2,5,10)),col=colorRamp2(c(-3,0,1,5,10),c('navyblue','grey','darkorange1','red1','red2')))
}

#SOME SUB-SUB-CLUSTERING WAS REQUIRED ON CERTAIN SUB-CLUSTERS AS BELOW
sub_sub=SCTransform(sub_sub)
sub_sub=RunPCA(sub_sub,dims=1:100)
DimHeatmap(sub_sub, dims = 21:40, cells = 500, balanced = TRUE)

sub_sub = RunUMAP(sub_sub,reduction="pca",dims=1:39)
sub_sub=FindNeighbors(sub_sub,reduction='pca',dims=1:39)
sub_sub = FindClusters(sub_sub,resolution=0.5)

sub_sub.marks=FindAllMarkers(sub_sub,verbose=TRUE)
View(sub_sub.marks[which(sub_sub.marks$avg_log2FC>0&sub_sub.marks$p_val_adj<0.05),])

heatclust=list()
n=0
sub_sub@meta.data$working=sub_sub@active.ident
for(i in unique(sub_sub@active.ident)){
  n=n+1
  cells=subset(sub_sub,subset=working==i)
  print(i)
  #ha=HeatmapAnnotation(Clust=cells@meta.data$cellidents,col=list(Clust=colorRamp2(c(1,20),c('maroon2','turquoise2'))),show_legend = FALSE,show_annotation_name = FALSE)
  ha=HeatmapAnnotation(Stage=cells@meta.data$orig.ident,col=list(Stage=c('E8'='blue4','E10'='cadetblue','E12'='blue3','E14'='royalblue1','E16'='brown2','Neo'='deeppink')),show_legend = FALSE,show_annotation_name = FALSE)
  heatclust[[as.numeric(n)]]=Heatmap(cells@assays$SCT@scale.data[c(intersect(unique(c(mesengenes)),rownames(cells@assays$SCT@scale.data))),],column_title=paste0('Cluster: ',i),column_title_rot = 90,top_annotation = ha,cluster_rows = FALSE,show_column_names = FALSE,show_heatmap_legend = FALSE,heatmap_legend_param = list(at=c(-3,0,1,2,5,10)),col=colorRamp2(c(-3,0,1,5,10),c('navyblue','grey','darkorange1','red1','red2')))
}

lgd = Legend(col_fun=colorRamp2(c(-3,0,1,5,10),c('navyblue','grey','darkorange1','red1','red2')),title='Expression')
png(paste0(OUTPUT_FILE_LOC,"/example_sub_sub_cluster.png"),height=2000,width=2000)
draw(heatclust[[1]]+heatclust[[2]]+heatclust[[3]]+heatclust[[4]]+heatclust[[5]]+heatclust[[6]]+heatclust[[7]]+heatclust[[8]]+heatclust[[9]]+heatclust[[10]])
dev.off()

FeaturePlot(sub_sub,features=c('Col1a1','Col1a2','Dcn','Myl7','Tnni3','Actc1','Ryr2','Nebl','Ttn','Hey2','Nr2f2','Myh7','Myh6'))
DimPlot(sub_sub,reduction = "umap", label = TRUE, pt.size = 0.5,repel=TRUE) + NoLegend()

#Then subcluster/sub-subcluster labels are assigned, and these are integrated back into the whole object.
#Again, many of these were performed using the console integrating sub-clustering and sub-subclustering labels into updated6.

pbmc@meta.data$updated7=as.character(pbmc@meta.data$updated6)

#The final clustering can be checked on a subclustering-basis as per the below code, using just the top expressed genes as a sense check
#For this please use the clustered object: pbmc = readRDS(#TODO)


for(n in c('Early Cardiomyocytes','Early Atrial Cardiomyocytes')){
  clust=pbmc@assays$SCT@counts[,which(pbmc@meta.data$updated8==n)]
  print(n)
  clust=as.matrix(clust)
  clust=clust[order(rowSums2(clust),decreasing=TRUE),]
  ha=HeatmapAnnotation(Stage=c(substr(colnames(clust),0L,3L)),col=list(Stage=c('E8.'='blue4','E10'='cadetblue','E12'='blue3','E14'='magenta4','E16'='maroon2','Neo'='red2')))
  if(n=='Brain/Spinal Cord'){
    n='Brain_Spinal Cord'
  }
  png(paste0(OUTPUT_FILE_LOC,"/",n,"_counts_split.png"),height=2000,width=2000)
  if(n=='Brain_Spinal Cord'){
    n='Brain/Spinal Cord'
  }
  print(Heatmap(clust[1:100,],show_column_names = FALSE,top_annotation = ha,column_split=pbmc@meta.data$updated6[which(pbmc@meta.data$updated8==n)]))
  dev.off()
}

###We then used Seurat FindAllMarkers as a final sense check for the clustering

pbmc.markers=FindAllMarkers(pbmc,only.pos=FALSE, min.pct=0.10,verbose=TRUE,assay='SCT',slot='counts',logfc.threshold=0.2,test.use='bimod')
saveRDS(pbmc.markers,paste0(INT_FILE_LOC,"/cell_types_markers_01_02_counts.csv"))
pbmc.markerz=pbmc.markers
View(pbmc.markerz %>%
       group_by(cluster) %>%
       top_n(n = 10, wt = avg_log2FC))

###########################
#NOW WE USE THE CLASSIFIED OBJECT FROM ABOVE TO EXTRACT ALL THE CELLS THAT CONTRIBUTE TO THE CARDIOMYOCYTE LINEAGE ACROSS ALL STAGES

cmlin=subset(pbmc,subset=updated8=='Atrial Cardiomyocytes'|updated8=='Cardiac Conduction System'|updated8=='Ventricular Cardiomyocytes_Myh6'|updated8=='Ventricular Cardiomyocytes'
             |updated8=='Early Cardiomyocytes'|updated8=='Early VMs'|updated8=='Early Atrial Cardiomyocytes'|updated8=='ECM Ventricular Cardiomyocytes'|updated8=='ECM Atrial Cardiomyocytes'|updated8=='Cardiac Progenitors')

#Now we renormalise to account for changes in sample gene expresssion
cmlin=SCTransform(cmlin,do.correct.umi=TRUE,ncells=10000,seed.use=804L,variable.features.n = 9000,do.scale=TRUE,do.center=TRUE,return.only.var.genes = FALSE,n_genes=9000,min_cells=10)

#Now we use PHATE for dimensionality reduction after using Seurat::DimHeatmap to establish npca 96 as acceptable drop off in biological signal

phatecm=phate(t(cmlin@assays$SCT@scale.data[VariableFeatures(cmlin),]),ndim=3L,mds.solver='smacof',npca=96,t=14,verbose=TRUE,knn=12,n.jobs=6,seed=804L,gamma=0.8)
saveRDS(phatecm,paste0(INT_FILE_LOC,"/CMlineagephateknn12t14npca96gamma0_8.csv"))
#Running MAGIC on approximated PHATE manifold, to aid with visual identification of cluster identity
magedgenes=magic(t(cmlin@assays$SCT@scale.data[VariableFeatures(cmlin),]),knn=12,t=14,n.jobs=3,seed=804L,npca=96)

AM=c('Nppa','Myl7','Sln','Myl4','Myh7','Myl2')
CCS=c('Ankrd1','Cacna2d2','Igfbp5','Shox2')
VM=c('Myl2','mt_Nd2','Ttn','Myl3','Myh7','Myl7','Myl4','Sln','Myh6')
VM_m=c('Myl2','mt_Nd2','Ttn','Myl3','Myh7','Myl7','Myl4','Sln','Myh6')
eCM=c('Myl3','Acta2','Tnnt2','Tnnt1','Atp2a2')
eAM=c('Myl3','Acta2','Tnnt2','Tnnt1','Tpm1','Atp2a2')
eVM=c('Myl2','Hey2','Atp2a2','Tnnt2','Tnnt1','Tpm1')
FB=c('Col1a1','Col1a2','Col3a1','Vim','Fbn1','Fbln2','Postn','Sost','Dcn','Ogn','Eln')
EAM=c('Ttn','Myl2','Myl3','Myh6','Col3a1','Dcn','Col1a1')
EVM=c('Sln','Myl7','Nppa','Col3a1','Pln','Myl3','Nppb','Dcn','Myh7','Ogn')
SKM=c('Myog','Myh11','Myh3','Ryr1','Myod1','Col1a1','Col1a2','Col3a1','Vim','Fbn1','Fbln2','Postn','Sost','Dcn','Ogn','Eln')
SMC=c('Rgs5','Kcnj8','Eln','Acta2','Cald1','Tagln','Tagln2','Actg2','Myog','Myh11','Myh3','Ryr1','Myod1','Col1a1','Col1a2','Col3a1','Vim','Fbn1','Fbln2','Postn','Sost','Dcn','Ogn','Eln')
EC=c('Emcn','Ecscr','Eng','Cdh5','Cav1','Vwf','Upk3b','Upk1b','Aldh1a2','Tmsb4x','Tmsb10')
EP=c('Upk3b','Upk1b','Aldh1a2','Rgs5','Kcnj8','Eln','Acta2','Cald1','Tagln','Tagln2','Actg2','Myog','Myh11','Myh3','Ryr1','Myod1','Col1a1','Col1a2','Col3a1','Vim','Fbn1','Fbln2','Postn','Sost','Dcn','Ogn','Eln','Emcn','Ecscr','Eng','Cdh5','Cav1','Vwf','Upk3b','Upk1b','Aldh1a2')
ED=c('Kdr','Eng','Tmsb4x','Tmsb10','Nfatc1','Etv2','Etsrp','Foxc1','Wnt2','Hand1','Hand2','Isl1','Hcn4','Nkx2-5','Bmp4','Bmp2','Emcn','Ecscr','Eng','Cdh5','Cav1','Vwf','Upk3b','Upk1b','Aldh1a2')
CPC=c('Wnt2','Hand1','Hand2','Isl1','Hcn4','Nkx2-5','Bmp4','Bmp2','Foxd1','Mest','T','Mesp1','Mesp2','Foxc1','Foxc2','Tcf21','Tcf15','Aldh1a2','Tbx1','Meox1','Tmsb10','Tmsb4x')
MS=c('Foxd1','Mest','T','Mesp1','Mesp2','Foxc1','Foxc2','Tcf21','Tcf15','Aldh1a2','Tbx1','Meox1')
Haem=c('Hemgn','Hhex','Tal1','Redrum','Hbb-bh1')
Plat=c('Pf4','Gp9','Gp1bb')
IC=c('Tyrobp','Fcer1g','Lyz2','C1qb')
Mixed=c('Snrpg','Rps28','Ptma')
NMP=c('T','Sox2','Hoxb5os','Hoxaas3')
NC=c('Foxd3','Sox10','Tfap2b','Ncam1','Mcam','Phox2b')
BS_SC=c('Onecut2','Dcx','Neurod4','Neurog1','Neurog2','Otx2','Tubb3','Dlx1','Isl1','Tlx1','Phox2a')
Neural=c('Hes5','Sox2','Pax6')
Endoderm=c('Ttr','Epcam','Alcam','Pyy','Alcam','Acer2','Foxa1','Afp')
FHF=c('Hand2','Tbx5','Tbx20')
AHF=c('Tbx1','Fgf8','Fgf10')
pSHF=c('Hoxb1','Raldh2','Tbx6')

goi2=list(AM,BS_SC,CCS,CPC,eAM,eCM,eVM,EAM,EVM,ED,Endoderm,EC,EP,FB,Haem,IC,MS,Mixed,Neural,NC,NMP,Plat,SKM,SMC,VM,VM_m)

#Plots a range of interactive 3D plots of smoothed gene expression using MAGIC - some are useful and some are less useful.
#Only used for visual cluster classification
for(j in unlist(goi2)){
  if(j %in% colnames(magedgenes$result)){
    figclust=plot_ly(as.data.frame(phatecm), type="scatter3d", mode="markers", x = ~PHATE1, y = ~PHATE2, z = ~PHATE3, color =magedgenes$result[,j] , colors = viridis((length(unique(magedgenes$result[,j])))),
                     marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1),
                     text = ~paste(j,'<br>PHATE1',PHATE1,'<br>PHATE2', PHATE2, '<br>PHATE3',PHATE3))
    htmlwidgets::saveWidget(as_widget(figclust), paste0(OUTPUT_FILE_LOC,"/magic",j,"_npc96_var9k_knn12_t.html"))
  }
}


phatecm=phate(t(cmlin@assays$SCT@scale.data[VariableFeatures(cmlin),]),ndim=60L,mds.solver='smacof',npca=96,t=14,verbose=TRUE,knn=12,n.jobs=2,seed=804L,gamma=0.8)

#So now cluster to 106 clusters
phategenesclust <- pyphate$cluster$kmeans(phatecm$operator,n_clusters=as.integer(106),random_state=804L)
saveRDS(phategenesclust,paste0(INT_FILE_LOC,"/CMlineageclusteredvars9kknn12t14npca96gamma0_8.csv"))

#Again, following cluster classification assignment, and manual reassignment of individual cells where necesssary, the following object was produced.
#This now represents the final clustered CM lineage Seurat object
#Load in the final cluster assignment here
cmlin=readRDS(cmlin,paste0(DATA_FILE_LOC,"/CMlineage.csv"))

#For representative code of the clustering classification assignment, please see the below.

library(RColorBrewer)
myColors <- viridis(106,option='H')
names(myColors) <- levels(as.character(sort(unique(phategenesclust))))
colScale <- scale_colour_manual(name = "Clusters",values = myColors)
library(ggrepel)
clustcolours=data.frame(sort(unique(phategenesclust)),viridis(106,option='H'))


#Plots clusters coloured specifically
gg=ggplot(data.frame(phatecm$embedding,Cluster=as.factor(phategenesclust))) +
  geom_point(aes(PHATE1, PHATE3, color=Cluster)) +
  labs(color=paste0("Clusters"),title=paste0('Cardiomyocyte Lineage'))+
  colScale +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(),axis.text.x = element_blank(),text=element_text(size=20),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position='none')





#This plots how the new clustering has split the original clustering from the bigger object
for(typeing in c(unique(as.character(cmlin@meta.data$updated8)))){
  coloz=phategenesclust[which(cmlin@meta.data$updated8==typeing)]
  png(paste0(OUTPUT_FILE_LOC,"/Clusters_",typeing,"_npc96_var9k_knn12_t14.png"),width=1250,height=850)
  if(typeing=='Brain_Spinal Cord'){
    typeing='Brain/Spinal Cord'
  }
  embedz=data.frame(phatecm$embedding[intersect(rownames(phatecm$embedding),colnames(cmlin@assays$SCT@counts)[which(cmlin@meta.data$updated8==typeing)]),],Clusters=coloz)
  finder=cbind(rownames(embedz),coloz)
  outputlabs=c()
  y=0
  for(i in c(unique(coloz))){
    y=y+1
    for(j in 1:nrow(finder)){
      if(finder[j,2]==i){
        outputlabs[y]=finder[j,1]
      }
    }
  }
  
  
  
  print(ggplot(embedz) +
          geom_point(aes(PHATE1, PHATE3, color=as.factor(coloz))) +
          labs(color=paste0("Clusters"),title=paste0(typeing,' Clusters'))+
          scale_colour_manual(values=c(clustcolours[levels(as.factor(coloz)),2])) +
          geom_label_repel(max.overlaps=20,data=as.data.frame(cbind(embedz[c(outputlabs),],unique(coloz)))[order(unique(coloz)),],aes(PHATE1,PHATE3),label=sort(unique(coloz)),fill=viridis(length(unique(coloz)),option="H"),color='white'
          )+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
                axis.text.x = element_blank(),axis.text.y = element_blank(),text=element_text(size=20),panel.background = element_blank(), legend.position = 'none',axis.line = element_line(colour = "black"))
  )
  dev.off()
}

#Store the clusters within the cmlin seurat object
cmlin@meta.data$phateyclusts=as.character(phategenesclust)
cmlin@meta.data$phateyclusts=as.factor(cmlin@meta.data$phateyclusts)

####### THE BELOW DOCUMENTS THE CLUSTERING AND SUBCLUSTERING OF THE CMLINEAGE THROUGH CLASSIFICATION OF THE PHATE CLUSTERS
#As above, this represents a vast corpus of work, some of which was carried out within the console and not written in this script and thus this section is representative and not exhaustive.
#Again, clustering was performing used a mixture of ViolinPlots, FeaturePlots, and heatmaps as above. We considered both scale.data and counts in our classification.

################This documents how it matches across to the PHATE cluster numbers assigned
newtypes=c('AM','eVM-trab','AM','AHF','VM','eVM','AM','AVN','PKJ','eAM-CCS','VM','VM','VM','VM','pHT','VM-trab','eAM','eVM-trab','VM-trab','AM','VM','AM','EDCM','VM',
           'eVM','AM','eVM','FHF','VM-trab','AM','VM-trab','AM','AM','DevCM','PKJ','VM-trab','PKJ','VM',
           'eVM','AM','eVM','AVN','AM','Apop','eVM','PKJ','AM','AM','eAM','VM','AM',
           'pSHF','VM-trab','EDCM','VM-trab','PKJ','AHF','eVM-trab','AVN','SAN','AM','AM',
           'VM','VM','VM','AM','eVM','AM','AM','pHT','AM','eAM','AM',
           'VM','pSHF','eAM','VM','pHT','AVN','DevCM','VM','PKJ','AM','eVM',
           'EDCM','AM','eVM-trab','AM','VM-trab','EDCM','AM','VM-trab','eVM-trab','AHF','eVM',
           'AM','VM-trab','AM','VM','AM-CCS','AM','AM','eVM-trab','AM','eVM','eAM-CCS')

names(newtypes)=levels(cmlin@meta.data$phateyclusts)
Idents(cmlin)=cmlin@meta.data$phateyclusts
cmlin=RenameIdents(cmlin,newtypes)
#Store the new clusters within the Seurat Obj
cmlin@meta.data$newtypez=cmlin@active.ident

#merged AM-CCS and eAM-CCS as indistinuishable
cmlin@meta.data$newtypez=as.character(cmlin@meta.data$newtypez)
cmlin@meta.data$newtypez[which(cmlin@meta.data$newtypez=='eAM-CCS')]='AM-CCS'
#Check and compare difference between VM-trab and early VM trabecular.
VlnPlot(subset(cmlin,subset=newtypez=='eVM-trab'|newtypez=='VM-trab'),group.by='newtypez',features=c('Ttn','Myl2','Vcan','Hyal2','Mmp2','Adamts1','Hey2','Notch1','Hopx','Bmp10','Myh10'))

#Finalise re-clustering
cmlin@meta.data$newtypez[which(cmlin@meta.data$phateyclusts=='52')]='eVM-trab'
cmlin@meta.data$newtypez[which(cmlin@meta.data$phateyclusts=='35')]='eVM-trab'

Idents(cmlin)=cmlin@meta.data$newtypez




#Can use this to confirm clustering

for(i in unique(cmlin@meta.data$newtypez)){
  print(paste0(i,": ",length(which(newtypes==i))))
  
  #for(cluster in unique(as.character(phategenesclust))){
  #take the stages and then convert into corrected names
  cells=subset(cmlin,subset=newtypez==i)
  
  stogy=substr(colnames(cells[["RNA"]]@counts),0L,4L)
  stogy[which(stogy=='E10.')]='E10.5'
  stogy[which(stogy=='E12.')]='E12.5'
  stogy[which(stogy=='E14.')]='E14.5'
  stogy[which(stogy=='E16.')]='E16.5'
  stogy[which(stogy=='Neo_')]='P3'
  
  #create the heatmap annotation for stages
  
  ha=HeatmapAnnotation(Stage=stogy,col=list(Stage=c('E8.5'=cividis(6)[1],'E10.5'=cividis(6)[2], 'E12.5'=cividis(6)[3],'E14.5'=cividis(6)[4],'E16.5'=cividis(6)[5],'P3'=cividis(6)[6])),show_legend = FALSE,show_annotation_name = FALSE)
  
  ##SCALE PLOT
  #create ordered data 
  ord=cells@assays$SCT@scale.data
  clust=as.matrix(ord)
  clust=clust[order(rowSums2(clust),decreasing=TRUE),]
  #clust=clust[intersect(earlygenes,rownames(clust)),]
  clust=clust[1:100,]
  col_fun = colorRamp2(c(-10, 0, 10),c(plasma(20)[1], plasma(20)[10], plasma(20)[20]))
  lgd = Legend(col_fun = col_fun, title = "Expression",at = c(-10,0, 10), 
               labels = c("Low",'Zero', "High"),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))
  
  lgd2 = Legend(labels = c('E8.5','E10.5','E12.5','E14.5','E16.5','P3'), title = "Stage", legend_gp = gpar(fill = cividis(6)),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))
  pd=packLegend(lgd,lgd2)
  #plot scale
  
  png(paste0(OUTPUT_FILE_LOC,"/",i,"_npc96_var9k_knn12_t14_scale.png"),width=2500,height=2000)
  try(draw(ComplexHeatmap::Heatmap(clust,use_raster=FALSE,column_split = cells@meta.data$Dbhnom,cluster_rows = FALSE,cluster_cols_fix=FALSE,show_column_names = FALSE,show_column_dend = FALSE,row_names_gp = gpar(fontsize = 10,fontface='bold.italic'),row_names_side = 'left',top_annotation = ha,show_heatmap_legend = FALSE,col =colorRamp2(c(-10, 0, 10), c(plasma(20)[1], plasma(20)[10], plasma(20)[20]))),annotation_legend_list=pd))
  dev.off()
  
  
  ord=cells@assays$SCT@counts
  clust=as.matrix(ord)
  clust=clust[order(rowSums2(clust),decreasing=TRUE),]
  clust=clust[1:100,]
  col_fun = colorRamp2(c(0, 5, 10,20,30),c(plasma(30)[1], plasma(30)[5], plasma(30)[10],plasma(30)[20],plasma(30)[30]))
  lgd = Legend(col_fun = col_fun, title = "Expression",at = c(0,5,10,20,30), 
               labels = c("0",'5',"10",'20','>=30'),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))
  
  
  pd=packLegend(lgd,lgd2)
  
  #plot counts
  png(paste0(OUTPUT_FILE_LOC,"/",i,"_npc96_var9k_knn12_t14_counts.png"),width=2500,height=2000)
  draw(ComplexHeatmap::Heatmap(clust,cluster_rows = FALSE,cluster_cols_fix=TRUE,show_column_names = FALSE,show_column_dend = FALSE,row_names_gp = gpar(fontsize = 10,fontface='bold.italic'),row_names_side = 'left',top_annotation = ha,show_heatmap_legend = FALSE,col =colorRamp2(c(0, 5, 10,20,30), c(plasma(50)[1], plasma(50)[5], plasma(50)[10],plasma(50)[20],plasma(30)[30]))),annotation_legend_list=pd)
  dev.off()    
  
  
}




#######################NOW IDENTIFY ALL THE DBH
Dbhnom=rep('',ncol(cmlin@assays$RNA@counts))
Dbhnom[which(cmlin@assays$RNA@counts['Dbh',]>0)]='Dbh_positive'
Dbhnom[which(cmlin@assays$RNA@counts['Dbh',]==0)]='Dbh_negative'
cmlin@meta.data$Dbhnom=as.character(Dbhnom)
cmlin@meta.data$Dbhnom=as.factor(cmlin@meta.data$Dbhnom)

Dbhp=subset(cmlin,subset=Dbhnom=='Dbh_positive')
Dbhp=SCTransform(Dbhp, assay='RNA',new.assay.name='SCT',ncells=2500,variable.features.n=9000,return.only.var.genes=FALSE,seed.use=804L, vars.to.regress = c('nCount_RNA'), verbose = TRUE, min_cells=3,do.scale=TRUE)

phatedbh=phate(t(Dbhp@assays$SCT@scale.data),ndim=60L,mds.solver='smacof',npca=20L,verbose=TRUE,n.jobs=1,seed=804L)

#The below represents the final processed R object of the Dbh +ve CMs.
saveRDS(Dbhp,paste0(OUTPUT_FILE_LOC,'/Dbh_CMs_data.csv'))
saveRDS(phatedbh,paste0(OUTPUT_FILE_LOC,'/Dbh_CMs_PHATE.csv'))

################# Cardiac development dataset construction for integration by BGI

#####now need to bring in the other cell types and go from there to produce the final seurat obj for Tianyi
#keep AM,eVM-trab,VM,eVM,AVN,PKJ,VM-trab,eAM,SAN,AM-CCS, add in fibroblast,EC,ED,EP,Immune,VSMC,Blood,Plat

pbmc@meta.data$identssplit=pbmc@meta.data$updated8

tianyi=merge(subset(cmlin,subset=newtypez=='AM'),y=c(subset(cmlin,subset=newtypez=='eVM-trab'),subset(cmlin,subset=newtypez=='VM'),
                                                     subset(cmlin,subset=newtypez=='eVM'),subset(cmlin,subset=newtypez=='AVN'),subset(cmlin,subset=newtypez=='PKJ'),
                                                     subset(cmlin,subset=newtypez=='VM-trab'),subset(cmlin,subset=newtypez=='eAM'),subset(cmlin,subset=newtypez=='SAN'),
                                                     subset(cmlin,subset=newtypez=='AM-CCS'),subset(pbmc,subset=identssplit=='Endocardium'),subset(pbmc,subset=identssplit=='Epicardium'),
                                                     subset(pbmc,subset=identssplit=='Endothelium'),subset(pbmc,subset=identssplit=='Fibroblast-like'),subset(pbmc,subset=identssplit=='Haematopoietic Progenitors'),
                                                     subset(pbmc,subset=identssplit=='Immune Cells'),subset(pbmc,subset=identssplit=='Neural Crest'),subset(pbmc,subset=identssplit=='Platelets'),
                                                     subset(pbmc,subset=identssplit=='Smooth Muscle-like')),add.cell.ids=c('AM__','eVM-trab__','VM__','eVM__','AVN__','PKJ__','VM-trab__','eAM__',
                                                                                                                           'SAN__','AM-CCS__','Endocardium__','Epicardium__','Endothelium__','Fibroblast-like__','Erythroid Cells__',
                                                                                                                           'Immune Cells__','Neural Crest__','Platelets__','Smooth Muscle-like__'))


###########now need to extract the names - added __ to cell ids to ensure can track them

tianyi@meta.data$newones=sub("__.*", "", colnames(tianyi@assays$RNA))
tianyi@meta.data$newones=as.factor(tianyi@meta.data$newones)

####now do SCT then save
tianyi=SCTransform(tianyi, assay='RNA',new.assay.name='SCT',ncells=12000,variable.features.n=9000,return.only.var.genes=FALSE,seed.use=804L, vars.to.regress = c('nCount_RNA'), verbose = TRUE, min_cells=3,do.scale=TRUE)

saveRDS(tianyi,paste0(INT_FILE_LOC,"/Final_cell_types.csv"))


#export to SMASHpy


#So want to export to SMASHpy to get list of genes
tianyi=readRDS("C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/Final_cell_types.csv")
#Want to check that it's giving the same clusters as smashpy
clusty=read.csv2("C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/plotlys/cmlineage/forST/clusters.csv")


bigo=t(tianyi@assays$RNA@counts)
rownames(bigo)=NULL
bigo1=bigo[,1:10000]
bigo1=as.matrix(bigo1)
write_csv2(as.data.frame(bigo1),"C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/plotlys/cmlineage/forST/forSMASH_1_10k.csv")
bigo2=bigo[,10001:20000]
bigo2=as.matrix(bigo2)
write_csv2(as.data.frame(bigo2),"C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/plotlys/cmlineage/forST/forSMASH_10001_20k.csv")
bigo3=bigo[,20001:27933]
bigo3=as.matrix(bigo3)
write_csv2(as.data.frame(bigo3),"C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/plotlys/cmlineage/forST/forSMASH_20001_28k.csv")
#Also export the cell types

write.csv2(as.character(tianyi@meta.data$newones),paste0(INT_FILE_LOC,"/clusters.csv"))

#export to RCTD

###############refine them now, so only want like the top 400, probably stopping around 600 max

tianyi=readRDS("C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/Final_cell_types.csv")
library(future)
plan("multiprocess", workers = 25)
options(future.globals.maxSize = 10000 * 1024^3) 
vargenesout=list()

cells=tianyi
subbed=seq(0,18)
Idents(cells)=cells@meta.data$cell_identities
names(subbed)=levels(cells@active.ident)
pbmz=RenameIdents(cells,subbed)
pbmz@meta.data$multimark=pbmz@active.ident
pbmz@
  pbmz=SCTransform(pbmz, assay='RNA',new.assay.name='SCT',ncells=ncol(pbmz@assays$RNA@counts),variable.features.n=4000,return.only.var.genes=FALSE,seed.use=804L,verbose = TRUE, min_cells=3,do.scale=FALSE,do.correct=FALSE)
pbmy=CreateSeuratObject(pbmz@assays$SCT@data[VariableFeatures(pbmz),])
vargenesout=VariableFeatures(pbmz)
pbmy@meta.data$multimark=pbmz@active.ident
pbmc.multimark=FindAllMarkers.multicore(pbmz,min_pct=0,logfc_threshold=0,only_pos=FALSE,resolution='multimark',nCores=25)
saveRDS(pbmc.multimark,paste0("C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/plotlys/cmlineage/genesout2d/Dbh/Tianyi_9kvarsct_markers.csv"))
saveRDS(vargenesout,paste0("C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/plotlys/cmlineage/genesout2d/Dbh/Tianyi_9kvargenesout.csv"))
names(pbmc.multimark)=levels(cells@active.ident)
pbmc.multimark$SAN[sort(c(which(rownames(pbmc.multimark$SAN)%in%genei))),]

#get all the DE genes
de=list()
for(i in names(pbmc.multimark)){
  de[[i]]=rownames(pbmc.multimark[[i]][which(pbmc.multimark[[i]]$p_val_adj<=0.05),])
}

#get all the top 100 genes
top100=list()
for(i in levels(pbmz)){
  cells=subset(pbmz,subset=multimark==i)
  ord=cells@assays$SCT@scale.data
  clust=as.matrix(ord)
  clust=clust[order(rowSums2(clust),decreasing=TRUE),]
  clust=clust[1:100,]
  top100[[i]]=rownames(clust)
}

#get the final gene lists
genesforST=sort(unique(intersect(unlist(top100),unlist(de))))

forST=tianyi@assays$RNA@counts[which(rownames(tianyi@assays$RNA@counts)%in%genesforST),]
forST=CreateSeuratObject(forST)
forST@meta.data$cell_identities=pbmz@meta.data$newones

forST=SCTransform(forST, assay='RNA',new.assay.name='SCT',ncells=ncol(forST@assays$RNA@counts)/5,variable.features.n=1000,return.only.var.genes=FALSE,seed.use=804L,verbose = TRUE, min_cells=3,do.scale=TRUE,do.correct=TRUE)

for(i in levels(forST@meta.data$cell_identities)){
  cells=subset(forST,subset=cell_identities==i)
  stogy=substr(colnames(cells[["RNA"]]@counts),0L,4L)
  stogy[which(stogy=='E10.')]='E10.5'
  stogy[which(stogy=='E12.')]='E12.5'
  stogy[which(stogy=='E14.')]='E14.5'
  stogy[which(stogy=='E16.')]='E16.5'
  stogy[which(stogy=='Neo_')]='P2'
  
  #create the heatmap annotation for stages
  
  ha=HeatmapAnnotation(Stage=stogy,col=list(Stage=c('E8.5'=cividis(6)[1],'E10.5'=cividis(6)[2], 'E12.5'=cividis(6)[3],'E14.5'=cividis(6)[4],'E16.5'=cividis(6)[5],'P2'=cividis(6)[6])),show_legend = FALSE,show_annotation_name = FALSE)
  
  ord=cells@assays$SCT@scale.data
  clust=as.matrix(ord)
  clust=clust[order(rowSums2(clust),decreasing=TRUE),]
  #clust=clust[intersect(earlygenes,rownames(clust)),]
  clust=clust[1:100,]
  print(i)
  print(rownames(clust))
  col_fun = colorRamp2(c(-10, 0, 10),c(plasma(20)[1], plasma(20)[10], plasma(20)[20]))
  lgd = Legend(col_fun = col_fun, title = "Expression",at = c(-10,0, 10), 
               labels = c("Low",'Zero', "High"),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))
  
  lgd2 = Legend(labels = c('E8.5','E10.5','E12.5','E14.5','E16.5','P2'), title = "Stage", legend_gp = gpar(fill = cividis(6)),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))
  pd=packLegend(lgd,lgd2)
  #plot scale
  
  png(paste0("C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/plotlys/cmlineage/forST/heatmaps/forST_",i,"_scaled_down_final.png"),width=1500,height=1000)
  try(draw(ComplexHeatmap::Heatmap(clust,use_raster=FALSE,cluster_rows = FALSE,cluster_cols_fix=FALSE,show_column_names = FALSE,show_column_dend = FALSE,row_names_gp = gpar(fontsize = 10,fontface='bold.italic'),row_names_side = 'left',top_annotation = ha,show_heatmap_legend = FALSE,col =colorRamp2(c(-10, 0, 10), c(plasma(20)[1], plasma(20)[10], plasma(20)[20]))),annotation_legend_list=pd))
  dev.off()
  
}
saveRDS(forST,"C:/Users/agrassamrowe/Documents/CMs/3_5_21_new/filt/seurat/alltypes/plotlys/cmlineage/forST/forST.csv")




######################## Thus ends the first script

