###This script is for reproduction of final figures, supplementary figures, and supplementary tables and data.
#We present each figure independently in 'chronological order'
#Loading in the data necessarily for each figure separately so you don't have to run each one.
#However, we attempt to use common names across figures so if you want to run multiple using the same objects:
#please only load the data once to save you having to wait for large files to read from disc.

#### Read in libraries and set-up global variables - please do even if you only want to run a single figure later on.

#### Main Figures

#### Supplementary Figures

#### Supplementary Tables

#### Supplementary Data


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


### Main Figures

#Figure 1b
#load in cmlin
#load in phatecm
colalt=c('royalblue3','tomato4','pink3','maroon2','lightsteelblue4','springgreen4','yellowgreen','goldenrod3','salmon4','lightskyblue3','darkkhaki','goldenrod1','chartreuse1','aquamarine','honeydew3')
niceones=c('steelblue1','darkgoldenrod2','seagreen3','dodgerblue3','darkorchid1','royalblue3','cornflowerblue','brown2','blueviolet','deeppink','lightsalmon','orchid4','olivedrab4','burlywood2','palegreen4','yellow','mediumslateblue','honeydew4','gold2','firebrick2','orange4','yellow4','cyan1','firebrick4','cadetblue','blue3')

chango=as.character(cmlin@meta.data$newtypez)
chango[which(chango=='eVM-trab')]='Early Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='AHF')]='Heart Fields'
chango[which(chango=='FHF')]='Heart Fields'
chango[which(chango=='pSHF')]='Heart Fields'
chango[which(chango=='VM')]='Ventricular Cardiomyocytes'
chango[which(chango=='AM')]='Atrial Cardiomyocytes'
chango[which(chango=='eVM')]='Immature Ventricular Cardiomyocytes'
chango[which(chango=='AVN')]='Atrioventricular Node'
chango[which(chango=='PKJ')]='Purkinje Fibres'
chango[which(chango=='pHT')]='Primary Heart Tube'
chango[which(chango=='VM-trab')]='Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='eAM')]='Immature Atrial Cardiomyocytes'
chango[which(chango=='EDCM')]='Endocardial gene-rich Cardiomyocytes'
chango[which(chango=='DevCM')]='Developing Cardiomyocytes'
chango[which(chango=='SAN')]='Sinoatrial Node'
chango[which(chango=='AM-CCS')]='Atrial Conduction System'
chango[which(chango=='Apop')]='Misc.'
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                              'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                              'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                              'Sinoatrial Node','Atrial Conduction System','Misc.'))
clustcolours2=data.frame(sort(unique(chango)),colalt)
rownames(clustcolours2)=clustcolours2[,1]

ggplot(data=as.data.frame(phatecm$embedding))+
  geom_point(aes(PHATE1,PHATE3,color=factor(chango)),show.legend=FALSE)+
  scale_colour_manual(values=c(clustcolours2[which(unique(as.factor(chango))%in%clustcolours2[,1]),2]))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),panel.background = element_blank(),axis.title = element_blank())

#Figure 1c

chango=as.character(Dbhp@meta.data$newtypez)
chango[which(chango=='eVM-trab')]='Early Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='AHF')]='Heart Fields'
chango[which(chango=='FHF')]='Heart Fields'
chango[which(chango=='pSHF')]='Heart Fields'
chango[which(chango=='VM')]='Ventricular Cardiomyocytes'
chango[which(chango=='AM')]='Atrial Cardiomyocytes'
chango[which(chango=='eVM')]='Immature Ventricular Cardiomyocytes'
chango[which(chango=='AVN')]='Atrioventricular Node'
chango[which(chango=='PKJ')]='Purkinje Fibres'
chango[which(chango=='pHT')]='Primary Heart Tube'
chango[which(chango=='VM-trab')]='Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='eAM')]='Immature Atrial Cardiomyocytes'
chango[which(chango=='EDCM')]='Endocardial gene-rich Cardiomyocytes'
chango[which(chango=='DevCM')]='Developing Cardiomyocytes'
chango[which(chango=='SAN')]='Sinoatrial Node'
chango[which(chango=='AM-CCS')]='Atrial Conduction System'
chango[which(chango=='Apop')]='Misc.'
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                              'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                              'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                              'Sinoatrial Node','Atrial Conduction System','Misc.'))

#Delete the colour for Misc.
colalt=c('royalblue3','tomato4','pink3','maroon2','lightsteelblue4','springgreen4','yellowgreen','goldenrod3','salmon4','lightskyblue3','darkkhaki','goldenrod1','chartreuse1','aquamarine')
clustcolours2=data.frame(sort(unique(chango)),colalt)

ggplot(data=as.data.frame(phatedbh$embedding))+
  geom_point(aes(PHATE1,PHATE3,color=factor(chango)),size=4.5,show.legend=FALSE)+
  scale_colour_manual(values=c(clustcolours2[which(levels(as.factor(chango))%in%clustcolours2[,1]),2]))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),panel.background = element_blank(),axis.title = element_blank())


#Figure 1d
chango=as.character(cmlin@meta.data$newtypez)
chango[which(chango=='eVM-trab')]='Early Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='AHF')]='Heart Fields'
chango[which(chango=='FHF')]='Heart Fields'
chango[which(chango=='pSHF')]='Heart Fields'
chango[which(chango=='VM')]='Ventricular Cardiomyocytes'
chango[which(chango=='AM')]='Atrial Cardiomyocytes'
chango[which(chango=='eVM')]='Immature Ventricular Cardiomyocytes'
chango[which(chango=='AVN')]='Atrioventricular Node'
chango[which(chango=='PKJ')]='Purkinje Fibres'
chango[which(chango=='pHT')]='Primary Heart Tube'
chango[which(chango=='VM-trab')]='Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='eAM')]='Immature Atrial Cardiomyocytes'
chango[which(chango=='EDCM')]='Endocardial gene-rich Cardiomyocytes'
chango[which(chango=='DevCM')]='Developing Cardiomyocytes'
chango[which(chango=='SAN')]='Sinoatrial Node'
chango[which(chango=='AM-CCS')]='Atrial Conduction System'
chango[which(chango=='Apop')]='Misc.'
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                              'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                              'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                              'Sinoatrial Node','Atrial Conduction System','Misc.'))


#Exclude Misc from heatmap
goiheatclusts=factor(as.character(chango[which(chango!='Misc.')]),levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                                                                           'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                                                                           'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                                                                           'Sinoatrial Node','Atrial Conduction System'))

neworder=c('Sln','Myh7','Hopx','Hyal2','Myl2','Myl7','Vsnl1','Shox2','Cacna2d2','Kcne1','Sema3a','Bmp4','Tmsb4x','Isl1','Tbx1')
goi=unlist(neworder)
goiheat=cmlin@assays$SCT@scale.data[which(rownames(cmlin@assays$SCT@scale.data)%in%unlist(goi)),]

clusts=as.matrix(Matrix(nrow=length(c(intersect(rownames(cmlin@assays$SCT@scale.data),unlist(goi)))),ncol=length(levels(goiheatclusts))))
colnames(clusts)=levels(goiheatclusts)
rownames(clusts)=c(intersect(rownames(cmlin@assays$SCT@scale.data),unlist(goi)))
stgmx=as.matrix(Matrix(nrow=14,ncol=6))
rownames(stgmx)=levels(goiheatclusts)
idents=c('E8.','E10','E12','E14','E16','Neo')
corrected=c('E8.5','E10.5','E12.5','E14.5','E16.5','P3')
colnames(stgmx)=corrected


out=0
for(j in c(levels(goiheatclusts))){
  #take average for the clusters, mean - popoulate the matrix
  #start with j for each cell type
  out = out +1
  type=goiheat[,which(goiheatclusts==j)]
  #print(rownames(type)[1:5])
  
  for(m in 1:nrow(type)){
    clusts[m,match(j,colnames(clusts))]=mean(type[m,])
  }
  #print(rownames(clusts)[1:5])
  stages=c(substr(colnames(type),0L,3L))
  #make a loop that cleans it for you
  for(jk in 1:6){
    stages[which(stages==idents[jk])]=corrected[jk]  
  }
  for(nm in 1:6){
    stgmx[j,nm]=length(which(stages==corrected[nm]))/length(stages)
  }
}

#################so we've normalised the rows now across the row, but now we need to normalise between the rows
#such that expression level isn't relevant, need to map it across a sort of, between 0 and 1 in each row
for(k in 1:nrow(clusts)){
  #divide it all by the average of the avg for each gene
  for(n in 1:length(clusts[k,])){
    x=clusts[k,n]
    #scale it between 0 and 1
    x=((x-min(clusts[k,]))/(max(clusts[k,])-min(clusts[k,])))
    clusts[k,n]=x
  }
}
#nice
clusts=clusts[match(intersect(intersect(unlist(goi),rownames(cmlin@assays$SCT@scale.data)),rownames(clusts)),rownames(clusts)),]
#here put the names of the order you want them in
typepick=c('Atrial Cardiomyocytes','Ventricular Cardiomyocytes','Trabecular Ventricular Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Immature Ventricular Cardiomyocytes','Immature Atrial Cardiomyocytes','Sinoatrial Node','Atrial Conduction System','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes','Heart Fields')

clustz=clusts[,intersect(typepick
                         ,colnames(clusts))]

col_fun = colorRamp2(c(0, 0.5, 1), c(plasma(20)[1], plasma(20)[10], plasma(20)[20]))
lgd = Legend(col_fun = col_fun, title = "Expression",at = c(0, 1), 
             labels = c("Low", "High"),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))


pd=packLegend(lgd)
draw(ComplexHeatmap::Heatmap(t(as.matrix(clustz)), row_names_max_width = max_text_width(
  rownames(t(as.matrix(clustz))), 
  gp = gpar(fontsize = 11)
),cluster_rows = FALSE,cluster_columns=FALSE,row_names_gp = gpar(fontsize = 11,fontface='bold'),column_names_side = 'top',column_names_rot = 70,row_names_side = 'left',show_heatmap_legend = FALSE,col =c(plasma(20)),column_names_centered = TRUE,column_names_gp = gpar(fontsize=12,fontface='bold.italic')),annotation_legend_list=pd)


####Figure 2

#Figure 2b
#Pre-processing and initial unsupervised clustering was performed by BGI collaborators using Seurat, as detailed in the methods, and we then assigned labels.
#This, as above, involved visualisation using a range of Seurat functions such as VlnPlot, SpatialDimPlot, and SpatialFeatureplot and downstream reclustering.
#Here we load in the R objects already classified with progress detailed in the meta.data from newclusts through to Identities as the final classification

E12y=readRDS(paste0(INT_FILE_LOC,"/E12_5_Unsupervised.csv"))
E14y=readRDS(paste0(INT_FILE_LOC,"/E14_5_Unsupervised.csv"))
Neoy=readRDS(paste0(INT_FILE_LOC,"/Neo_Unsupervised.csv"))

colz=c()
colz['Trab-VM']='black'
colz['VM']='grey'
colz['AM']='springgreen4'
colz['AVJ']='moccasin'
colz['Epicardium']='orchid4'
colz['SAN']='darkorange3'
colz['Immune Cells']='yellow'
colz['Low quality']='honeydew4'
colz['AVN']='deeppink'
colz['Vessels - internal']='goldenrod1'
colz['RBC']='firebrick4'
colz['VM']='grey'
colz['Vessels']='springgreen2'
colz['His-BB']='darkslateblue'
colz['Compact VMs']='yellowgreen'
colz['Septum VM']='chartreuse1'
colz=colz[match(levels(E14y@meta.data$Identities),names(colz))]

p=0
stage=c('E12','E14','P3')
for(x in c(E12y,E14y,Neoy)){
  p=p+1
  stage_out=stage[p]
  png(paste0(OUTPUT_FILE_LOC,"/Unsupervised_ST_clustering_",i,".png"),width=1100,height=800)
  
  if(p==1){
    i='slice1_A1_R'
  }
  if(p==2){
    i='slice1_LU_p'  
  }
  if(p==3){
    i='slice1_4LB_p'  
  }
  
  Idents(x)='Identities'
  #Depending on your settings, you might want to play around with pt.size.factor +/- stroke to get the best image.
  #We did not use pt.size.factor in a manner correlated with true 'dot' size
  
  SpatialDimPlot(x,images=c(i),cols=colz,pt.size.factor=1.7,stroke=0.001)+NoLegend()  
  dev.off()
}


#Figure 2ci

merged=merge(E12y,y=c(E14y,Neoy),add.cell.ids = c('E12.5','E14.5','P3'))

merged@meta.data$Identities=c(as.character(E12y@meta.data$Identities),as.character(E14y@meta.data$Identities),as.character(Neoy@meta.data$Identities))
merged@meta.data$Identities=factor(merged@meta.data$Identities,levels=c('SAN','AM','AVN','AVJ','His-BB','Trab-VM','Septum VM','VM','Compact VMs'))

ploting=subset(merged,subset=Identities=='SAN'|Identities=='AM'|Identities=='AVN'|Identities=='AVJ'|Identities=='His-BB'|Identities=='Trab-VM'|Identities=='Compact VMs'|Identities=='Septum VM'|Identities=='VM')
ploting@meta.data$Identities=factor(ploting@meta.data$Identities,levels=rev(c('SAN','AM','AVN','AVJ','His-BB','Trab-VM','Septum VM','VM','Compact VMs')))

print(DotPlot(ploting,features=c('Dbh','WPRE','Hcn4','Cacna2d2','Sln','Myl7','Cntn2','Epha4','Nptn','Ntm','Slit2','Slc22a1',
                                 'Rspo3','Tbx3','Nkx2-5','Kcnj5','Nppa','Bmp10','Myl2','Hey2'),col.min=0,col.max=1,dot.min=0,scale.min=0,group.by = 'Identities',cols=c('RdBu'),assay='SCT')+
        scale_color_gradient2(low=plasma(10)[1],mid=plasma(10)[5],midpoint=0.5,high=plasma(10)[10],aesthetics='color',breaks=c(0.001,1),n.breaks=2,labels=c('Low','High'))+
        scale_size_continuous(range=c(3,10.5),breaks=c(1,5,10,15,20,25,50,75,100),labels = c('1','5','10','15','20','25','50','75','100'))+
        theme(axis.text.x=element_text(angle=45,size=15,vjust=0.6,face='bold.italic')))


#Figure 2cii


i='slice1_LU_p'    
mergco=T
mergcoords=Neoy@images[[i]]@coordinates
coord1=unlist(mergcoords$row)
coord2=unlist(mergcoords$col)

for(gene in unique(c('Myl7','Gja5','Hcn4','Cacna2d2','Cntn2','Id2','WPRE','Dbh'))){
  if(gene %in% rownames(Neoy@assays$SCT@counts)){
    if(mergco){
      genefor1=log(Neoy@assays$SCT@counts[gene,which(colnames(Neoy@assays$SCT@counts)%in%rownames(mergcoords))]+1)
    }
    ploty=data.frame(Coord1=coord1,Coord2=coord2*-1,Gene1=genefor1)
    
    
    if(max(ploty$Gene1)!=min(ploty$Gene1[which(ploty$Gene1>0)])){
      png(paste0(OUTPUT_FILE_LOC,"/",gene,'_P3_',i,'.png'),height=746,width=1156)  
      print(ggplot(as.data.frame(ploty),aes(Coord1, Coord2)) +
              geom_point(colour=c('grey'),size=1.5,show.legend = FALSE,alpha=c(0.6)) +
              geom_point(aes(color=Gene1),data=~subset(.,Gene1>0),size=3,alpha=1,show.legend=TRUE) +
              scale_color_gradientn(colours=c(viridis(11)[1],viridis(11)[11]),name=gene,labels=c('Low','High'),breaks=c(min(ploty$Gene1[which(ploty$Gene1>0)]),max(ploty$Gene1)))+
              labs(title=paste0(gene,' Log Expression ','Slice ',i))+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),axis.title = element_blank(),
                    axis.text.x = element_blank(),axis.text.y = element_blank(),text=element_text(size=20),panel.background = element_blank(), legend.position = 'right',legend.key.size=unit(1,'cm'),legend.title=element_text(size=19,face='bold'),legend.text=element_text(size=15,face='bold'),axis.line = element_blank()))
      
      
      dev.off()
    }
  }
}


####Now we have completed the panels for the main Figures

#We now address the panels for the Supplementary Figures

#Figure S2 panels
#Left-hand panels
pbmc@meta.data$lognumi=log10(pbmc@meta.data$nUMI)
FeaturePlot(pbmc,features='lognumi',raster=FALSE)
FeaturePlot(pbmc,features='mitoRatio',raster=FALSE)
FeaturePlot(pbmc,features='nGene',raster=FALSE)
#Right-hand panels
VlnPlot(pbmc,features='nUMI',group.by='orig.ident',raster=FALSE)+scale_y_continuous(breaks=c(seq(0,550000,5000)))+NoLegend()
VlnPlot(pbmc,features='mitoRatio',group.by='orig.ident',raster=FALSE)+scale_y_continuous(breaks=c(seq(0,0.20,0.05)))+NoLegend()
VlnPlot(pbmc,features='nGene',group.by='orig.ident',raster=FALSE)+scale_y_continuous(breaks=c(seq(0,6000,2000)))+NoLegend()

#Figure S3

#Panel A i
niceones=c('steelblue1','darkgoldenrod2','seagreen3','dodgerblue3','darkorchid1','royalblue3','cornflowerblue','brown2','blueviolet','lightyellow3','lightsalmon','burlywood2','olivedrab4','orchid4','darkslateblue','yellow','mediumslateblue','honeydew4','gold2','firebrick2','orange4','navyblue','cyan1','lightpink3','cadetblue','blue3')
colordf=data.frame(colours=niceones,row.names=levels(pbmc@meta.data$updated8))
ggplot(data=as.data.frame(pbmc@reductions$umap@cell.embeddings))+
  geom_point(aes(UMAP_1,UMAP_2,color=factor(pbmc@meta.data$updated8)),show.legend=FALSE)+
  scale_colour_manual(values=c(colordf[which(rownames(colordf)%in%levels(pbmc@meta.data$updated8)),1]))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),panel.background = element_blank())
#Panel A ii

pbmc@meta.data$stage2=factor(as.character(pbmc@meta.data$orig.ident),levels=c('E8','E10','E12','E14','E16','Neo'))
DimPlot(pbmc,group.by='stage2',raster=FALSE,cols=cividis(6))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y=element_blank(),
                                                                   axis.title.x = element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x = element_blank(),text=element_text(size=20),panel.background = element_blank(), axis.line = element_line(colour='black',size=3))

#Panel A iii
s.genes=cc.genes.updated.2019$s.genes
g2m.genes=cc.genes.updated.2019$g2m.genes
pbmc=CellCycleScoring(pbmc,s.features=s.genes,g2m.features = g2m.genes,set.ident=FALSE)
pbmc@meta.data$Phase=factor(pbmc@meta.data$Phase,levels=c('G1','G2M','S'))
DimPlot(pbmc,group.by='Phase',raster=FALSE,cols=viridis(3))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y=element_blank(),
                                                                  axis.title.x = element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.text.x = element_blank(),text=element_text(size=20),panel.background = element_blank(), axis.line = element_line(colour='black',size=3))


###Figure S5
Genes=c('Sln','Myh7','Myh6','Shox2','Hcn4','Myl7','Atp2a2','Tnni1','Tnni3','Fbn1','Dcn','Col1a1','Myh11','Rgs4','Col1a1','Upk3b','Pecam1',
        'Kdr','Fli1','Isl1','Hbb-y','Myog','Hand1','Pdgfra','Sox2','Tbx6','Otx2','Crabp2','Foxd3','Phox2b','Olig3','Pax6','Tubb3','Irx1','Ttr',
        'Snrpg','C1qb','Pf4')
for(gene in Genes){
  png(paste0(OUTPUT_FILE_LOC,"/",gene,"_UMAP_PLOT.png"))
  FeaturePlot(pbmc,features=gene,raster=FALSE,slot='counts')
  dev.off()
}

#Figure S6

#Panel a
chango=as.character(pbmc@meta.data$updated8)
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Cardiac Conduction System','Ventricular Cardiomyocytes_Myh6','Ventricular Cardiomyocytes','Early Cardiomyocytes',
                              'Early VMs','Early Atrial Cardiomyocytes','ECM Ventricular Cardiomyocytes','ECM Atrial Cardiomyocytes','Fibroblast-like',
                              'Smooth Muscle-like','Epicardium','Endothelium','Endocardium','Cardiac Progenitors','Skeletal Muscle Progenitors','Mesoderm','Haematopoietic Progenitors',
                              'Platelets','Immune Cells','Mixed Replicating','Neuromesodermal Progenitors','Neural','Neural Crest','Brain/Spinal Cord','Endoderm'))

goiheatclusts=chango
neworder=c('Nppa','Myl7','Sln','Hcn4','Cacna2d2','Tbx3','Shox2','Cacna1g','Ttn','Myh6','Myl2','Myl7','Tpm1','Tnnt2','Myl4','Eef1a1','Tpt1','Actc1','Gja1','Ryr3','Kcnj3','Vsnl1','Fbn1','Fbln2','Postn','Dcn','Col1a1','Upk3b','Upk1b','Aldh1a2','Ecscr','Emcn','Cdh5','Kdr','Tbx1','Etv2','Bmp4','Hand1','Wnt2','Gata6','Myf6','Myog','Foxd1','Prrx2','Redrum','Hemgn','Hbb-bh1','Tal1','Pf4','Lyz2','Fcer1g','C1qb','Snrpg','Rps28','Ptma','T','Sox2','Otx2','Pax6','Sox10','Foxd3','Isl1','Neurog1','Tubb3','Oncecut2','Pyy','Epcam','Afp','Ttr')

goi=unlist(neworder)
goiheat=pbmc@assays$SCT@scale.data[which(rownames(pbmc@assays$SCT@scale.data)%in%unlist(goi)),]

clusts=as.matrix(Matrix(nrow=length(c(intersect(rownames(pbmc@assays$SCT@scale.data),unlist(goi)))),ncol=length(levels(goiheatclusts))))
colnames(clusts)=levels(goiheatclusts)
rownames(clusts)=c(intersect(rownames(pbmc@assays$SCT@scale.data),unlist(goi)))
stgmx=as.matrix(Matrix(nrow=26,ncol=6))
rownames(stgmx)=levels(goiheatclusts)
idents=c('E8.','E10','E12','E14','E16','Neo')
corrected=c('E8.5','E10.5','E12.5','E14.5','E16.5','P3')
colnames(stgmx)=corrected


out=0
for(j in c(levels(goiheatclusts))){
  #take average for the clusters, mean - popoulate the matrix
  #start with j for each cell type
  out = out +1
  type=goiheat[,which(goiheatclusts==j)]
  #print(rownames(type)[1:5])
  
  for(m in 1:nrow(type)){
    clusts[m,match(j,colnames(clusts))]=mean(type[m,])
  }
  #print(rownames(clusts)[1:5])
  stages=c(substr(colnames(type),0L,3L))
  #make a loop that cleans it for you
  for(jk in 1:6){
    stages[which(stages==idents[jk])]=corrected[jk]  
  }
  for(nm in 1:6){
    stgmx[j,nm]=length(which(stages==corrected[nm]))/length(stages)
  }
}

#################so we've normalised the rows now across the row, but now we need to normalise between the rows
#such that expression level isn't relevant, need to map it across a sort of, between 0 and 1 in each row
for(k in 1:nrow(clusts)){
  #divide it all by the average of the avg for each gene
  for(n in 1:length(clusts[k,])){
    x=clusts[k,n]
    #scale it between 0 and 1
    x=((x-min(clusts[k,]))/(max(clusts[k,])-min(clusts[k,])))
    clusts[k,n]=x
  }
}
#nice
clusts=clusts[match(intersect(intersect(unlist(goi),rownames(pbmc@assays$SCT@scale.data)),rownames(clusts)),rownames(clusts)),]
#here put the names of the order you want them in
typepick=c('Atrial Cardiomyocytes','Cardiac Conduction System','Ventricular Cardiomyocytes_Myh6','Ventricular Cardiomyocytes','Early Cardiomyocytes',
           'Early VMs','Early Atrial Cardiomyocytes','ECM Ventricular Cardiomyocytes','ECM Atrial Cardiomyocytes','Fibroblast-like',
           'Smooth Muscle-like','Epicardium','Endothelium','Endocardium','Cardiac Progenitors','Skeletal Muscle Progenitors','Mesoderm','Haematopoietic Progenitors',
           'Platelets','Immune Cells','Mixed Replicating','Neuromesodermal Progenitors','Neural','Neural Crest','Brain/Spinal Cord','Endoderm')

clustz=clusts[,intersect(typepick
                         ,colnames(clusts))]
#now need to make a thing where it takes the orig and makes into histogram
ha=rowAnnotation(Stages = anno_barplot(stgmx[intersect(typepick
                                                       ,rownames(stgmx)),], gp=gpar(fill=cividis(6)), 
                                       bar_width = 1, width = unit(2, "cm")),show_annotation_name=FALSE,show_legend=FALSE)
col_fun = colorRamp2(c(0, 0.5, 1), c(plasma(20)[1], plasma(20)[10], plasma(20)[20]))
lgd = Legend(col_fun = col_fun, title = "Expression",at = c(0, 1), 
             labels = c("Low", "High"),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))

lgd2 = Legend(labels = c('E8.5','E10.5','E12.5','E14.5','E16.5','P3'), title = "Stage", legend_gp = gpar(fill = cividis(6)),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))

pd=packLegend(lgd,lgd2)
draw(ComplexHeatmap::Heatmap(t(as.matrix(clustz)), row_names_max_width = max_text_width(
  rownames(t(as.matrix(clustz))), 
  gp = gpar(fontsize = 11)
),right_annotation = ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_gp = gpar(fontsize = 11,fontface='bold'),column_names_side = 'top',column_names_rot = 70,row_names_side = 'left',show_heatmap_legend = FALSE,col =c(plasma(20)),column_names_centered = TRUE,column_names_gp = gpar(fontsize=12,fontface='bold.italic')),annotation_legend_list=pd)




#Panel b
chango=as.character(cmlin@meta.data$newtypez)
chango[which(chango=='eVM-trab')]='Early Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='AHF')]='Heart Fields'
chango[which(chango=='FHF')]='Heart Fields'
chango[which(chango=='pSHF')]='Heart Fields'
chango[which(chango=='VM')]='Ventricular Cardiomyocytes'
chango[which(chango=='AM')]='Atrial Cardiomyocytes'
chango[which(chango=='eVM')]='Immature Ventricular Cardiomyocytes'
chango[which(chango=='AVN')]='Atrioventricular Node'
chango[which(chango=='PKJ')]='Purkinje Fibres'
chango[which(chango=='pHT')]='Primary Heart Tube'
chango[which(chango=='VM-trab')]='Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='eAM')]='Immature Atrial Cardiomyocytes'
chango[which(chango=='EDCM')]='Endocardial gene-rich Cardiomyocytes'
chango[which(chango=='DevCM')]='Developing Cardiomyocytes'
chango[which(chango=='SAN')]='Sinoatrial Node'
chango[which(chango=='AM-CCS')]='Atrial Conduction System'
chango[which(chango=='Apop')]='Misc.'
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                              'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                              'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                              'Sinoatrial Node','Atrial Conduction System','Misc.'))


#Exclude Misc from heatmap
goiheatclusts=factor(as.character(chango[which(chango!='Misc.')]),levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                                                                           'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                                                                           'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                                                                           'Sinoatrial Node','Atrial Conduction System'))

neworder=c('Myl7','Sln','Pln','Myh7','Myl2','Myl3','Myh6','Dcn','Ckm','Casq2','Hopx','Notch1','Notch3','Nppb','Nppa','Mmp2','Adamts1','Vcan','Hyal2','Actc1','Acta2','Tagln','Eef1a1','Shox2','Hcn4','Tbx3','Kcne1','Cacna2d2','Cacna1g','Bmp10','Ephb3','Sema3a','Tmsb10','Tmsb4x','Flt1','Emcn','Isl1','Pmp22','Krt18','Krt8','Hand1','Sfrp1','Gata5','Gata6','Irx5','Tcf15','Meox1','Lef1','Fgf10')

goi=unlist(neworder)
goiheat=cmlin@assays$SCT@scale.data[which(rownames(cmlin@assays$SCT@scale.data)%in%unlist(goi)),]

clusts=as.matrix(Matrix(nrow=length(c(intersect(rownames(cmlin@assays$SCT@scale.data),unlist(goi)))),ncol=length(levels(goiheatclusts))))
colnames(clusts)=levels(goiheatclusts)
rownames(clusts)=c(intersect(rownames(cmlin@assays$SCT@scale.data),unlist(goi)))
stgmx=as.matrix(Matrix(nrow=14,ncol=6))
rownames(stgmx)=levels(goiheatclusts)
idents=c('E8.','E10','E12','E14','E16','Neo')
corrected=c('E8.5','E10.5','E12.5','E14.5','E16.5','P3')
colnames(stgmx)=corrected


out=0
for(j in c(levels(goiheatclusts))){
  #take average for the clusters, mean - popoulate the matrix
  #start with j for each cell type
  out = out +1
  type=goiheat[,which(goiheatclusts==j)]
  #print(rownames(type)[1:5])
  
  for(m in 1:nrow(type)){
    clusts[m,match(j,colnames(clusts))]=mean(type[m,])
  }
  #print(rownames(clusts)[1:5])
  stages=c(substr(colnames(type),0L,3L))
  #make a loop that cleans it for you
  for(jk in 1:6){
    stages[which(stages==idents[jk])]=corrected[jk]  
  }
  for(nm in 1:6){
    stgmx[j,nm]=length(which(stages==corrected[nm]))/length(stages)
  }
}

#################so we've normalised the rows now across the row, but now we need to normalise between the rows
#such that expression level isn't relevant, need to map it across a sort of, between 0 and 1 in each row
for(k in 1:nrow(clusts)){
  #divide it all by the average of the avg for each gene
  for(n in 1:length(clusts[k,])){
    x=clusts[k,n]
    #scale it between 0 and 1
    x=((x-min(clusts[k,]))/(max(clusts[k,])-min(clusts[k,])))
    clusts[k,n]=x
  }
}
#nice
clusts=clusts[match(intersect(intersect(unlist(goi),rownames(cmlin@assays$SCT@scale.data)),rownames(clusts)),rownames(clusts)),]
#here put the names of the order you want them in
typepick=c('Atrial Cardiomyocytes','Ventricular Cardiomyocytes','Trabecular Ventricular Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Immature Ventricular Cardiomyocytes','Immature Atrial Cardiomyocytes','Sinoatrial Node','Atrial Conduction System','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes','Heart Fields')

clustz=clusts[,intersect(typepick
                         ,colnames(clusts))]

ha=rowAnnotation(Stages = anno_barplot(stgmx[intersect(typepick
                                                       ,rownames(stgmx)),], gp=gpar(fill=cividis(6)), 
                                       bar_width = 1, width = unit(2, "cm")),show_annotation_name=FALSE,show_legend=FALSE)
col_fun = colorRamp2(c(0, 0.5, 1), c(plasma(20)[1], plasma(20)[10], plasma(20)[20]))
lgd = Legend(col_fun = col_fun, title = "Expression",at = c(0, 1), 
             labels = c("Low", "High"),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))

lgd2 = Legend(labels = c('E8.5','E10.5','E12.5','E14.5','E16.5','P3'), title = "Stage", legend_gp = gpar(fill = cividis(6)),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))

pd=packLegend(lgd,lgd2)
draw(ComplexHeatmap::Heatmap(t(as.matrix(clustz)), row_names_max_width = max_text_width(
  rownames(t(as.matrix(clustz))), 
  gp = gpar(fontsize = 11)
),right_annotation = ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_gp = gpar(fontsize = 11,fontface='bold'),column_names_side = 'top',column_names_rot = 70,row_names_side = 'left',show_heatmap_legend = FALSE,col =c(plasma(20)),column_names_centered = TRUE,column_names_gp = gpar(fontsize=12,fontface='bold.italic')),annotation_legend_list=pd)


#Panel c
embedz=data.frame(phatedbh$embedding)
for(i in c('Tnnt2','Cacna2d2','Cpne5','Hcn4','Shox2')){
  tiff(paste0(OUTPUT_FILE_LOC,"/Dbh_bold_trans_",i,".tiff"),width=700,height=600)
  try(print(
    ggplot(embedz) +
      geom_point(aes(PHATE1, PHATE3), colour=c('grey'),size=6,show.legend = FALSE) +
      geom_point(data=embedz[which(Dbhp@assays$SCT@counts[i,]>0),],show.legend=TRUE,aes(PHATE1, PHATE3, colour=Dbhp@assays$SCT@counts[i,which(Dbhp@assays$SCT@counts[i,]>0)]),size=6) +
      labs(colour=paste0(i," Expression"),title=paste0(i, ' Expression'))+
      scale_colour_gradientn(colours=c(magma(max(Dbhp@assays$SCT@counts[i,])))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),text=element_text(size=20),panel.background = element_blank(), legend.position = 'right',legend.key.size=unit(1.5,'cm'),legend.title=element_text(size=19,face='bold'),legend.text=element_text(size=15,face='bold'),axis.line = element_line(colour = "black",size=3))
    
  ))
  dev.off()
}


#Figure S8
S8genes=c('Sln','Myh7','Myl2','Myl7','Hopx','Bmp10','Nppb','Actc1','Shox2','Hcn4','Cacna2d2','Cacna1g','Bmp2','Tbx3','Kcne1','Sema3a','Ephb3','Hand1','Ttn','Gata5','Tmsb4x','Emcn','Pmp22','Krt18','Tcf15','Meox1','Fgf10','Lef1')

for(i in S8genes){
  png(paste0(OUTPUT_FILE_LOC,"/",i,"_FigureS8_out.png"))
  try(print(
    ggplot(as.data.frame(phatecm$embedding)) +
      geom_point(aes(PHATE1, PHATE3), colour=c('grey'),size=1,show.legend = FALSE) +
      geom_point(data=as.data.frame(phatecm$embedding[which(cmlin@assays$RNA@counts[i,]>0),]),show.legend=TRUE,aes(PHATE1, PHATE3, colour=cmlin@assays$SCT@data[i,which(cmlin@assays$RNA@counts[i,]>0)]),size=1) +
      labs(colour=paste0(i," Expression"),title=paste0(i, ' Expression'))+
      scale_colour_gradientn(colours=c(magma(max(cmlin@assays$SCT@data[i,])))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),text=element_text(size=20),panel.background = element_blank(), legend.position = 'right',legend.key.size=unit(1.5,'cm'),legend.title=element_text(size=19,face='bold'),legend.text=element_text(size=15,face='bold'),axis.line = element_line(colour = "black",size=3))
    
  ))
  dev.off()
}

#Figure S9
chango=as.character(cmlin@meta.data$newtypez)
chango[which(chango=='eVM-trab')]='Early Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='AHF')]='Heart Fields'
chango[which(chango=='FHF')]='Heart Fields'
chango[which(chango=='pSHF')]='Heart Fields'
chango[which(chango=='VM')]='Ventricular Cardiomyocytes'
chango[which(chango=='AM')]='Atrial Cardiomyocytes'
chango[which(chango=='eVM')]='Immature Ventricular Cardiomyocytes'
chango[which(chango=='AVN')]='Atrioventricular Node'
chango[which(chango=='PKJ')]='Purkinje Fibres'
chango[which(chango=='pHT')]='Primary Heart Tube'
chango[which(chango=='VM-trab')]='Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='eAM')]='Immature Atrial Cardiomyocytes'
chango[which(chango=='EDCM')]='Endocardial gene-rich Cardiomyocytes'
chango[which(chango=='DevCM')]='Developing Cardiomyocytes'
chango[which(chango=='SAN')]='Sinoatrial Node'
chango[which(chango=='AM-CCS')]='Atrial Conduction System'
chango[which(chango=='Apop')]='Misc.'
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                              'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                              'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                              'Sinoatrial Node','Atrial Conduction System','Misc'))

goiheatclusts=factor(as.character(chango[which(chango!='Misc.')]),levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                                                                           'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                                                                           'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                                                                           'Sinoatrial Node','Atrial Conduction System'))


goi=c('Smoc2','Igfbp5','Cpne5','Pcdh17','Shox2','Slitrk5','Rgs6','Hcn4','Tbx3','Tbx2','Bmp2','Bmp4','Slc22a1','Ntm','Rspo3','Gfra2','Hcn1','Nkx2-5','Dkk3','Epha4','Nptn','Hcn2','Isl1','Tbx18','Actn2','Scn5a','Scn10a','Gja5','Kcne1','Etv1','Irx5','Irx3','Sema3a','Slit2','Bmp10','Tbx5','Tbx20','Cacna2d2','Cacna1g','Myl7','Myl4','Myh7','Myl2','Myl3','Hey1')


goiheat=cmlin@assays$SCT@counts[which(rownames(cmlin@assays$SCT@counts)%in%unlist(goi)),]

clusts=as.matrix(Matrix(nrow=length(c(intersect(rownames(cmlin@assays$SCT@counts),unlist(goi)))),ncol=length(levels(goiheatclusts))))
colnames(clusts)=levels(goiheatclusts)
rownames(clusts)=c(intersect(rownames(cmlin@assays$SCT@counts),unlist(goi)))
stgmx=as.matrix(Matrix(nrow=14,ncol=6))
rownames(stgmx)=levels(goiheatclusts)
idents=c('E8.','E10','E12','E14','E16','Neo')
corrected=c('E8.5','E10.5','E12.5','E14.5','E16.5','P2')
colnames(stgmx)=corrected


out=0
for(j in c(levels(goiheatclusts))){
  #take average for the clusters, mean - popoulate the matrix
  #start with j for each cell type
  out = out +1
  type=goiheat[,which(goiheatclusts==j)]
  #print(rownames(type)[1:5])
  
  for(m in 1:nrow(type)){
    clusts[m,match(j,colnames(clusts))]=mean(type[m,])
  }
  #print(rownames(clusts)[1:5])
  stages=c(substr(colnames(type),0L,3L))
  #make a loop that cleans it for you
  for(jk in 1:6){
    stages[which(stages==idents[jk])]=corrected[jk]  
  }
  for(nm in 1:6){
    stgmx[j,nm]=length(which(stages==corrected[nm]))/length(stages)
  }
}


for(i in 1:nrow(clusts)){
  for(j in 1:ncol(clusts)){
    if(isZero(clusts[i,j])){
      print(paste0(rownames(clusts)[i],' and ',colnames(clusts)[j],' is 0'))
      
    }
  }}
#need to remove any that are 0,
clustscal=apply(clusts, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))
#clustscal is now continously scaled ebtween 0 and 1
#################so we've normalised the rows now across the row, but now we need to normalise between the rows
#such that expression level isn't relevant, need to map it across a sort of, between 0 and 1 in each row
for(k in 1:nrow(clusts)){
  #divide it all by the average of the avg for each gene
  for(n in 1:length(clusts[k,])){
    x=clusts[k,n]
    #scale it between 0 and 1
    x=((x-min(clusts[k,]))/(max(clusts[k,])-min(clusts[k,])))
    clusts[k,n]=x
  }
}
#nice
clustscal2=clustscal[,match(intersect(intersect(unlist(goi),rownames(cmlin@assays$SCT@scale.data)),colnames(clustscal)),colnames(clustscal))]
#here put the names of the order you want them in

typepick=c('Sinoatrial Node','Atrioventricular Node','Atrial Conduction System','Purkinje Fibres')
clustz=clustscal2[intersect(typepick
                            ,rownames(clustscal2)),]

ha=rowAnnotation(Stages = anno_barplot(stgmx[intersect(typepick
                                                       ,rownames(stgmx)),], gp=gpar(fill=cividis(6)), 
                                       bar_width = 1, width = unit(2, "cm")),show_annotation_name=FALSE,show_legend=FALSE)
col_fun = colorRamp2(c(0, 0.5, 1), c(plasma(20)[1], plasma(20)[10], plasma(20)[20]))
lgd = Legend(col_fun = col_fun, title = "Expression",at = c(0, 1), 
             labels = c("Low", "High"),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))

lgd2 = Legend(labels = c('E8.5','E10.5','E12.5','E14.5','E16.5','P3'), title = "Stage", legend_gp = gpar(fill = cividis(6)),labels_gp = gpar(fontsize=11,fontface='bold.italic'),title_gp = gpar(fontsize = 12, fontface = "bold"))

pd=packLegend(lgd,lgd2)
draw(ComplexHeatmap::Heatmap(as.matrix(clustz), row_names_max_width = max_text_width(
  rownames(as.matrix(clustz)), 
  gp = gpar(fontsize = 11)
),right_annotation = ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_gp = gpar(fontsize = 11,fontface='bold'),column_names_side = 'top',column_names_rot = 70,row_names_side = 'left',show_heatmap_legend = FALSE,col =c(plasma(20)),column_names_centered = TRUE,column_names_gp = gpar(fontsize=12,fontface='bold.italic')),annotation_legend_list=pd)

#Figure S10
#UMAP of all cells
png(paste0(OUTPUT_FILE_LOC,"/",gene,"_UMAP_PLOT_DBH_FigS10.png"))
FeaturePlot(pbmc,features='Dbh',raster=FALSE,slot='counts',pt.size=4,order=TRUE)
dev.off()

#PHATE of CMLIN
i='Dbh'
png(paste0(OUTPUT_FILE_LOC,"/",i,"_PHATE_DBH_FigS10_out.png"))
try(print(
  ggplot(as.data.frame(phatecm$embedding)) +
    geom_point(aes(PHATE1, PHATE3), colour=c('grey'),size=1,show.legend = FALSE) +
    geom_point(data=as.data.frame(phatecm$embedding[which(cmlin@assays$RNA@counts[i,]>0),]),show.legend=TRUE,aes(PHATE1, PHATE3, colour=cmlin@assays$SCT@data[i,which(cmlin@assays$RNA@counts[i,]>0)]),size=1) +
    labs(colour=paste0(i," Expression"),title=paste0(i, ' Expression'))+
    scale_colour_gradientn(colours=c(magma(max(cmlin@assays$SCT@data[i,])))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),axis.text.y = element_blank(),text=element_text(size=20),panel.background = element_blank(), legend.position = 'right',legend.key.size=unit(1.5,'cm'),legend.title=element_text(size=19,face='bold'),legend.text=element_text(size=15,face='bold'),axis.line = element_line(colour = "black",size=3))
  
))
dev.off()


#Figure S11
chango=as.character(cmlin@meta.data$newtypez)
chango[which(chango=='eVM-trab')]='Early Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='AHF')]='Heart Fields'
chango[which(chango=='FHF')]='Heart Fields'
chango[which(chango=='pSHF')]='Heart Fields'
chango[which(chango=='VM')]='Ventricular Cardiomyocytes'
chango[which(chango=='AM')]='Atrial Cardiomyocytes'
chango[which(chango=='eVM')]='Immature Ventricular Cardiomyocytes'
chango[which(chango=='AVN')]='Atrioventricular Node'
chango[which(chango=='PKJ')]='Purkinje Fibres'
chango[which(chango=='pHT')]='Primary Heart Tube'
chango[which(chango=='VM-trab')]='Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='eAM')]='Immature Atrial Cardiomyocytes'
chango[which(chango=='EDCM')]='Endocardial gene-rich Cardiomyocytes'
chango[which(chango=='DevCM')]='Developing Cardiomyocytes'
chango[which(chango=='SAN')]='Sinoatrial Node'
chango[which(chango=='AM-CCS')]='Atrial Conduction System'
chango[which(chango=='Apop')]='Misc.'
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                              'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                              'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                              'Sinoatrial Node','Atrial Conduction System','Misc.'))

cmlin@meta.data$chango=chango

perce=c()
for(i in unique(Dbhp@meta.data$newtypez)){
  perce[i]=(paste0(i," ",((length(which(Dbhp@meta.data$newtypez==i))/length(which(cmlin@meta.data$newtypez==i)))*100),"%"))
  print(paste0(i," ",length(which(Dbhp@meta.data$newtypez==i))))
  print(paste0(i," ",length(which(cmlin@meta.data$newtypez==i))))
  print(paste0(i," ",length(which(cmlin@meta.data$newtypez==i))-length(which(Dbhp@meta.data$newtypez==i))))
  print(paste0(i," ",(length(which(Dbhp@meta.data$newtypez==i))/length(which(cmlin@meta.data$newtypez==i))*100)))
}


#get the %s
Dbh=list()
nomcap=list()
Pos=list()
tot=list()
x=0
for(i in levels(cmlin@meta.data$chango)){
  x=x+1
  Dbh[[x]]=length(which(cmlin@assays$RNA@counts['Dbh',which(cmlin_new@meta.data$chango==i)]>0))
  nomcap[[x]]=i
  Pos[[x]]='+'
  tot[[x]]=length(cmlin@assays$RNA@counts['Dbh',which(cmlin@meta.data$chango==i)])
  x=x+1
  Dbh[[x]]=length(which(cmlin@assays$RNA@counts['Dbh',which(cmlin@meta.data$chango==i)]==0))
  nomcap[[x]]=i
  Pos[[x]]='-'
  tot[[x]]=length(cmlin@assays$RNA@counts['Dbh',which(cmlin@meta.data$chango==i)])
}
unlist(Dbh)
unlist(nomcap)
outy=data.frame('Number'=unlist(Dbh),'Cell'=unlist(nomcap),'Pos_Neg'=unlist(Pos),'Total'=unlist(tot))

#let's do it as percentages with n numbers given
outy$Percent=round((outy$Number/outy$Total)*100,1)

dodge <- position_dodge(width = 0.9)

outy$Cell=factor(outy$Cell,levels=c('Purkinje Fibres','Atrioventricular Node','Sinoatrial Node',
                                    'Trabecular Ventricular Cardiomyocytes','Ventricular Cardiomyocytes',
                                    'Atrial Conduction System','Early Trabecular Ventricular Cardiomyocytes',
                                    'Atrial Cardiomyocytes','Immature Atrial Cardiomyocytes','Immature Ventricular Cardiomyocytes',
                                    'Primary Heart Tube','Developing Cardiomyocytes','Heart Fields','Endocardial gene-rich Cardiomyocytes',
                                    'Misc.'))

ggplot(outy, aes(x = Cell, y = Percent, fill = factor(Pos_Neg))) +
  geom_bar(stat = "identity", position = position_dodge())+theme(axis.text=element_text(family='Helvetica',size=20),axis.text.x=element_text(angle=90))



#Figure S13
#All the spatial plots
#Run for E12.5,E14.5, P3

for(slice in c('slice1_A1_R','slice1_RB_p','slice1_4LB_p')){
  mergco=T
  if(slice=='slice1_A1_R'){
    usy=E12y
    stg='E12_5'
  }
  elif(slice=='slice1_RB_p'){
    usy=E14y
    stg='E14_5'
  }
  if(slice=='slice1_4LB_p'){
    usy=Neoy
    stg='P3'
  }
  mergcoords=usy@images[[slice]]@coordinates
  coord1=unlist(mergcoords$row)
  coord2=unlist(mergcoords$col)
  
  for(gene in unique(c('Myh6','Myh7','Nppa','Myl7','Myl2','Myl4','Id2','Bmp10','Tbx18','Hcn4','Shox2','Smoc2','Igfbp5','Cpne5','Cntn2','Cacna2d2','Slit2'))){
    if(gene %in% rownames(usy@assays$SCT@counts)){
      if(mergco){
        genefor1=log(usy@assays$SCT@counts[gene,which(colnames(usy@assays$SCT@counts)%in%rownames(mergcoords))]+1)
      }
      ploty=data.frame(Coord1=coord1,Coord2=coord2*-1,Gene1=genefor1)
      
      
      if(max(ploty$Gene1)!=min(ploty$Gene1[which(ploty$Gene1>0)])){
        png(paste0(OUTPUT_FILE_LOC,"/",gene,'_',stg,'_',i,'.png'),height=746,width=1156)  
        print(ggplot(as.data.frame(ploty),aes(Coord1, Coord2)) +
                geom_point(colour=c('grey'),size=1.5,show.legend = FALSE,alpha=c(0.6)) +
                geom_point(aes(color=Gene1),data=~subset(.,Gene1>0),size=3,alpha=1,show.legend=TRUE) +
                scale_color_gradientn(colours=c(viridis(11)[1],viridis(11)[11]),name=gene,labels=c('Low','High'),breaks=c(min(ploty$Gene1[which(ploty$Gene1>0)]),max(ploty$Gene1)))+
                labs(title=paste0(gene,' Log Expression ','Slice ',slice))+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),axis.title = element_blank(),
                      axis.text.x = element_blank(),axis.text.y = element_blank(),text=element_text(size=20),panel.background = element_blank(), legend.position = 'right',legend.key.size=unit(1,'cm'),legend.title=element_text(size=19,face='bold'),legend.text=element_text(size=15,face='bold'),axis.line = element_blank()))
        
        
        dev.off()
      }
    }
  }
}


#Now run for adult
slicelist=list()
for(sliceo in c('C5')){
  slicelist[[sliceo]]=readRDS(paste0(INT_FILE_LOC,"/Adult_",sliceo,"_bin_seurat.rds"))
  mergco=T
  mergcoords=slicelist[[sliceo]]@images[[names(slicelist[[sliceo]]@images)[1]]]@coordinates
  coord1=unlist(mergcoords$row)
  coord2=unlist(mergcoords$col)
  
  
  for(gene in unique(c('Myh6','Myh7','Nppa','Myl7','Myl2','Myl4','Id2','Bmp10','Tbx18','Hcn4','Shox2','Smoc2','Igfbp5','Cpne5','Cntn2','Cacna2d2','Slit2'))){
    if(gene %in% rownames(slicelist[[sliceo]]@assays$SCT@counts)){
      if(mergco){
        genefor1=log(slicelist[[sliceo]]@assays$SCT@counts[gene,]+1)
      }
      ploty=data.frame(Coord1=coord1,Coord2=coord2*-1,Gene1=genefor1)
      
      
      if(max(ploty$Gene1)!=min(ploty$Gene1[which(ploty$Gene1>0)])){
        png(paste0(OUTPUT_FILE_LOC,"/Adult_",sliceo,'_',gene,'.png'),height=746,width=1156)  
        print(ggplot(as.data.frame(ploty),aes(Coord1, Coord2)) +
                geom_point(colour=c('grey'),size=1.5,show.legend = FALSE,alpha=c(0.6)) +
                geom_point(aes(color=Gene1),data=~subset(.,Gene1>0),size=3,alpha=1,show.legend=TRUE) +
                scale_color_gradientn(colours=c(viridis(11)[1],viridis(11)[11]),name=gene,labels=c('Low','High'),breaks=c(min(ploty$Gene1[which(ploty$Gene1>0)]),max(ploty$Gene1)))+
                labs(title=paste0(gene,' Log Expression ','Slice ',sliceo))+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),axis.title = element_blank(),
                      axis.text.x = element_blank(),axis.text.y = element_blank(),text=element_text(size=20),panel.background = element_blank(), legend.position = 'right',legend.key.size=unit(1,'cm'),legend.title=element_text(size=19,face='bold'),legend.text=element_text(size=15,face='bold'),axis.line = element_blank()))
        
        
        dev.off()
      }
    }
  }
}

#Figure S14
#Dotplots for all things
merged=merge(E12y,y=c(E14y,Neoy),add.cell.ids = c('E12.5','E14.5','P3'))

merged@meta.data$Identities=c(as.character(E12y@meta.data$Identities),as.character(E14y@meta.data$Identities),as.character(Neoy@meta.data$Identities))
merged@meta.data$Identities=factor(merged@meta.data$Identities,levels=rev(c('SAN','AM','AVN','AVJ','His-BB','Trab-VM','Septum VM','VM','Compact VMs','Epicardium','Vessels','Vessels - internal','RBC','Immune Cells','Low quality')))


ploting=subset(merged,subset=Identities=='SAN'|Identities=='AM'|Identities=='AVN'|Identities=='AVJ'|Identities=='His-BB'|Identities=='Trab-VM'|Identities=='Septum VM'|Identities=='VM'|Identities=='Epicardium'|Identities=='Vessels'|Identities=='Vessels - internal'|Identities=='RBC'|Identities=='Immune Cells'|Identities=='Low quality')
ploting@meta.data$Identities=factor(ploting@meta.data$Identities,levels=rev(c('SAN','AM','AVN','AVJ','His-BB','Trab-VM','Septum VM','VM','Compact VMs','Epicardium','Vessels','Vessels - internal','RBC','Immune Cells','Low quality')))

print(DotPlot(ploting,features=unique(c('Smoc2','Pcdh17','Gfra2','Igfbp5','Hcn4','Cacna2d2','Sln','Myl7','Cntn2','Epha4','Nptn','Ntm','Slit2','Slc22a1','Rspo3','Slitrk5','Tbx3','Nkx2-5','Kcnj5','Nppa','Bmp10','Myl2','Hey2','Cpne5','Shox2','Bmp2','Bmp4','Tbx3','Tbx2','Bmp10','Rgs6','Rgs5','Hcn4','Isl1','Tbx18','Actn2','Scn5a','Gja5','Gja1','Dkk3','Nkx2-5','Hcn1','Hcn2','Kcne1','Etv1','Irx5','Irx3','Sema3a','Sema3c','Slit2','Tbx5','Hcn4','Rspo3','Gja5','Bmp10','Cacna2d2','Myl7','Myl4','Myh6','Myh7','Myl2','Myl3','Hey1','Upk3b','Cxcl12','Rgs5','Hemgn','Fcer1g','Mgp','Eln','Sfrp2','Tagln','Fbln5','Tyrobp','Lyz2','Cd52')),col.min=0,col.max=1,dot.min=0,scale.min=0,group.by = 'Identities',cols=c('RdBu'),assay='SCT')+
        scale_color_gradient2(low=plasma(10)[1],mid=plasma(10)[5],midpoint=0.5,high=plasma(10)[10],aesthetics='color',breaks=c(0.001,1),n.breaks=2,labels=c('Low','High'))+
        scale_size_continuous(range=c(3,10.5),breaks=c(1,5,10,15,20,25,50,75,100),labels = c('1','5','10','15','20','25','50','75','100'))+
        theme(axis.text.x=element_text(angle=45,size=15,vjust=0.6,face='bold.italic')))

#Figure S15
#All stages DBH and WPRE
#Run for E12.5,E14.5, P3

for(slice in c('slice1_A1_R','slice1_RB_p','slice1_4LB_p')){
  mergco=T
  if(slice=='slice1_A1_R'){
    usy=E12y
    stg='E12_5'
  }
  elif(slice=='slice1_RB_p'){
    usy=E14y
    stg='E14_5'
  }
  if(slice=='slice1_4LB_p'){
    usy=Neoy
    stg='P3'
  }
  mergcoords=usy@images[[slice]]@coordinates
  coord1=unlist(mergcoords$row)
  coord2=unlist(mergcoords$col)
  
  for(gene in unique(c('Dbh','WPRE'))){
    if(gene %in% rownames(usy@assays$SCT@counts)){
      if(mergco){
        genefor1=log(usy@assays$SCT@counts[gene,which(colnames(usy@assays$SCT@counts)%in%rownames(mergcoords))]+1)
      }
      ploty=data.frame(Coord1=coord1,Coord2=coord2*-1,Gene1=genefor1)
      
      
      if(max(ploty$Gene1)!=min(ploty$Gene1[which(ploty$Gene1>0)])){
        png(paste0(OUTPUT_FILE_LOC,"/",gene,'_',stg,'_',i,'.png'),height=746,width=1156)  
        print(ggplot(as.data.frame(ploty),aes(Coord1, Coord2)) +
                geom_point(colour=c('grey'),size=1.5,show.legend = FALSE,alpha=c(0.6)) +
                geom_point(aes(color=Gene1),data=~subset(.,Gene1>0),size=3,alpha=1,show.legend=TRUE) +
                scale_color_gradientn(colours=c(viridis(11)[1],viridis(11)[11]),name=gene,labels=c('Low','High'),breaks=c(min(ploty$Gene1[which(ploty$Gene1>0)]),max(ploty$Gene1)))+
                labs(title=paste0(gene,' Log Expression ','Slice ',slice))+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),axis.title = element_blank(),
                      axis.text.x = element_blank(),axis.text.y = element_blank(),text=element_text(size=20),panel.background = element_blank(), legend.position = 'right',legend.key.size=unit(1,'cm'),legend.title=element_text(size=19,face='bold'),legend.text=element_text(size=15,face='bold'),axis.line = element_blank()))
        
        
        dev.off()
      }
    }
  }
}


#Now run for adult
slicelist=list()
for(sliceo in c('C5')){
  slicelist[[sliceo]]=readRDS(paste0(INT_FILE_LOC,"/Adult_",sliceo,"_bin_seurat.rds"))
  mergco=T
  mergcoords=slicelist[[sliceo]]@images[[names(slicelist[[sliceo]]@images)[1]]]@coordinates
  coord1=unlist(mergcoords$row)
  coord2=unlist(mergcoords$col)
  
  
  for(gene in unique(c('Dbh','WPRE'))){
    if(gene %in% rownames(slicelist[[sliceo]]@assays$SCT@counts)){
      if(mergco){
        genefor1=log(slicelist[[sliceo]]@assays$SCT@counts[gene,]+1)
      }
      ploty=data.frame(Coord1=coord1,Coord2=coord2*-1,Gene1=genefor1)
      
      
      if(max(ploty$Gene1)!=min(ploty$Gene1[which(ploty$Gene1>0)])){
        png(paste0(OUTPUT_FILE_LOC,"/Adult_",sliceo,'_',gene,'.png'),height=746,width=1156)  
        print(ggplot(as.data.frame(ploty),aes(Coord1, Coord2)) +
                geom_point(colour=c('grey'),size=1.5,show.legend = FALSE,alpha=c(0.6)) +
                geom_point(aes(color=Gene1),data=~subset(.,Gene1>0),size=3,alpha=1,show.legend=TRUE) +
                scale_color_gradientn(colours=c(viridis(11)[1],viridis(11)[11]),name=gene,labels=c('Low','High'),breaks=c(min(ploty$Gene1[which(ploty$Gene1>0)]),max(ploty$Gene1)))+
                labs(title=paste0(gene,' Log Expression ','Slice ',sliceo))+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),axis.title = element_blank(),
                      axis.text.x = element_blank(),axis.text.y = element_blank(),text=element_text(size=20),panel.background = element_blank(), legend.position = 'right',legend.key.size=unit(1,'cm'),legend.title=element_text(size=19,face='bold'),legend.text=element_text(size=15,face='bold'),axis.line = element_blank()))
        
        
        dev.off()
      }
    }
  }
}


#Figure S16
#RCTD clustering
#RCTD confusion matrix
#RCTD score



#Then BGI ran RCTD

library(RCTD)
E12_RCTD=readRDS(paste0(INT_FILE_LOC,"/E12.5_RCTD.rds"))
E14_RCTD=readRDS(paste0(INT_FILE_LOC,"/E14.5_RCTD.rds"))
Neo_RCTD=readRDS(paste0(INT_FILE_LOC,"/Neo_RCTD.rds")

E12_RCTD@results$results_df$first_type

E12_RCTD@spatialRNA

E12y@meta.data$RCTDd=E12_RCTD@results$results_df$first_type
E14y@meta.data$RCTDd=E14_RCTD@results$results_df$first_type
Neoy@meta.data$RCTDd=Neo_RCTD@results$results_df$first_type



coloursz=c('royalblue3',
'darkturquoise',
'springgreen4',
'lightskyblue3',
'lightyellow3',
'burlywood2',
'olivedrab4',
'darkslateblue',
'lightsteelblue4',
'cyan4',
'orchid4',
'yellow',
'firebrick2',
'tan2',
'navyblue',
'chartreuse1',
'lightpink3',
'maroon2',
'salmon4')

names(coloursz)=c('AM','AM-CCS','AVN','eAM','Endocardium','Endothelium','Epicardium','Erythroid Cells','eVM','eVM-trab','Fibroblast-like','Immune Cells','Neural Crest','PKJ','Platelets','SAN','Smooth Muscle-like','VM','VM-trab')


colalt=c('royalblue3','tomato4','pink3','maroon2','lightsteelblue4','springgreen4','yellowgreen','goldenrod3','salmon4','lightskyblue3','darkkhaki','goldenrod1','chartreuse1','aquamarine','honeydew3')
coloursz['eAM']='moccasin'
coloursz['Endocardium']='darkolivegreen4'
coloursz['Endothelium']='lightcyan3'
coloursz['PKJ']='darkorange3'


colsz=coloursz[match(levels(E12y@meta.data$RCTDd),names(coloursz))]

Idents(E12y)='RCTDd'
SpatialDimPlot(E12y,cols = c(colsz),images='slice1_A1_R',pt.size.factor=1.7,stroke=0.001)+NoLegend()
Idents(E14y)='RCTDd'
SpatialDimPlot(E14y,cols = c(colsz),images='slice1_LU_p',pt.size.factor=1.7,stroke=0.001)+NoLegend()
Idents(Neoy)='RCTDd'
SpatialDimPlot(Neoy,cols = c(colsz),images='slice1_4LB_p',pt.size.factor=1.7,stroke=0.001)+NoLegend()


Neoy@meta.data$pkjscore=Neo_RCTD@results$weights[,'PKJ']
Neoy@meta.data$sanscore=Neo_RCTD@results$weights[,'SAN']
Neoy@meta.data$avnscore=Neo_RCTD@results$weights[,'AVN']


for(j in c('pkjscore','sanscore','avnscore')){
  print(SpatialFeaturePlot(Neoy,images='slice1_4LU_p',features=c(j),pt.size.factor=1.7,stroke=c(0),alpha=c(0.1,1))+scale_fill_gradient(low='white',high='red'))
}


####Table S1
#Clean up merged_all_res2
pbmc=readRDS(paste0(INT_FILE_LOC,"/merged_all_res2_clustered_dim36.csv"))

#Table S1


chango=as.character(cmlin@meta.data$newtypez)
chango[which(chango=='eVM-trab')]='Early Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='AHF')]='Heart Fields'
chango[which(chango=='FHF')]='Heart Fields'
chango[which(chango=='pSHF')]='Heart Fields'
chango[which(chango=='VM')]='Ventricular Cardiomyocytes'
chango[which(chango=='AM')]='Atrial Cardiomyocytes'
chango[which(chango=='eVM')]='Immature Ventricular Cardiomyocytes'
chango[which(chango=='AVN')]='Atrioventricular Node'
chango[which(chango=='PKJ')]='Purkinje Fibres'
chango[which(chango=='pHT')]='Primary Heart Tube'
chango[which(chango=='VM-trab')]='Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='eAM')]='Immature Atrial Cardiomyocytes'
chango[which(chango=='EDCM')]='Endocardial gene-rich Cardiomyocytes'
chango[which(chango=='DevCM')]='Developing Cardiomyocytes'
chango[which(chango=='SAN')]='Sinoatrial Node'
chango[which(chango=='AM-CCS')]='Atrial Conduction System'
chango[which(chango=='Apop')]='Misc.'
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                              'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                              'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                              'Sinoatrial Node','Atrial Conduction System','Misc.'))

cmlin@meta.data$chango=chango
Idents(cmlin)='chango'
subbed=seq(0,14)
names(subbed)=levels(cmlin@meta.data$chango)
cmlin=RenameIdents(cmlin,subbed)
cmlin@meta.data$multimark=cmlin@active.ident
pbmc.multimark=FindAllMarkers.multicore(cmlin,min_pct=0.15,logfc_threshold=0.2,only_pos=FALSE,resolution='multimark',nCores=25)
saveRDS(pbmc.multimark,paste0(OUTPUT_FILE_LOC,"/CMlineage_DE_015_02.csv"))
for(i in 1:length(pbmc.multimark)){
  write.csv(pbmc.multimark[[i]][pbmc.multimark[[i]]$p_val_adj<=0.05,][order(pbmc.multimark[[i]][pbmc.multimark[[i]]$p_val_adj<=0.05,]$avg_log2FC,decreasing=TRUE),],paste0("F:/ML/Revision_NC/TableS1/CMlineage_DE_015_02_excel_",levels(chango)[i],".csv"))
}


##### Supplementary Data
#Data S1
corr_stg=pbmc@meta.data$orig.ident
corr_stg[which(corr_stg=='E8')]='E8.5'
corr_stg[which(corr_stg=='E10')]='E10.5'
corr_stg[which(corr_stg=='E12')]='E12.5'
corr_stg[which(corr_stg=='E14')]='E14.5'
corr_stg[which(corr_stg=='E16')]='E16.5'
corr_stg[which(corr_stg=='Neo')]='P3'
corr_stg=factor(corr_stg,levels=c('E8.5','E10.5','E12.5','E14.5','E16.5','P3'))
niceones=c('steelblue1','darkgoldenrod2','seagreen3','dodgerblue3','darkorchid1','royalblue3','cornflowerblue','brown2','blueviolet','lightyellow3','lightsalmon','burlywood2','olivedrab4','orchid4','darkslateblue','yellow','mediumslateblue','honeydew4','gold2','firebrick2','orange4','navyblue','cyan1','lightpink3','cadetblue','blue3')
colordf=data.frame(colours=niceones,row.names=levels(pbmc@meta.data$updated8))
figclust=plot_ly(as.data.frame(pbmc@reductions$umap@cell.embeddings), type="scatter3d", mode="markers", x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color =pbmc@meta.data$updated8 , colors = colordf[which(rownames(colordf)%in%unique(pbmc@meta.data$updated8)),1],
                 marker = list(symbol = 'circle', sizemode = 'diameter'),alpha=1, sizes = c(1),
                 text = ~paste('<br>Cell identity:',pbmc@meta.data$updated8,'<br>Stage:', corr_stg, '<br>nUMI:',pbmc@meta.data$nUMI,'<br>nGene:',pbmc@meta.data$nGene,'<br>MitoRatio:',round(pbmc@meta.data$mitoRatio,3)))
htmlwidgets::saveWidget(as_widget(figclust), paste0(OUTPUT_FILE_LOC,"/DataS1.html"))


#Data S2
cmlin=readRDS(paste0(INT_FILE_LOC,"/CMlineage_denovo.csv"))
phatecm=readRDS(paste0(INT_FILE_LOC,"/CMlineagephateknn12t14npca96gamma0_8.csv"))

corr_stg=cmlin@meta.data$orig.ident
corr_stg[which(corr_stg=='E8')]='E8.5'
corr_stg[which(corr_stg=='E10')]='E10.5'
corr_stg[which(corr_stg=='E12')]='E12.5'
corr_stg[which(corr_stg=='E14')]='E14.5'
corr_stg[which(corr_stg=='E16')]='E16.5'
corr_stg[which(corr_stg=='Neo')]='P3'
corr_stg=factor(corr_stg,levels=c('E8.5','E10.5','E12.5','E14.5','E16.5','P3'))


chango=as.character(cmlin@meta.data$newtypez)
chango[which(chango=='eVM-trab')]='Early Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='AHF')]='Heart Fields'
chango[which(chango=='FHF')]='Heart Fields'
chango[which(chango=='pSHF')]='Heart Fields'
chango[which(chango=='VM')]='Ventricular Cardiomyocytes'
chango[which(chango=='AM')]='Atrial Cardiomyocytes'
chango[which(chango=='eVM')]='Immature Ventricular Cardiomyocytes'
chango[which(chango=='AVN')]='Atrioventricular Node'
chango[which(chango=='PKJ')]='Purkinje Fibres'
chango[which(chango=='pHT')]='Primary Heart Tube'
chango[which(chango=='VM-trab')]='Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='eAM')]='Immature Atrial Cardiomyocytes'
chango[which(chango=='EDCM')]='Endocardial gene-rich Cardiomyocytes'
chango[which(chango=='DevCM')]='Developing Cardiomyocytes'
chango[which(chango=='SAN')]='Sinoatrial Node'
chango[which(chango=='AM-CCS')]='Atrial Conduction System'
chango[which(chango=='Apop')]='Misc.'
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                              'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                              'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                              'Sinoatrial Node','Atrial Conduction System','Misc.'))


colalt=c('royalblue3','tomato4','pink3','maroon2','lightsteelblue4','springgreen4','yellowgreen','goldenrod3','salmon4','lightskyblue3','darkkhaki','goldenrod1','chartreuse1','aquamarine','honeydew3')
clustcolours2=data.frame(sort(unique(chango)),colalt)
rownames(clustcolours2)=clustcolours2[,1]

figclust=plot_ly(as.data.frame(phatecm), type="scatter3d", mode="markers", x = ~PHATE1, y = ~PHATE2, z = ~PHATE3, color =as.factor(chango), colors = clustcolours2[which(unique(as.factor(chango))%in%clustcolours2[,1]),2],
                 marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1),
                 text = ~paste('<br>Cell identity:',chango,'<br>Stage:', corr_stg, '<br>nUMI:',cmlin@meta.data$nUMI,'<br>nGene:',cmlin@meta.data$nGene,'<br>MitoRatio:',round(cmlin@meta.data$mitoRatio,3)))
htmlwidgets::saveWidget(as_widget(figclust), paste0(OUTPUT_FILE_LOC,"/DataS2.html"))


corr_stg=Dbhp@meta.data$orig.ident
corr_stg[which(corr_stg=='E8')]='E8.5'
corr_stg[which(corr_stg=='E10')]='E10.5'
corr_stg[which(corr_stg=='E12')]='E12.5'
corr_stg[which(corr_stg=='E14')]='E14.5'
corr_stg[which(corr_stg=='E16')]='E16.5'
corr_stg[which(corr_stg=='Neo')]='P3'
corr_stg=factor(corr_stg,levels=c('E8.5','E10.5','E12.5','E14.5','E16.5','P3'))


chango=as.character(Dbhp@meta.data$newtypez)
chango[which(chango=='eVM-trab')]='Early Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='AHF')]='Heart Fields'
chango[which(chango=='FHF')]='Heart Fields'
chango[which(chango=='pSHF')]='Heart Fields'
chango[which(chango=='VM')]='Ventricular Cardiomyocytes'
chango[which(chango=='AM')]='Atrial Cardiomyocytes'
chango[which(chango=='eVM')]='Immature Ventricular Cardiomyocytes'
chango[which(chango=='AVN')]='Atrioventricular Node'
chango[which(chango=='PKJ')]='Purkinje Fibres'
chango[which(chango=='pHT')]='Primary Heart Tube'
chango[which(chango=='VM-trab')]='Trabecular Ventricular Cardiomyocytes'
chango[which(chango=='eAM')]='Immature Atrial Cardiomyocytes'
chango[which(chango=='EDCM')]='Endocardial gene-rich Cardiomyocytes'
chango[which(chango=='DevCM')]='Developing Cardiomyocytes'
chango[which(chango=='SAN')]='Sinoatrial Node'
chango[which(chango=='AM-CCS')]='Atrial Conduction System'
chango[which(chango=='Apop')]='Misc.'
chango=factor(chango,levels=c('Atrial Cardiomyocytes','Early Trabecular Ventricular Cardiomyocytes','Heart Fields','Ventricular Cardiomyocytes',
                              'Immature Ventricular Cardiomyocytes','Atrioventricular Node','Purkinje Fibres','Primary Heart Tube','Trabecular Ventricular Cardiomyocytes',
                              'Immature Atrial Cardiomyocytes','Endocardial gene-rich Cardiomyocytes','Developing Cardiomyocytes',
                              'Sinoatrial Node','Atrial Conduction System','Misc.'))




colalt=c('royalblue3','tomato4','pink3','maroon2','lightsteelblue4','springgreen4','yellowgreen','goldenrod3','salmon4','lightskyblue3','darkkhaki','goldenrod1','chartreuse1','aquamarine')
clustcolours=data.frame(sort(unique(chango)),colalt)
names(colalt)=clustcolours[,1]
figclust=plot_ly(as.data.frame(phatedbh$embedding), type="scatter3d", mode="markers", x = ~PHATE1, y = ~PHATE2, z = ~PHATE3, color =chango, colors = colalt,
                 marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(1),
                 text=~paste('<br>Cell identity:',Dbhp@meta.data$newtypez,'<br>Stage:', corr_stg, '<br>nUMI:',Dbhp@meta.data$nUMI,'<br>nGene:',Dbhp@meta.data$nGene,'<br>MitoRatio:',round(Dbhp@meta.data$mitoRatio,3)))
htmlwidgets::saveWidget(as_widget(figclust), paste0(OUTPUT_FILE_LOC,"/DataS3.html"))

