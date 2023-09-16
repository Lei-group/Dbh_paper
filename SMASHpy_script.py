# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Where intermediate files stored as per other script
INT_FILE_LOC=XXXX


import scanpy as sc


import datatable
import pandas as pd
import numpy as np
import anndata
import tensorflow as tf


gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.experimental.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)


from tensorflow.keras.optimizers import Adam
import smashpy
sm=smashpy.smashpy()

#Due to size, found loading in chunkwise was easiest in terms of stability

sms1=datatable.fread(INT_FILE_LOC+"/forSMASH_1_10k.csv")
sms1=sms1.to_pandas()

sms2=datatable.fread(INT_FILE_LOC+"/forSMASH_10001_20k.csv")
sms2=sms2.to_pandas()

sms3=datatable.fread(INT_FILE_LOC+"/forSMASH_20001_28k.csv")
sms3=sms3.to_pandas()


smshy=pd.concat([sms1,sms2],axis=1)
smshy=pd.concat([smshy,sms3],axis=1)

#read cell type classifications
clusts=pd.read_csv(INT_FILE_LOC+"/clusters.csv")

clusty=[None]*len(clusts.values)
for i in range(0,len(clusts.values)):
   clusty[i]=str((str(np.array(clusts)[i])[len(str(i))+4:len(str(np.array(clusts)[i]))-3]))

#If any problems with the clusters - try this. Sometimes it seems pd.read doesn't clear starting punct.
for i in range(0,len(clusty)):
    if clusty[i]=='"AM':
        clusty[i]='AM'
    if clusty[i]=='"VM-trab':
        clusty[i]='VM-trab'


adata=anndata.AnnData(np.array(smshy),obs=np.array(clusty),var=np.array(smshy.columns))
adata.var.index=adata.var[0]
########now ready to start with SMASHPY
adata.obs["annotation"]=np.array(clusty)

del sms1
del sms2
del sms3
del smshy

#list(set(adata.obs['annotation']))
objsav=adata
adata=objsav

sm.data_preparation(adata)


#These can be downloaded from SMASHpy github
adata=sm.remove_general_genes(adata,species="mouse",path_1=INT_FILE_LOC+"/housekeeping_fibrobblast_mus.txt",
							  path_2=INT_FILE_LOC+"/housekeeping_ventricles_mus.txt")

###need to remove housekeeping genes too

"""
#This function kept generating OOM errors, so we just manually run through it below
adata=sm.remove_features_pct(adata,group_by="annotation",pct=0.3)
"""

group_by="annotation"
pct=0.3

adata.X   = adata.layers["counts"]
adata.raw = adata.copy()
		
list_keep_genes = []
		
df = pd.DataFrame(data=False, 
				  index=adata.var.index.tolist(),
				  columns=set(adata.obs[group_by]))
for g in set(adata.obs[group_by]): 
	reduced = adata[adata.obs[group_by]==g]
	boolean, values = sc.pp.filter_genes(reduced, min_cells = reduced.n_obs*pct, inplace=False)
	df[g] = boolean
dfT = df.T
for g in dfT.columns:
    if True in dfT[g].tolist():
        list_keep_genes.append(True)
    else:
        list_keep_genes.append(False)
		
adata.var["general"] = list_keep_genes
		
adata = adata[:, adata.var["general"]]
		
adata.X   = adata.layers["log"]
adata.raw = adata.copy()
adata.X   = adata.layers["scale"]
		
del reduced
del df
del dfT


#Now run next filter step
adata=sm.remove_features_pct_2groups(adata,group_by="annotation",pct1=0.75,pct2=0.5)

#Inverse PCA filter
adata = sm.scale_filter_features(adata, n_components=None, filter_expression=True)

###### Ensemble learning method - this step takes a while

adata.var.index.name=None

# Top 20 genes as a final dictionary, for each annotation (class) provided
# Calculate the importances of each gene using the Gini impurity
#Produces several plots - including Figure S16 panel a - please note there is a degree of stochasity to ensemble learning methods so results may vary 
clf=sm.ensemble_learning(adata,group_by='annotation',classifier='XGBoost',balance=True,verbose=True)
selectedGenes, selectedGenes_dict = sm.gini_importance(adata, clf, group_by="annotation", verbose=True, restrict_top=("local", 100))
sm.run_classifiers(adata, group_by="annotation", genes=selectedGenes, classifier="KNN", balance=True, title="XGBoost-KNN") 
axs, selectedGenes_top_dict = sm.sort_and_plot(adata, selectedGenes, group_by="annotation", top=50, figsize=(13,50),restricted=False)


pd.DataFrame.from_dict(selectedGenes_top_dict).to_csv(INT_FILE_LOC+"/XGboost_KNN_50_out.csv")

listygenes=[]
for i in selectedGenes_top_dict:
	listygenes.extend(selectedGenes_top_dict[i])
listygenes=set(listygenes)

pd.DataFrame(listygenes).to_csv(INT_FILE_LOC+"/XGboost_KNN_50_out_as_gene_list.csv")
