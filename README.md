# Dbh_paper
This includes code and data relevant to produce the figures shown in the Lei group paper 2023 NC.

Plots for the main figures, supplementary figures, and additional diagnostic and plots of interest are produced throughout the code. However, we also produce distinct scripts that collate this into a Figure-specific basis for ease of reproduction. Panels were using these scripts below, with additional annotation or arrangement manually.

#Pre_process.R
The code (mostly in R), starts by loading in the raw count data and performing preliminary QC and filtering on it.
Then the code includes our preliminary engagement with the data across examining batch effects, 'dropout' rates, and clustering and subclustering.
The code then extracts the cell types within the cardiomyocyte lineage and performs clustering and subclustering on these to classify these cells.
The clustering and subclustering included a series of iterative re-use of the same code for different clusters thus this portion of the code is representative not exhaustive.
The final objects are produced within this script, for use for downstream analysis.

The code then varies between downstream analyses, plot production, and exporting data to Python.
#Figure1.R
  #Associated supplementary figures for Figure 1
  #Figure S2
  #Figure S3
  #Figure S5
  #Figure S6
  #Figure S8
  #Figure S9
  #Figure S10
  
#Figure2.R
  #Associated supplementary figures for Figure 2
  #Figure S13
  #Figure S14
  #Figure S15
  #Figure S16

Figure S16 involved work in python, which we include here.
#FigureS16.py
  

