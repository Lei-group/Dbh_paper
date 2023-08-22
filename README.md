# Dbh_paper
This includes code and data relevant to produce the figures shown in the Lei group paper 2023 NC.

Plots for the main figures, supplementary figures, and additional diagnostic and plots of interest are produced throughout the code.

The code (mostly in R), starts by loading in the raw count data and performing preliminary QC and filtering on it.
Then the code includes our preliminary engagement with the data across examining batch effects, 'dropout' rates, and clustering and subclustering. The clustering and subclustering included a series of iterative re-use of the same code for different clusters which is why there are only certain clusters named in this portion of the code.
The code then extracts the cell types within the cardiomyocyte lineage and performs clustering and subclustering on these to classify these cells.
The code then varies between downstream analyses, plot production, and exporting data to Python.

