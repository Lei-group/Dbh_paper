####This is based on code from https://github.com/vertesy with minor modifications such as specifying the 'bimod' option by default etc


##############VERTSEY MULTICORE SEURAT
#install.packages('tictoc')
#install.packages("doSNOW")

#install.packages("doParallel") 

#BiocManager::install('doMC')
#install.packages("doMPI")
library(tictoc)
library(doSNOW)
library(doParallel)
library(doMPI)
library(doMC)
######################################################################
# Seurat2.Multicore.Functions.R
######################################################################
# source('~/GitHub/Packages/Seurat.multicore/Seurat2.Multicore.Functions.R')
# Source: self + web

try(require(Seurat), silent = F)
try(require(doMC), silent = F)


# FindAllMarkers.multicore ------------------------------

FindAllMarkers.multicore <- function(obj = org, min_pct = 0.2, logfc_threshold=0.5, only_pos=F, wait=10, resolution=Clust_by, nCores =6){ # Multicore version of FindAllMarkers.
  tictoc::tic()
  nrClusters=length(unique(obj@meta.data[,resolution]))
  N=nrClusters-1
  
  j=seq(from = 0, by = wait, length.out = nCores)
  j=rep(j,nrClusters)
  ls.DE <- foreach(i=0:N) %dopar% {
    Sys.sleep(j[i+1])
    FindMarkers(object = obj,slot="data", ident.1=i, only.pos = only_pos, min.pct=min_pct, test.use='bimod', logfc.threshold = logfc_threshold, random.seed=804L)
  };
  tictoc::toc()
  return(ls.DE)
}



# Work in progress ------------------------------------------------------------




# ------------------------------------------------------------------------------------------
## Seurat.Functions.other.R
# ------------------------------------------------------------------------------------------
# Source: self + web
# source("~/GitHub/Packages/Seurat.multicore/Seurat.Functions.other.R")

# ------------------------------


MergeDuplicateGenesSeurat <- function(seu=ls.Seurat[[i]]){ # How many entries are duplicated
  duplicated(rownames(seu))
  if (summarize & y){
    x = table(vec); x= x[x>1]-1;
    print("The following elements have >1 extra copies:")
    print(x) # table formatting requires a separate entry
  }
  return(y)
}


# quick umap ---------------
umap <- function(gene='DLX5', obj =org, pt_size =1) {
  FeaturePlot(object = obj, features.plot = gene, reduction.use = 'umap', pt.size = pt_size)
}

# FeaturePlot with different defaults ------------------------------------------------------------------
aFeaturePlot <- function(object=org, features.plot, min.cutoff = 'q1', max.cutoff = 'q99',
                         dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1
                         , cols.use = c("yellow", "red"), pch.use = 16, overlay = FALSE, do.hover = FALSE,
                         data.hover = "ident", do.identify = FALSE, reduction.use = "umap",
                         use.imputed = FALSE, nCol = NULL, no.axes = T, no.legend = F,
                         coord.fixed = FALSE, dark.theme = FALSE, do.return = FALSE,
                         vector.friendly = TRUE, png.file = NULL, png.arguments = c(10, 10, 100))
{
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  if (is.null(x = nCol)) {
    nCol <- 2
    if (length(x = features.plot) == 1) {
      nCol <- 1
    }
    if (length(x = features.plot) > 6) {
      nCol <- 3
    }
    if (length(x = features.plot) > 9) {
      nCol <- 4
    }
  }
  num.row <- floor(x = length(x = features.plot)/nCol - 1e-05) +
    1
  if (overlay | do.hover) {
    num.row <- 1
    nCol <- 1
  }
  par(mfrow = c(num.row, nCol))
  dim.code <- GetDimReduction(object = object, reduction.type = reduction.use,
                              slot = "key")
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- as.data.frame(GetCellEmbeddings(object = object,
                                               reduction.type = reduction.use, dims.use = c(dim.1, dim.2),
                                               cells.use = cells.use))
  x1 <- paste0(dim.code, dim.1)
  x2 <- paste0(dim.code, dim.2)
  data.plot$x <- data.plot[, x1]
  data.plot$y <- data.plot[, x2]
  data.plot$pt.size <- pt.size
  names(x = data.plot) <- c("x", "y")
  data.use <- t(x = FetchData(object = object, vars.all = features.plot,
                              cells.use = cells.use, use.imputed = use.imputed))
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    ifelse(test = is.na(x = cutoff), yes = min(data.use[feature,
    ]), no = cutoff)
  }, cutoff = min.cutoff, feature = features.plot)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    ifelse(test = is.na(x = cutoff), yes = max(data.use[feature,
    ]), no = cutoff)
  }, cutoff = max.cutoff, feature = features.plot)
  check_lengths = unique(x = vapply(X = list(features.plot,
                                             min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check_lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  if (overlay) {
    pList <- list(BlendPlot(data.use = data.use, features.plot = features.plot,
                            data.plot = data.plot, pt.size = pt.size, pch.use = pch.use,
                            cols.use = cols.use, dim.codes = dim.codes, min.cutoff = min.cutoff,
                            max.cutoff = max.cutoff, coord.fixed = coord.fixed,
                            no.axes = no.axes, no.legend = no.legend, dark.theme = dark.theme))
  }
  else {
    pList <- mapply(FUN = SingleFeaturePlot, feature = features.plot,
                    min.cutoff = min.cutoff, max.cutoff = max.cutoff,
                    coord.fixed = coord.fixed, MoreArgs = list(data.use = data.use,
                                                               data.plot = data.plot, pt.size = pt.size, pch.use = pch.use,
                                                               cols.use = cols.use, dim.codes = dim.codes, no.axes = no.axes,
                                                               no.legend = no.legend, dark.theme = dark.theme,
                                                               vector.friendly = vector.friendly, png.file = png.file,
                                                               png.arguments = png.arguments), SIMPLIFY = FALSE)
  }
  if (do.hover) {
    if (length(x = pList) != 1) {
      stop("'do.hover' only works on a single feature or an overlayed FeaturePlot")
    }
    if (is.null(x = data.hover)) {
      features.info <- NULL
    }
    else {
      features.info <- FetchData(object = object, vars.all = data.hover)
    }
    return(HoverLocator(plot = pList[[1]], data.plot = data.plot,
                        features.info = features.info, dark.theme = dark.theme,
                        title = features.plot))
  }
  else if (do.identify) {
    if (length(x = pList) != 1) {
      stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
    }
    return(FeatureLocator(plot = pList[[1]], data.plot = data.plot,
                          dark.theme = dark.theme))
  }
  else {
    print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
  }
  ResetPar()
  if (do.return) {
    return(pList)
  }
}



# Save multiple FeaturePlot from a list of genes on A4 jpeg ------------------------
# "Does not work!"
# amultiFeaturePlot.A4 <- function(list.of.genes, object = org, plot.reduction='umap'
#                                  , colors=c("grey", "red"), nr.Col=2, nr.Row =4, cex = ceiling(10/(nr.Col*nr.Row))
#                                  , gene.min.exp = 'q01', gene.max.exp = 'q99'
#                                  , jpeg.res = 225, jpeg.q = 90) {
#   tictoc::tic()
#   list.of.genes = check.genes(list.of.genes, obj = object)
#   lsG = iterBy.over(1:length(list.of.genes), by=nr.Row*nr.Col)
#
#   ls.plots <- foreach(i=1:length(lsG)) %dopar% {
#     genes = list.of.genes[lsG[[i]]]
#     FeaturePlot(object, features.plot =genes, reduction.use = plot.reduction
#                 , nCol = nr.Col, cols.use = colors, no.axes = T, no.legend = F, vector.friendly = T
#                 , min.cutoff = gene.min.exp, max.cutoff = gene.max.exp, do.return = T
#                 , pt.size = cex)
#   }
#   print ("Plots are calculated.")
#   for (i in 1:length(lsG)) { print(i )
#     plotname = kpp(c(plot.reduction,i, genes, 'jpg' ))
#     jjpegA4(plotname, r = jpeg.res, q = jpeg.q)
#     print(ls.plots[[i]])
#     try.dev.off()
#   }
#   tictoc::toc()
# };
# "Does not work!"



# LabelPoint ------------------------
LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0,
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t,
                            label = i, size = text.size)
    plot <- plot + annotate("segment",
                            x = x1 + adj.x.s,
                            xend = x1,
                            y = y1 +  adj.y.s,
                            yend = y1,
                            size = segment.size)
  }
  return(plot)
}

#  ------------------------
LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05,
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(
    plot,
    genes,
    exp.mat,
    adj.y.t = adj.u.t,
    adj.x.t = adj.r.t,
    adj.y.s = adj.u.s,
    adj.x.s = adj.r.s,
    ...
  ))
}

#  ------------------------
LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05,
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(
    plot,
    genes,
    exp.mat,
    adj.y.t = adj.u.t,
    adj.x.t = -adj.l.t,
    adj.y.s = adj.u.s,
    adj.x.s = -adj.l.s,
    ...
  ))
}

#  ------------------------
LabelBR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05,
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(
    plot,
    genes,
    exp.mat,
    adj.y.t = -adj.u.t,
    adj.x.t = adj.l.t,
    adj.y.s = adj.u.s,
    adj.x.s = adj.l.s,
    ...
  ))
}


#  ------------------------
LabelBL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05,
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(
    plot,
    genes,
    exp.mat,
    adj.y.t = -adj.u.t,
    adj.x.t = -adj.l.t,
    adj.y.s = adj.u.s,
    adj.x.s = -adj.l.s,
    ...
  ))
}


# ------------------------------------------------------------

if (F) {
  "Very slow for some reason"
  # extended Seurat object
  setClass(
    "xseurat",
    contains="seurat",
    slots=c(p="list",
            all.genes = "list")
  ) -> xseurat
  # x <- as(org.discarded, "xseurat")
  
  
  extended.Seurat.class <- function(obj) {
    as(obj, "xseurat")
  }
  
  xz = extended.Seurat.class(org)
  is(xz)
  x@p = p
  x@all.genes = all.genes
  is(x)
  
} # if

# ------------------------------------------------------------
# check.genes <- function(list.of.genes = ClassicMarkers, obj = seu3) { # check if genes exist in your dataset
#   missingGenes = setdiff(list.of.genes, rownames(obj))
#   if (length(missingGenes)>0) {iprint("Genes not found in the data:", missingGenes)}
#   intersect(list.of.genes, rownames(obj))
# }



# ------------------------------------------------------------

# ------------------------------------------------------------

# ------------------------------------------------------------

# ------------------------------------------------------------
######################################################################
# Seurat3.Multicore.Generic.Functions.R
######################################################################
# source('~/GitHub/Packages/Seurat.multicore/Seurat3.Multicore.Generic.Functions.R')

# NOTE:
# Seurat v3 uses the 'future' framework for parallelization.
# https://satijalab.org/seurat/v3.0/future_vignette.html

# These provide an alternative way, advantages/disadvantages are untested.

try(require(Seurat), silent = F)
try(require(doMC), silent = F)
# try(require(future), silent = T)

# try(source("~/GitHub/pseudoBulk/barcode.export.from.Seurat/Seurat3.Write.Out.CBCs.for.subset-bam.R"), silent = T)

# ------------------------------------------------------------------------

parallel.computing.by.future <- function(workers_ = 6, maxMemSize = 4000 * 1024^2) { # Run gc(), load multi-session computing and extend memory limits.
  # https://satijalab.org/seurat/v3.0/future_vignette.html
  cat(
    "1. If you load futures before you finished using foreach loops,
    NormalizeData inside a foreach loop fails (Error writing to connection)
    -> I assume 'future' and 'doMC' are not compatible
    2. If you setup computing on e.g. six cores, it runs 6 instances of R with the entire memory space copied.
    If you run out of memory, the system starts using the SSD as memory, and it slows you down extremely extremely extremely.
    -> Therefore it is important to clean up the memory space before setting up multicore computation.
    Loaded: library(future), workers set to 6 (def),set Max mem size to 2GB (def)."   )
  
  gc(full = T)
  try(memory.biggest.objects())
  
  library(future)
  plan("multiprocess", workers = workers_)
  # plan("multisession", workers = workers_)
  # So to set Max mem size to 2GB, you would run :
  options(future.globals.maxSize = maxMemSize)
}

######################################################################
# Seurat3.Multicore.Read.Write.R
######################################################################
# source('~/GitHub/Seurat.multicore/Seurat3.Multicore.Read.Write.R')

"Multicore read / write (I/O) functions are https://github.com/vertesy/Seurat.multicore"
"Single core read / write (I/O) functions are in https://github.com/vertesy/Seurat.utils/"

try(require(MarkdownReportsDev))
try(require(tictoc))
try(require(readr))
try (source ('https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R'),silent= F)

# Save an object -----------------------------------------------
isave.RDS.pigz <- function(object, prefix =NULL, suffix=NULL, showMemObject=T, saveParams =T){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  dir.create(path_rdata)
  
  if ( "seurat" %in% is(object) & saveParams) {
    try(object@misc$p <- p, silent = T)
    try(object@misc$all.genes  <- all.genes, silent = T)
  }
  if (showMemObject) { memory.biggest.objects() }
  fnameBase = kppu(prefix, substitute(object), suffix, idate())
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",fnameBase , ".Rds")
  tictoc::tic()
  ssaveRDS(object = object, file = fname)
  tictoc::toc()
  say()
}

# seuSaveRds ------------------------------------------------------------------------

seuSaveRds <- function(object = ls.Seurat, tags = setupFlags, use_Original_OutDir = F) { # Save a compressed Seurat Object, with parallel gzip by pgzip
  if (use_Original_OutDir) create_set_Original_OutDir()
  fname.comb.rds = ppp(substitute(object), tags, idate(), ".Rds")
  iprint( fname.comb.rds)
  ssaveRDS(object = object, filename = fname.comb.rds)
  say()
}

# Wrapper layer 2 (top) -----------------------------------------------------------------------------------

# read / load multiple objects
rrRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org"), ...) { # Load a list of RDS files with parallel ungzip by pgzip.
  tictoc::tic()
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  iprint("Looking for files under:" , path_rdata)
  iprint("(Base on OutDir)")
  dir.exists(path_rdata)
  # lsX <- foreach(obj = list_of_objectnames) %dopar% {
  for (obj in list_of_objectnames) {
    fname = MarkdownReportsDev::kollapse(path_rdata , "/", obj)
    tmp.obj = rreadRDS(filename = fname, ..., cores=2)
    # pigzx multi-core decompression is not so useful, with 2 cores ist actually faster:
    
    assigned.name = strsplit(obj, split = '\\.201*')[[1]][1]
    print(paste(" --> Loaded as object: ",assigned.name))
    assign(x = assigned.name,value = tmp.obj, envir = .GlobalEnv)
  }
  tictoc::toc()
}

# require(clipr)
# clip2clip.vector()
# x = c("org.2019.03.12_10h.Rds", "ls2.2019.03.12_10h.Rds", "ls.Seurat.2019.03.12_10h.Rds")
# rrRDS(x)

# Save multiple objects using pigz by default ---------------------------------------------
sssRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org.ALL", "org"), name.suffix =NULL, ...) { # Save multiple objects into a list of RDS files using parallel gzip by pgzip (optional).
  tictoc::tic()
  base_name <- character()
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  iprint("Files saved under:" , path_rdata)
  
  try(dir.create(path_rdata))
  for (obj in list_of_objectnames) {
    print(obj)
    if ( "seurat" %in% is(obj)) {
      obj@misc$p = p
      obj@misc$all.genes = all.genes
    }
    base_name[i] = paste0(obj, '.', name.suffix, '.', idate(),".Rds", collapse = "")
    fname = paste0(path_rdata , "/", base_name[i], collapse = "")
    ssaveRDS( object = get(obj), filename = fname, ...)
  }
  tictoc::toc()
  dput(base_name)
}

### Wrapper layer 1 -----------------------------------------------------------------------------------
ssaveRDS <- function(object, filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], func_type = c("pipe", "builtin")[1], ...) { # Save an object with parallel gzip by pgzip.
  tictoc::tic()
  if (func_type == "builtin") {
    con <- con_func(filename)
  } else {
    con <- con_func(filename, mode="write")
  }
  on.exit(close(con))
  saveRDS(object, con)
  tictoc::toc()
}

### rreadRDS -----------------------------------------------------------------------------------
rreadRDS <- function(filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], func_type = c("pipe", "builtin")[1], ...) {  # Read an object with parallel ungzip by pgzip.
  tictoc::tic()
  if (func_type == "builtin") {
    con <- con_func(filename)
  } else {
    con <- con_func(filename, mode="read", ...)
  }
  on.exit(close(con))
  obj=readr::read_rds(con)
  tictoc::toc()
  return(obj)
}

# snappy_pipe -----------------------------------------------------------------------------------
snappy_pipe <- function(filename, mode="read") { # Alternative, fast compression. Low compression rate, lightning fast.
  if (mode == "read") {
    con <- pipe(paste0("cat ", filename, " | snzip -dc"), "rb")
  } else {
    con <- pipe(paste0("snzip > ", filename), "wb")
  }
  con
}

# pigz_pipe -----------------------------------------------------------------------------------
pigz_pipe <- function(filename, mode="read", cores=6) { # Alternative: normal gzip output (& compression rate), ~*cores faster in zipping.
  if (mode == "read") {
    con <- pipe(paste0("cat ", filename, " | pigz -dcp ", cores), "rb")
  } else {
    con <- pipe(paste0("pigz -p ", cores, " > ", filename), "wb")
  }
  con
}
