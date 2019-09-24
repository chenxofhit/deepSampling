rm(list=ls())
library(SingleCellExperiment)
library(rhdf5)

setwd('/Users/chenx/Documents/Exp/DeepSampling/')

# Import necessary libraries and functions required for preporcessing
suppressMessages(source("./pre-processing/functions.R"))
suppressMessages(source("./pre-processing/libraries.R"))

#process 1: deal  with raw data
RAW_DIR <- file.path(getwd(), "raw/")       # SPECIFY HERE
DATA_DIR  <- file.path(getwd(),"data/")

#worm_neuron_cell.h5
h5name <- 'worm_neuron_cell.h5'
exprfilename <- paste0(h5name, ".expr.csv")
annotationfilename <- paste0(h5name, ".label.csv")

sce.y <- h5read(file.path(RAW_DIR, h5name), "Y")
sce.x <- h5read(file.path(RAW_DIR, h5name), "X")

exprfile <- file.path(DATA_DIR, exprfilename)
annotationfile <- file.path(DATA_DIR, annotationfilename)

texpr.logcounts <-  t(as.matrix((sce.x + 1)))

# Pre-process
l <- normalize_by_umi_2(texpr.logcounts)
texpr.logcounts.processed <- matrix.subset(l, 3000, log=FALSE)

write.csv(as.matrix(texpr.logcounts.processed), file=exprfile, sep = ",", row.names=FALSE, col.names=FALSE)
write.csv(sce.y, file=annotationfile, sep = ",", row.names=FALSE, col.names=FALSE)

print(paste0(h5name, ' write finished!'))

#10X_PBMC.h5
h5name <- '10X_PBMC.h5'
exprfilename <- paste0(h5name, ".expr.csv")
annotationfilename <- paste0(h5name, ".label.csv")

sce.y <- h5read(file.path(RAW_DIR, h5name), "Y")
sce.x <- h5read(file.path(RAW_DIR, h5name), "X")

exprfile <- file.path(DATA_DIR, exprfilename)
annotationfile <- file.path(DATA_DIR, annotationfilename)

texpr.logcounts <-  t(as.matrix((sce.x + 1)))

# Pre-process
l <- normalize_by_umi_2(texpr.logcounts)
texpr.logcounts.processed <- matrix.subset(l, 3000, log=FALSE)

write.csv(as.matrix(texpr.logcounts.processed), file=exprfile, sep = ",", row.names=FALSE, col.names=FALSE)
write.csv(sce.y, file=annotationfile, sep = ",", row.names=FALSE, col.names=FALSE)

print(paste0(h5name, ' write finished!'))


#mouse_ES_cell.h5
h5name <- 'mouse_ES_cell.h5'
exprfilename <- paste0(h5name, ".expr.csv")
annotationfilename <- paste0(h5name, ".label.csv")

sce.y <- h5read(file.path(RAW_DIR, h5name), "Y")
sce.x <- h5read(file.path(RAW_DIR, h5name), "X")

exprfile <- file.path(DATA_DIR, exprfilename)
annotationfile <- file.path(DATA_DIR, annotationfilename)

texpr.logcounts <-  t(as.matrix((sce.x + 1)))

# Pre-process
l <- normalize_by_umi_2(texpr.logcounts)
texpr.logcounts.processed <- matrix.subset(l, 3000, log=FALSE)

write.csv(as.matrix(texpr.logcounts.processed), file=exprfile, sep = ",", row.names=FALSE, col.names=FALSE)
write.csv(sce.y, file=annotationfile, sep = ",", row.names=FALSE, col.names=FALSE)

print(paste0(h5name, ' write finished!'))



#mouse_bladder_cell.h5
h5name <- 'mouse_bladder_cell.h5'
exprfilename <- paste0(h5name, ".expr.csv")
annotationfilename <- paste0(h5name, ".label.csv")

sce.y <- h5read(file.path(RAW_DIR, h5name), "Y")
sce.x <- h5read(file.path(RAW_DIR, h5name), "X")

exprfile <- file.path(DATA_DIR, exprfilename)
annotationfile <- file.path(DATA_DIR, annotationfilename)

texpr.logcounts <-  t(as.matrix((sce.x + 1)))

# Pre-process
l <- normalize_by_umi_2(texpr.logcounts)
texpr.logcounts.processed <- matrix.subset(l, 3000, log=FALSE)

write.csv(as.matrix(texpr.logcounts.processed), file=exprfile, sep = ",", row.names=FALSE, col.names=FALSE)
write.csv(sce.y, file=annotationfile, sep = ",", row.names=FALSE, col.names=FALSE)

print(paste0(h5name, ' write finished!'))
