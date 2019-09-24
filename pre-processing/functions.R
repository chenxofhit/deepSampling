# ----------------------------------------------------------------------
# Normalization by library size
# ----------------------------------------------------------------------

normalize_by_umi_2 <-function(x, min.count=2, min.cell=3,verbose=TRUE) {
  tsmessage("[!normalization] Onset...", verbose=verbose)
  mat  = x
  gene_symbols = colnames(x)
  cs <- colSums(mat>min.count)
  x_use_genes <- which(cs > min.cell)
  
  x_filt<-mat[,x_use_genes]
  
  gene_symbols = gene_symbols[x_use_genes]

  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  tsmessage("[!normalization] Complete.", verbose=verbose)
  
  list(m=x_norm,use_genes=gene_symbols)
}

normalize_by_umi <- function(M,doc_proportion_max = 1,doc_proportion_min = 0.05, normalizeCounts=FALSE,verbose=TRUE){

  ix = Matrix::rowSums(M!=0)
  M = M[ix>ncol(M)*doc_proportion_min & ix<=ncol(M)*doc_proportion_max,]
  gficf <- tf(M, verbose = verbose)
  w <- getIdfW(gficf, verbose = verbose)
  gficf.new <- idf(gficf, w, verbose = verbose)
  gficf.norm <- l.norm(t(gficf.new), norm = "l2", verbose = verbose)
  gficf.final <- scaleMatrix(gficf.norm, rescale = T, centre = T, verbose = verbose)
  
  print(paste("[!info] Log normalized matrix dimensions :", dim(gficf.final)[1], " x ", dim(gficf.final)[2]))
  print("[!log-transform] Complete.")
  
  return(gficf.final)
}

# ----------------------------------------------------------------------
# Dispersion Genes and umi_2bsetting (Gene-selection and log transform)
# ----------------------------------------------------------------------
matrix.subset<-function(normalized_data, ngenes_keep = 1000, verbose=TRUE, log=FALSE){
  tsmessage("[!geneselection] Onset...", verbose=verbose)
  df<-get_variable_gene(normalized_data$m)
  gc()
  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[ngenes_keep]
  df$used<-df$dispersion_norm >= disp_cut_off
  
  features = head(order(-df$dispersion_norm),ngenes_keep)
  system("mkdir ./tmp/",ignore.stderr = T)
  system("rm ./tmp/genes",ignore.stderr = T)
  write.csv(features, file = "./tmp/genes", quote = F,row.names = F)
  write.csv(normalized_data$use_genes[features], file = "./tmp/genes_used_all", quote = F,row.names = F)
  tsmessage("[!normalization] Complete.", verbose=verbose)
  
  print("[!log-transform] Onset...")
  genes = read.csv(file = "./tmp/genes")
  features = genes$x
  
  # Log transformation
  m_n_whole <- normalized_data$m[,features]

  m_filt <- Matrix(m_n_whole, sparse = T)
  if(log){
    m_filt <- Matrix(log2(m_n_whole + 1), sparse = T)    
  }
  
  m_n_whole <- scaleMatrix(m_filt, rescale = F, centre = T, verbose = T)
  
  return(m_n_whole)
}

get_variable_gene<-function(m) {
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion), bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}


#' Gene Frequency - Inverse Cell Frequency (GF-ICF)
#'
#' R implementation of the GF-ICF (https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract)
#' Thanks to 3’-end scRNA-seq approaches, we can now have an accurate estimation of gene expression without having to account for gene length,
#' thus the number of transcripts (i.e. UMI) associated to each gene, strictly reflects the frequency of a gene in a cell, exactly like a word in a document.
#' GFICF (Gene Frequency - Inverce Cell Frequency) is analugous of TF-IDF scoring method as defined for tex dada. With GFICF we consider a cell to be analogous
#' to a document, genes analogous to words and gene counts to be analogous of the word’s occurrence in a document.
#' 
#' @param M Matrix; UMI cell count matrix
#' @param cell_proportion_max integer; Remove genes present in more then to the specifided proportion (0,1). Default 1.
#' @param cell_proportion_min integer; Remove genes present in less then or equal to the specifided proportion (0,1). Default is 0.05 (i.e. 5 percent).
#' @param storeRaw logical; Store UMI counts.
#' @param  normalize logical; Rescale UMI counts before applay GFICF. Recaling is done using EdgeR normalization.
#' @param verbose boolean; Increase verbosity.
#' @return The updated gficf object.
#' @export
gficf = function(M,cell_proportion_max = 1,cell_proportion_min = 0.05,storeRaw=TRUE,normalize=TRUE,verbose=TRUE)
{
  data = list()
  M = normCounts(M,doc_proportion_max = cell_proportion_max,doc_proportion_min = cell_proportion_min,normalizeCounts=normalize,verbose=verbose)
  data$gficf = tf(M,verbose = verbose)
  if (storeRaw) {data$rawCounts=M;rm(M)}
  data$w = getIdfW(data$gficf,verbose = verbose)
  data$gficf = idf(data$gficf,data$w,verbose = verbose)
  data$gficf = t(l.norm(t(data$gficf),norm = "l2",verbose = verbose))
  gc()
  
  data$param <- list()
  data$param$cell_proportion_max = cell_proportion_max
  data$param$cell_proportion_min = cell_proportion_min
  data$param$normalized = normalize
  return(data)
}

#' @import Matrix
#' @importFrom edgeR DGEList calcNormFactors cpm
#' 
normCounts = function(M,doc_proportion_max = 1,doc_proportion_min = 0.01,normalizeCounts=FALSE,verbose=TRUE)
{
  ix = Matrix::rowSums(M!=0)
  M = M[ix>ncol(M)*doc_proportion_min & ix<=ncol(M)*doc_proportion_max,]
  
  if (normalizeCounts) 
  {
    require(edgeR) 
    require(DGEList)
    require(calcNormFactors)
    require(cpm)
    tsmessage("[!normalization] Normalize counts..",verbose = verbose)
    M <- Matrix::Matrix(cpm(calcNormFactors(DGEList(counts=M),normalized.lib.sizes = T)),sparse = T) 
  } 
  
  return(M)
}


#' @import Matrix
#' 
tf = function(M,verbose)
{
  
  tsmessage("[!normalization] Apply GF transformation..",verbose = verbose)
  M =t(t(M) / Matrix::colSums(M))
  
  return(M)
}

#' @import Matrix
#' 
idf = function(M,w,verbose)
{
  tsmessage("[!normalization] Apply ICF..",verbose = verbose)
  M = M[rownames(M) %in% names(w),]
  if(nrow(M)<length(w))
  {
    g = names(w)[!names(w)%in%rownames(M)]
    tmp = Matrix::Matrix(data = 0,nrow = length(g),ncol = ncol(M))
    rownames(tmp) = g
    colnames(tmp) = colnames(M)
    M = rbind(M,tmp)
  }
  M = M[names(w),]
  M = M * w
  return(M)
}

#' @import Matrix
#' 
getIdfW = function(M, type="classic",verbose)
{
  tsmessage("[!normalization] Compute ICF weigth..",verbose = verbose)
  nt = Matrix::rowSums(M!=0)
  if (type == "classic") {w = log( (ncol(M)+1) / (nt+1) );rm(nt)}
  if (type == "prob") {w = log( (ncol(M) - nt) / nt );rm(nt)}
  if (type == "smooth") {w = log( 1 + ncol(M)/nt );rm(nt)}
  return(w)
}

l.norm = function (m, norm = c("l1", "l2"),verbose) 
{
  tsmessage(paste("[!normalization] Apply ",norm), verbose = verbose)
  norm_vec = switch(norm, l1 = 1/(rowSums(m)), l2 = 1/sqrt(rowSums(m^2)))
  norm_vec[is.infinite(norm_vec)] = 0
  if (inherits(m, "sparseMatrix")) 
    Diagonal(x = norm_vec) %*% m
  else m * norm_vec
}

#' @import Matrix
#' 
scaleMatrix = function(x,rescale,centre,verbose)
{
  if (rescale)
  {
    tsmessage("[!normalization] Apply Rescaling...", verbose = verbose)
    
     bc_tot <- Matrix::rowSums(x)
     median_tot <- stats::median(bc_tot)
     x <- base::sweep(x, 1, median_tot/bc_tot, '*')
    
    # library(scales)
    # x <- apply(x, 2, rescale)
  }
  
  if (centre)
  {
    tsmessage("[!normalization] Apply Centering...", verbose = verbose)
    
    x <- base::sweep(x, 2, Matrix::colMeans(x), '-')
    x <- base::sweep(x, 2, base::apply(x, 2, sd), '/')
  }
  return(x)
}


stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
tsmessage <- function(..., domain = NULL, appendLF = TRUE, verbose = TRUE,time_stamp = TRUE) {
  if (verbose) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(stime(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}