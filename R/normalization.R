#' \code{normFreqRPKM} 
#' @title normFreqRPKM
#' @export
#' @description converts a raw count matrix onto a frequency matrix using the RPKM normalization method.
#'        This method consists of two consecutive steps, first dividing the raw counts by the length of 
#'        the gene sequence and the second shrinking the signal per column to a sum of 1
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param dat : raw counts data matrix with gene_ids as rownames
#' @param cat : the current working catalogue where the reads are mapped and counted, (i.e. hs_3.3_metahit, hs_3.9_metahit)
#'        This can also be a vector of genelength values that correspond to the number of rows in the dat matrix and are 
#'        ordered respectively
#' @param HQ_reads : To take into account a vector of high quality reads that are not mapped.
#' @return a normalized frequency matrix
normFreqRPKM <- function (dat, cat = NULL, HQ_reads=NULL) 
{
  if (is.null(cat)) {
    stop("cat should be provided as a vector containing the gene length of the catalog with the same names as the profile matrix.")
  }
  else {
    if (length(cat) != nrow(dat)) {
      stop("cat should have the same length as the number of rows in the dat profile matrix.")
    }
    else {
      genesize <- cat[match(rownames(dat), names(cat))]
    }
  }
  if (is.matrix(dat)) 
    print("The dataset is a matrix")
  else {
    dat <- as.matrix(dat)
  }
  res <- dat
  mapped <- colSums(res)
  for (i in 1:ncol(dat)) {
    res[, i] <- dat[, i]/genesize
  }
  if (is.null(HQ_reads)){
    for (i in 1:ncol(res)) {
      res[, i] <- res[, i]/sum(res[, i])
    }
  } else {
    if (length(HQ_reads) != ncol(dat))
      stop("a vector with the number of HQ_reads used for the mapping of each sample is expected")
    mapped_genesize <- colSums(res)
    cor_factor <- mapped_genesize * HQ_reads / mapped
    for (i in 1:ncol(res)) { 
      res[, i] <- res[, i]/cor_factor[i]   
    }    
  }
  return(res)
}


#' \code{normFreqTC} 
#' @title normFreqTC
#' @export
#' @description converts a raw count matrix onto a frequency matrix using the TC (total count) normalization method.
#'        This method consists of scaling the signal by the total counts per each sample
#' @author Edi Prifti
#' @param dat : raw counts data matrix with gene_ids as rownames
#' @param HQ_reads : To take into account a vector of high quality reads that are not mapped.
#' @return a normalized frequency matrix
normFreqTC <- function (dat, HQ_reads=NULL) 
{
  if (is.matrix(dat)) 
    print("The dataset is a matrix")
  else {
    dat <- as.matrix(dat)
  }
  res <- dat
  if (is.null(HQ_reads)){
    for (i in 1:ncol(res)) {
      res[, i] <- res[, i]/sum(res[, i])
    }
  } else {
    if (length(HQ_reads) != ncol(dat))
      stop("a vector with the number of HQ_reads used for the mapping of each sample is expected")
    for (i in 1:ncol(res)) { 
      res[, i] <- res[, i]/HQ_reads[i]
    }    
  }
  return(res)
}

#' \code{rdn} 
#' @title rdn
#' @export
#' @description fixes an R base issue with rounding
#' @author Edi Prifti
#' @param x : data frame or vector
#' @return rounded integer object
rdn <- function(x) 
{
  # round the shared counts following the right non biaised way (a problem with R)
  # http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
  return(trunc(x+sign(x)*0.5))
}
  
#' \code{matrix2Integer} 
#' @title matrix2Integer
#' @export
#' @description converts numeric matrix into an integer for memory optimization
#' @author Edi Prifti
#' @param data : numerical data frame
#' @param verbose : print running information (default = TRUE)
#' @return Integer matrix
matrix2Integer <- function(data, verbose=TRUE){
  dataInt <- data
  for(i in 1:ncol(data)){
    if(i%%10==0 & verbose) print(i)
    dataInt[,i] <- as.integer(data[,i])
  }; gc()
  return(dataInt)
}


#' ADD other normalization algorithms
#' End of section and file
