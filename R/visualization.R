

#' \code{plotBarcode} 
#' @title plotBarcode
#' @export
#' @description plots the intensity of a frequency matrix with a 4-fold color step
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param data : a frequency matrix to be visualized
#' @param main : the main title of the plot empty by default
#' @return nothing
#' @note this may be slightly affected by the size of the catalogue when comparing different studies
plotBarcode <- function(data, main=""){
  cols <- list()
  cols$val <- c(0, 0, 0.0000001, 0.0000004, 0.0000016, 0.0000064, 0.0000256, 0.0001024, 0.0004096, 0.0016384)
  cols$col <- c("white", "deepskyblue", "blue", "green3", "yellow", "orange", "red", "orangered2", "darkred")
  #cols$col <- c("white", "skyblue", "deepskyblue3", "green3", "yellow", "orange", "red", "darkred", "black")
  image(t(data[nrow(data):1,]), breaks=cols$val, col=cols$col, axes=FALSE, main=main)
  box()
}

#' \code{plotBarcode2} 
#' @title plotBarcode2
#' @export
#' @description plots the intensity of a frequency matrix with a 4-fold color step. Usually used
#'      for complex figures where different MGS are overlapped and annotated with different data.
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param data : a frequency matrix to be visualized
#' @param main : the main title of the plot empty by default
#' @param ylabl : label for the left y axis
#' @param ylabr : label for the right y axis
#' @param col.axisl : color for the left y axis
#' @param col.axisr : color for the right y axis
#' @param box : default FALSE.
#' @return nothing
#' @note this may be slightly affected by the size of the catalogue when comparing different studies
plotBarcode2 <- function(data, main="",ylabl="", ylabr="", col.axisl="white", col.axisr="white", box=FALSE){
  #par(mai = par()$mai-c(0,1,0,0))
  cols <- list()
  cols$val <- c(0, 0, 0.0000001, 0.0000004, 0.0000016, 0.0000064, 0.0000256, 0.0001024, 0.0004096, 0.0016384)
  cols$col <- c("white", "deepskyblue", "blue", "green3", "yellow", "orange", "red", "orangered2", "darkred")
  image(t(data[nrow(data):1,]), breaks=cols$val, col=cols$col, axes=FALSE, main=main)
  axis(2,at=0.5, labels=ylabl, las=1, cex.axis=2, col.axis = col.axisl, col.ticks="white")
  axis(4,at=0.5, labels=ylabr, las=1, cex.axis=2, col.axis = col.axisr, col.ticks="white")
  if(box){
    #par(mai = par()$mai+c(0,1,0,0))
    box(col="darkgray",cex=0.5, which="plot")
  }
}

#' \code{plotBarcodeBW} 
#' @title plotBarcodeBW
#' @export
#' @description plots in black when a signal is different from zero and white otherwise
#' @author Edi Prifti
#' @param data : a frequency matrix to be visualized
#' @param main : the main title of the plot empty by default
#' @return nothing
plotBarcodeBW <- function(data, main=""){
  cols <- list()
  cols$val <- c(0, 0, 1)
  cols$col <- c("white", "black")
  image(t(data[nrow(data):1,]), breaks=cols$val, col=cols$col, axes=F, main=main)
  box()
}

#' \code{plotPvals} 
#' @title plotPvals
#' @export
#' @description plots a heatmap of a matrix composed of p-values. Gray is not significant at p=0.05
#'      and significance decrases from skyblue to darkred
#' @author Edi Prifti
#' @param data : a frequency matrix to be visualized
#' @param main : the main title of the plot empty by default
#' @return nothing
plotPvals <- function(data, main=""){
  cols <- list()
  cols$val <- c(0, 0, 1e-20, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 5e-2,1)
  cols$col <- c("white", "skyblue", "deepskyblue3", "green3", "yellow", "orange", "red", "darkred", "black","gray")
  image(t(data[nrow(data):1,]), breaks=cols$val, col=cols$col, axes=F, main=main)
  box()
}

#' \code{plotCors} 
#' @title plotCors
#' @export
#' @description plots a heatmap of a matrix composed of correlation values from -1 to 1. The blue 
#'      colors are negative correlations while the red are positive
#' @author Edi Prifti
#' @param data : a frequency matrix to be visualized
#' @param main : the main title of the plot empty by default
#' @return nothing
plotCors <- function(data, main=""){
  cols <- list()
  cols$val <- c(-1, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1)
  cols$col <- c("darkblue", "blue3", "blue2", "lightblue", "white", "white", "orange","red3" ,"darkred", "black")
  image(t(data[nrow(data):1,]), breaks=cols$val, col=cols$col, axes=F, main=main)
  box()
}


#' \code{plotCors2} 
#' @title plotCors2
#' @export
#' @description similar to plotCors but with more levels of colors
#' @author Edi Prifti
#' @param data : a frequency matrix to be visualized
#' @param main : the main title of the plot empty by default
#' @param grid : add the grid on top of the image (default=TRUE)
#' @param axes : plot the axes (default=TRUE)
#' @param topdown : reverse order of rows in the matrix (default=TRUE)
#' @return nothing
plotCors2 <- function (data, main = "", grid=TRUE, axes=TRUE, topdown=TRUE) {
  col <- list()
  col$val <- seq(-1, 1, 0.1)
  col$col <- gplots::colorpanel(20, low = "blue", mid = "white", high = "red")
  if(topdown) image(t(data[nrow(data):1, ]), breaks = col$val, col = col$col, axes = F, main = main)
  else image(t(data), breaks = col$val, col = col$col, axes = F, main = main)
  box()
  if(grid) grid(nx = ncol(data), ny = nrow(data), col = "black")
  if(axes){
    axis(1, at = seq(0,1,length.out=ncol(data)), labels = colnames(data), las = 2)
    if(topdown) axis(2, at = seq(0,1,length.out=nrow(data)), labels = rev(rownames(data)), las = 2)
    else axis(2, at = seq(0,1,length.out=nrow(data)), labels = rownames(data), las = 2)
  }
}

#' \code{plotCol} 
#' @export
#' @title plotCol: a matrix of colored cells when debugging
#' @import grid
#' @description a matrix of colored cells when debugging
#' @author Edi Prifti
#' @param col : a vector of colors
#' @param nrow : number of rows in of the matrix (1 = default)
#' @param ncol : number of columns in of the matrix. If nothing is used it will be computed automatically to fit the cols.
#' @param txt.col : font color.
#' @return nothing
plotCol <- function(col, nrow=1, ncol=ceiling(length(col) / nrow), txt.col="black") {
  stopifnot(nrow >= 1, ncol >= 1)
  if(length(col) > nrow*ncol)
    warning("some colors will not be shown")
  grid.newpage()
  gl <- grid.layout(nrow, ncol)
  pushViewport(viewport(layout=gl))
  ic <- 1
  for(i in 1:nrow) {
    for(j in 1:ncol) {
      pushViewport(viewport(layout.pos.row=i, layout.pos.col=j))
      grid.rect(gp= gpar(fill=col[ic]))
      grid.text(col[ic], gp=gpar(col=txt.col))
      upViewport()
      ic <- ic+1
    }
  }
  upViewport()
  invisible(gl)
}

#===========================================================================================
# This section computes and visualizes scores of data variation within a given MGS cluster.
# Some of the scores were inspired from # http://en.wikipedia.org/wiki/Coefficient_of_variation
#===========================================================================================

#' \code{presenceScore} 
#' @title presenceScore
#' @export
#' @description Computes the percentage [0,1]of values of a vector that are aboove a given threshold
#' @author Edi Prifti
#' @param vect : a numerical vector
#' @param th : the threshold to be applied, default is 0
#' @return percentage
presenceScore <- function(vect, th=0){
  if(!is.numeric(vect) & !is.integer(vect)) warning("Please provide a numerical vector.")
  return(sum((vect>th)+0.0)/length(vect))
}

#' \code{abondanceScore} 
#' @title abondanceScore
#' @export
#' @description Computes the sum of the vectors divided by the prevalence
#' @author Edi Prifti
#' @param vect : a numerical vector
#' @param th : the threshold to be applied, default is 0
#' @return a mean abundance 
abondanceScore <- function (vect, th = 0) 
{
  return(sum(vect > th)/sum((vect > th) + 0)) # error corrected elc
}


#' \code{computeSignalMetrics} 
#' @title computeSignalMetrics
#' @export
#' @description Computes scores of data variation within a given MGS cluster.
#' @author Edi Prifti
#' @param dat : a matrix where operations will be performed on the columns. Please transpose if operations
#' are needed in the rows.
#' @return a matrix of results where scores are in the columns 
computeSignalMetrics <- function(dat){
  # elements for the scores
  moy<- apply(dat, 2, mean)
  std <- apply(dat, 2, sd)
  var <- apply(dat, 2, var)
  q1 <- apply(dat, 2, quantile, 0.25)
  q2 <- apply(dat, 2, quantile, 0.5)
  q3 <- apply(dat, 2, quantile, 0.75)
  
  # Variance to mean ratio
  variance_to_mean <- (var/moy)
  variance_to_mean[is.nan(variance_to_mean)] <- 0
  
  # signal to noise ratio
  signal_to_noise <- (moy/std)
  signal_to_noise[is.nan(signal_to_noise)] <- 0
  
  # Coefficient of variation
  variation_coefficient <- (std/moy)
  variation_coefficient[is.nan(variation_coefficient)] <- 0
  
  # Effiency
  efficiency <- (std^2/moy^2)
  efficiency[is.nan(efficiency)] <- 0
  
  # Quartile coefficient of dispersion
  quartile_dispertion <- (q3-q1)/q2
  quartile_dispertion[is.nan(quartile_dispertion) | is.infinite(quartile_dispertion)] <- 0
  
  res <- as.data.frame(cbind(variance_to_mean, signal_to_noise, variation_coefficient, efficiency, quartile_dispertion))
  return(res)
}


#' \code{plotMGSQuality} 
#' @title plotMGSQuality
#' @export
#' @importFrom gplots heatmap.2
#' @description Visualized scores of data variation within a given MGS cluster as well as the barcode of the MGS.
#' A subset of 50 most connected genes is also plotted the same way.
#' @author Edi Prifti
#' @param dat : a matrix where operations will be performed on the columns. Please transpose if operations
#' are needed in the rows. This is typically an MGS frequency matrix with genes in the rows.
#' @param main : the name of the plot
#' @param return.scores : weather the function should return the computed scores or not (default=TRUE).
#' @return a matrix of results where scores are in the columns 
plotMGSQuality <- function(dat, main = "mgs", return.scores = TRUE) 
{
  par(mfcol = c(6, 2), xaxs = "i", yaxs = "i", mar = c(2, 2, 2, 1))
  colpan <- c("black", colorpanel(n = 20, low = "darkblue", mid = "darkorchid", high = "red"))
  size <- 50
  if (nrow(dat) < size) {
    size <- nrow(dat)
    warning("need more than 50 genes")
  }
  dat50 <- dat[1:size, ]
  plotBarcode(dat50, main = paste(main, nrow(dat50), "genes"))
  scores50 <- computeSignalMetrics(dat50)
  for (i in 1:ncol(scores50)) {
    if (all(scores50[, i] == 0)) {
      cols = "black"
    }
    else {
      cols <- colpan[as.numeric(cut(scores50[, i], breaks = 20))]
    }
    plot(scores50[, i], pch = 19, cex = 0.4, main = colnames(scores50)[i],
         ylab = "", xlab = "individual index", col = cols)
  }
  plotBarcode(dat, main = paste(main, nrow(dat), "genes"))
  scores <- computeSignalMetrics(dat)
  for (i in 1:ncol(scores)) {
    if (all(scores[, i] == 0)) {
      cols = "black"
    }
    else {
      cols <- colpan[as.numeric(cut(scores[, i], breaks = 20))]
    }
    plot(scores[, i], pch = 19, cex = 0.4, main = colnames(scores)[i], 
         ylab = "", xlab = "individual index", col = cols)
  }
  if (return.scores) 
    return(scores)
}
