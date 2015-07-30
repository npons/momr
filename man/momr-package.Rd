\name{momr-package}
\alias{momr-package}
\alias{momr}
\docType{package}
\title{
  Mining Metaomics Data In R
}
\description{
  momr also known as MetaOMineR is a R package that offers many functions and modules needed for 
  the analyses of quantitative metagenomics data. It is conceived for the analyses of whole NGS
  data but can be used for 16S datasets as well or other type of omics data. Developed since the 
  very beginning of the field the package has evolved and is structured around different modules 
  such as preprocessing, analysis, visualization, etc. This package contains the different algorithms
  and routines as well as some test data objects. It is used along with other data packages that
  contain the needed information to describe a given catalogue developped in the same series.
}

\details{
\tabular{ll}{
  Package: \tab momr\cr
  Type: \tab Package\cr
  Version: \tab 1.101\cr
  Date: \tab 2015-07-30\cr
  License: \tab Artistic-2.0\cr
}
  The starting point of the analyses starts with a read count matrix that has been mapped onto a gene 
  catalog. This raw read count matrix can be preprocessed through downsizing, normalization and filtering
  steps to obtain the abundance frequencies. The samples can then be clustered in different ways to 
  check for similarity and outliers. Genes can than be statistically related to a given phenotype in 
  order to select those that are of most interest (the biomarkers). Genes of interest can be projected
  onto the MGS catalog to obtain a reduced dataset of microbial entities that is to be further annotated.
  
  Updates
  \tabular{ll}{
    2015/06/24: \tab First version for CRAN. One year long changes implemented, map reduce, normalization etc...\cr
    2014/03/24: \tab First official release, licence added, and map-reduce procedures\cr
    2014/02/03: \tab added new functions phenoPairwiseRelations, extractSignificant and lmp\cr
    2014/02/03: \tab testRelations modified to give the direction of a correlation\cr
    2014/02/03: \tab MGS sample catalog udate. This version has genes sorted based on the whole metahit cohort\cr
  }
}
\author{
  \tabular{ll}{
    Authors: \tab Edi Prifti and Emmanuelle Le Chatelier\cr
    Maintainer: \tab Edi Prifti <edi.prifti [at] gmail.com>\cr
  }
}
\references{
  Le Chatelier, Emmanuelle, Trine Nielsen, Junjie Qin, Edi Prifti, Falk Hildebrand, Gwen Falony, Mathieu Almeida, et al. "Richness of Human Gut Microbiome Correlates with Metabolic Markers." Nature 500, no. 7464: 541???546.
}
\keyword{ 
  package
  shotgun metagenomics
  MGS
  MetaGenomicSpecies
  data mining
  biomarker selection
  metagenomics
}
\examples{
%% ~~ simple examples of the most important functions ~~
  
# load the package
library(momr)

#' all the data in the package
# data(package="momr")

#' load the raw and frequency test dataset
data("hs_3.3_metahit_sample_dat_raw")
data("hs_3.3_metahit_sample_dat_freq")


#' NORMALIZATION
#' This should be performed with the whole dataset (complete catalogue). 
#' But here is an exemple with the subset of the data for illustration purposes
data(hs_3.3_metahit_genesize)
norm.data <- normFreqRPKM(dat=hs_3.3_metahit_sample_dat_raw, cat=hs_3.3_metahit_genesize)


#' CLUSTERING OF SAMPLES
hc.data <- hierClust(data=hs_3.3_metahit_sample_dat_freq[,1:5], side="col")
clust.order <- hc.data$mat.hclust$order
#' order samples followin the hierarchical clustering
ordered.samples <- colnames(hs_3.3_metahit_sample_dat_freq[,1:5])[clust.order]
#' how close are the two first samples (spearman, rho)
hc.data$mat.rho[ordered.samples[1], ordered.samples[2]]
# select the samples closely related together
close.samples <- filt.hierClust(hc.data$mat.rho, hclust.method = "ward", plot = TRUE, filt = 0.5)

#' CLUSTER GENES ON THE MGS CATALOG
#' load the curated mgs data for the hs_3.3_metahit catalog
data("mgs_hs_3.3_metahit_sup500")

#' project a list of genes onto the mgs
genebag <- rownames(hs_3.3_metahit_sample_dat_freq)
mgs <- projectOntoMGS(genebag=genebag, list.mgs=mgs_hs_3.3_metahit_sup500)

#' extract the profile of a list of genes from the whole dataset
mgs.dat <- extractProfiles(mgs, hs_3.3_metahit_sample_dat_freq, silent=FALSE)

#' plot the barcodes
par(mfrow=c(5,1), mar=c(1,0,0,0))
for(i in 1:5){ 
  plotBarcode(mgs.dat[[i]])
}

#' compute the filtered vectors
mgs.mean.vect <- computeFilteredVectors(profile=mgs.dat, type="mean")


#' TEST RELATIONS
#' for the first 1000 genes
res.test <- testRelations(data=hs_3.3_metahit_sample_dat_freq[1:500,],
                          trait=c(rep(1,150),rep(2,142)),type="wilcoxon")
head(res.test)
print(paste("There are",sum(res.test$p<0.05, na.rm=TRUE),"significant genes and",
            sum(res.test$q<0.05, na.rm=TRUE), "after adjustment for multiple testing"))

res.test.mgs <- testRelations(data=mgs.mean.vect,trait=c(rep(1,150),rep(2,142)),type="wilcoxon")


#' DOWNSIZING
#' downsize the matrix
data.downsized <- downsizeMatrix(data=hs_3.3_metahit_sample_dat_raw[,1:5],level=600,repetitions=1)
colSums(data.downsized, na.rm=TRUE)

#' downsize the genecount
data.genenb <- downsizeGC(data=hs_3.3_metahit_sample_dat_raw[,1:5], level=600, repetitions=3)
par(mfrow=c(1,1), mar=c(4,4,4,4))
plot(density(colMeans(data.genenb, na.rm=TRUE)), main="density of downsized gene richness")

#' End of test file


}
