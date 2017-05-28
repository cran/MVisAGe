#' Gene expression, DNA copy number, and human papillomavirus infection status for 
#' head and neck squamous cell carcinoma (HNSC)
#'
#' Gene expression data (RSEM values) and DNA copy number values (log2 ratios) from
#' The Cancer Genome Atlas (TCGA) study of HNSC were obtained from the Broad Institute's
#' Firehose GDAC (https://gdac.broadinstitute.org/),  Human papillomavirus (HPV) infection
#' status for n = 279 samples in the 2015 Nature publication are included, as are hg38
#' gene annotation data obtained from Galaxy.
#'
#' @docType data
#'
#' @usage data(mVisAGe)
#'
#' @format  Four tab-delimited text files.
#'
#' @keywords datasets
#'
#' @references TCGA (2015) Nature 517: 576–582
#' (\href{http://www.nature.com/nature/journal/v517/n7536/full/nature14129.html})
#'
#' @source \href{https://gdac.broadinstitute.org/},
#' href{http://http://www.nature.com/nature/journal/v517/n7536/full/nature14129.html#supplementary-information},
#' href{https://galaxyproject.org/}
#'
#' @examples
#' data(mVisAGe)
#' ls()