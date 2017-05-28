#' Gene expression data from 100 head and neck squamous cell carcinoma (HNSC) patients
#'
#' RSEM gene expression measurements for 100 samples from The Cancer Genome Atlas (TCGA)
#' HNSC cohort after restricting to genes that lie in chromosomes 11 and 12.  Genes
#' appear in rows; samples appear in columns (other than the first two columns described 
#' below).  Gene symbols are used as row names and sample identifiers are used as 
#' column names (other than the first two columns).
#' @format A matrix with 2161 rows and 100 columns
#' \describe{
#'  \item{Columns}{RSEM gene expression measurements}
#'	\item{Column names}{sample identifiers}
#'	\item{Row names}{gene symbols}
#'	}
#' @source \url{https://gdc.cancer.gov/}
"exp.mat"