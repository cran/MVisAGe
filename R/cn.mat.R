#' DNA copy number data from 98 head and neck squamous cell carcinoma (HNSC) patients
#'
#' Quantitative gene-level DNA copy number measurements for 98 samples from The
#' Cancer Genome Atlas (TCGA) HNSC cohort.  The all_data_by_genes.txt dataset from
#' the GISTIC2 output was restricted to the first 100 columns and genes that lie
#' on chromosomes 11 and 12.  Genes appear in rows; samples appear in columns (other
#' than the first two columns described below).  Gene symbols are used as row
#' names and sample identifiers are used as column names (other than the first two columns).
#'
#' @format A matrix with 2719 rows and 100 columns
#' \describe{
#'	\item{Locus.ID}{gene identifier}
#'	\item{Cytoband}{cytoband containing the gene of interest}
#'	\item{Remaining Columns}{quantitative DNA copy number}
#'	\item{Column names}{sample identifiers (other than the first two columns)}
#'	\item{Row names}{gene symbols}
#'	}
#' @source \url{https://gdac.broadinstitute.org/}
"cn.mat"