#' @title A Function for Reformatting TCGA DNA Copy Number Matrices
#'
#' @description This function reformats DNA copy number matrices obtained from the Broad Institute's Firehose GDAC 
#'  (\url{https://gdac.broadinstitute.org/}) so they can be used as input for mVisAGe functions.
#'
#' @param cn.mat A matrix of DNA copy number data included in the GISTIC2 output.  Typically all_data_by_genes.txt, or a subset thereof,
#'  including the Locus.ID and Cytoband columns.
#'
#' @return A matrix of DNA copy number data (rows = genes, columns = samples) that is suitable for input to 
#'  mVisAGe functions.
#'
#' @examples cn.mat = tcga.cn.convert(cn.mat)
#'
#' @export
tcga.cn.convert = function(
	cn.mat
	)
	{
	#Remove the first two columns of cn.mat and convert to a data matrix
	cn.mat = cn.mat[,c(3:ncol(cn.mat))]
	cn.mat = as.matrix(cn.mat)
	num.mat = matrix(as.numeric(cn.mat), nrow(cn.mat), ncol(cn.mat))
	rownames(num.mat) = rownames(cn.mat)
	colnames(num.mat) = colnames(cn.mat)
	cn.mat = num.mat
	rm(num.mat)
	
	#Shorten the names of the barcodes in cn.mat
	sep.char = substring(colnames(cn.mat)[1], 5, 5)
	barcode.mat = matrix(unlist(strsplit(colnames(cn.mat), sep.char, fixed = T)),
		nrow = ncol(cn.mat), ncol = 7, byrow = T)
	colnames(cn.mat) = apply(barcode.mat[,1:4], 1, paste, collapse = sep.char)
	
	#Return the output
	return(cn.mat)
	}