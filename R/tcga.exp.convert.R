#' @title A Function for Reformatting TCGA mRNA Expression Matrices
#'
#' @description This function reformats mRNA expression matrices obtained from the Broad Institute's Firehose GDAC 
#'  (\url{https://gdac.broadinstitute.org/}) so they can be used as input for mVisAGe functions.
#'
#' @param exp.mat A matrix of mRNA expression data.  Typically illuminahiseq_rnaseqv2-RSEM_genes_normalized, or a subset thereof,
#'  including the header rows.
#'
#' @return A matrix of mRNA expression data (rows = genes, columns = samples) that is suitable for input to mVisAGe functions.
#'
#' @examples exp.mat = tcga.exp.convert(exp.mat)
#'
#' @export
tcga.exp.convert = function(
	exp.mat
	)
	{
	#Remove the top row
	exp.mat = exp.mat[c(2:nrow(exp.mat)),]
	num.mat = matrix(as.numeric(as.matrix(exp.mat)), nrow(exp.mat), ncol(exp.mat))
	rownames(num.mat) = rownames(exp.mat)
	colnames(num.mat) = colnames(exp.mat)
	exp.mat = num.mat
	rm(num.mat)
	
	#Shorten the names of the barcodes in exp.mat
	sep.char = substring(colnames(exp.mat)[1], 5, 5)
	barcode.mat = matrix(unlist(strsplit(colnames(exp.mat), sep.char, fixed = T)),
		nrow = ncol(exp.mat), ncol = 7, byrow = T)
	colnames(exp.mat) = apply(barcode.mat[,1:4], 1, paste, collapse = sep.char)
	
	#Rewrite the gene names using gene symbols, if necessary, and remove
	#repeated gene names
	gene.symbols = rownames(exp.mat)
	gene.symbols = unlist(strsplit(gene.symbols, split = "|", fixed = T))[seq(1, (2 * nrow(exp.mat) - 1), by = 2)]
	rownames(exp.mat) = gene.symbols
	repeat.gene.rows = which(rownames(exp.mat) %in% names(which(table(gene.symbols) > 1)))
	if (length(repeat.gene.rows) > 0)
		{
		exp.mat = exp.mat[-repeat.gene.rows,]
		}
		
	#Return the output
	return(exp.mat)
	}
	