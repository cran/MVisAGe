#' @title A Function for Computing a Vector of Pearson Correlation Coefficients
#'
#' @description This function computes Pearson correlation coefficients on a row-by-row basis for two numerical input matrices of the same size.
#'
#' @param exp.mat A matrix of gene-level expression data (rows = genes, columns = samples).  Missing values are not permitted.
#'
#' @param cn.mat A matrix of gene-level DNA copy number data (rows = genes, columns = samples).  Both genes and samples should
#'	appear in the same order as exp.mat.  Missing values are not permitted.
#'
#' @param gene.annot A three-column matrix containing gene position information.  Column 1 = chromosome number written in 
#'	the form 'chr1' (note that chrX and chrY should be written chr23 and chr24), Column 2 = position (in base pairs), Column 3 = cytoband.
#'	Genes should appear in the same order as exp.mat and cn.mat.
#'
#' @return Returns a five-column matrix.  The first three columns are the same as gene.annot.  The fourth column contains
#'	gene-specific Pearson correlation coefficients based on the entries in each row of exp.mat and cn.mat, respectively (column name = "R").
#' 	The fifth column contains squared Pearson correlation coefficients (column name = "R^2").  Genes with constant gene expression
#'	or DNA copy number are removed.
#'
#' @examples corr.results = exp.mat = tcga.exp.convert(exp.mat)
#'
#'  cn.mat = tcga.cn.convert(cn.mat)
#'
#'  prepped.data = data.prep(exp.mat, cn.mat, gene.annot, sample.annot, log.exp = FALSE)
#'
#'  corr.compute(prepped.data[["exp"]], prepped.data[["cn"]], prepped.data[["gene.annot"]])
#'
#' @export
corr.compute = function(
	exp.mat, 
	cn.mat, 
	gene.annot 
	)
	{
	#Remove rows of exp.mat and cn.mat that are constant
	zero.var.exp.rows = which(apply(exp.mat, 1, var) == 0)
	zero.var.cn.rows = which(apply(cn.mat, 1, var) == 0)
	zero.var.drop.rows = unique(c(zero.var.exp.rows, zero.var.cn.rows))
	
	if (length(zero.var.drop.rows) == nrow(exp.mat))
		{
		stop("Exiting.  At least one input matrix has no variance.")
		} else if (length(zero.var.drop.rows) > 0)
			{
			exp.mat = exp.mat[-zero.var.drop.rows,]
			cn.mat = cn.mat[-zero.var.drop.rows,]
			gene.annot = gene.annot[-zero.var.drop.rows,]
			}
	
	#Compute R^2 values for all genes
	cn.mat = cn.mat - rowMeans(cn.mat)
	cn.mat = cn.mat/sqrt(rowSums(cn.mat^2))
	exp.mat = exp.mat - rowMeans(exp.mat)
	exp.mat = exp.mat/sqrt(rowSums(exp.mat^2))
	r.vec = rowSums(cn.mat * exp.mat)

	#Create and return output
	output = cbind(gene.annot, r.vec, r.vec^2)
	colnames(output)[ncol(output) - 1] = "R"
	colnames(output)[ncol(output)] = "R^2"
	rownames(output) = rownames(exp.mat)
	return(output)
	}