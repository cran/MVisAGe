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
#' @param method A character string (either "pearson" or "spearman") specifying the method used to calculate the correlation coefficient
#'    (default = "pearson").
#'
#' @param digits Used with signif() to specify the number of significant digits (default = 5).
#'
#' @param alternative A character string ("greater" or "less") that specifies the direction of the alternative hypothesis, 
#'	either rho > 0 or rho < 0 (default = "greater").
#'
#' @return Returns a eight-column matrix.  The first three columns are the same as gene.annot.  The fourth column contains
#'	gene-specific Pearson or Spearman correlation coefficients based on the entries in each row of exp.mat and cn.mat, 
#'	respectively (column name = "R").  The fifth column contains squared Pearson correlation coefficients (column name = "R^2").  
#'	The sixth column contains t statistics corresponding to the correlation coefficients (column name = "tStat").  The
#'	seventh column contains the right-tailed p-value based on the t statistic (column name = "pValue").  The eighth column
#'	contains Benjamini-Hochberg q-values corresponding to the p-values.  Genes with constant gene expression
#'	or DNA copy number are removed because they have zero variance.
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
	gene.annot,
	method = "pearson",
	digits = 5,
	alternative = "greater"
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
	
	#Convert exp.mat and cn.mat to ranks if method = "spearman"
	if (method == "spearman")
		{
		exp.mat = t(apply(exp.mat, 1, rank))
		cn.mat = t(apply(cn.mat, 1, rank))
		}
	
	#Compute correlation coefficients for all genes
	cn.mat = cn.mat - rowMeans(cn.mat)
	cn.mat = cn.mat/sqrt(rowSums(cn.mat^2))
	exp.mat = exp.mat - rowMeans(exp.mat)
	exp.mat = exp.mat/sqrt(rowSums(exp.mat^2))
	r.vec = rowSums(cn.mat * exp.mat)

	#Compute t statitics, p-values, and q-values
	t.vec = r.vec * sqrt((ncol(exp.mat) - 2)/(1 - r.vec^2))
	p.vec = (as.numeric(alternative == "greater") * (1 - pt(t.vec, ncol(exp.mat) - 2))) +
		(as.numeric(alternative == "less") * (pt(t.vec, ncol(exp.mat) - 2)))
	q.vec = p.adjust(p.vec, method = "BH")

	#Create and return output
	output = cbind(gene.annot, signif(r.vec, digits), signif(r.vec^2, digits),
		signif(t.vec, digits), signif(p.vec, digits), signif(q.vec, digits))
	colnames(output)[4] = "R"
	colnames(output)[5] = "R^2"
	colnames(output)[6] = "tStat"
	colnames(output)[7] = "pValue"
	colnames(output)[8] = "qValue"
	rownames(output) = rownames(exp.mat)
	return(output)
	}

