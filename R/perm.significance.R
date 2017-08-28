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
#' @param num.perms Number of permutations used to assess significance (default = 1e3).
#'
#' @param random.seed Random seed (default = NULL).
#'
#' @param alternative A character string ("greater" or "less") that specifies the direction of the alternative hypothesis, 
#'	either rho > 0 or rho < 0 (default = "greater").
#'
#' @return Returns a five-column matrix.  The first three columns are the same as gene.annot.  The fourth column contains
#'	gene-specific Pearson or Spearman correlation coefficients based on the entries in each row of exp.mat and cn.mat, 
#'	respectively (column name = "R").  The fifth column contains squared Pearson correlation coefficients (column name = "R^2").  
#'	The sixth column contains the permutation-based right-tailed p-value of the correlation coefficient (column name = "perm_pValue").
#'	The seventh column contains Benjamini-Hochberg q-values corresponding to the p-values.  Genes with constant gene expression
#'	or DNA copy number are removed because they have zero variance.
#'
#' @examples exp.mat = tcga.exp.convert(exp.mat)
#'
#'  cn.mat = tcga.cn.convert(cn.mat)
#'
#'  genes = c("MYEOV", "CCND1", "ORAOV1", "FGF19", "FGF4", "FGF3", "ANO1", "PPFIA1")
#'
#'  pd = data.prep(exp.mat, cn.mat, gene.annot, sample.annot, log.exp = FALSE, gene.list = genes)
#'
#'  pd.exp = pd[["exp"]]
#'
#'  pd.cn = pd[["cn"]]
#'
#'  pd.ga = pd[["gene.annot"]]
#'
#'  perm.significance(pd.exp, pd.cn, pd.ga)
#'
#' @export
perm.significance = function(
	exp.mat, 
	cn.mat, 
	gene.annot,
	method = "pearson",
	digits = 5,
	num.perms = 1e3,
	random.seed = NULL,
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
	cn.mat2 = cn.mat - rowMeans(cn.mat)
	cn.mat2 = cn.mat2/sqrt(rowSums(cn.mat2^2))
	exp.mat2 = exp.mat - rowMeans(exp.mat)
	exp.mat2 = exp.mat2/sqrt(rowSums(exp.mat2^2))
	r.vec = rowSums(cn.mat2 * exp.mat2)

	#Set random seed, if appropriate
	if (!is.null(random.seed))
		{
		set.seed(random.seed)
		}

	#Recompute correlation coefficients after permutation
	perm.results = matrix(NA, nrow(exp.mat), num.perms)
	for (i in (1:num.perms))
		{
		temp.perm = sample(c(1:ncol(exp.mat)))
		perm.cn.mat = cn.mat[,temp.perm]
		perm.cn.mat = perm.cn.mat - rowMeans(perm.cn.mat)
		perm.cn.mat = perm.cn.mat/sqrt(rowSums(perm.cn.mat^2))
		perm.exp.mat = exp.mat - rowMeans(exp.mat)
		perm.exp.mat = perm.exp.mat/sqrt(rowSums(perm.exp.mat^2))
		perm.results[,i] = rowSums(perm.cn.mat * perm.exp.mat)
		}

	#Assess the significance of the entries in r.vec using rows of perm.results
	perm.results = cbind(r.vec, perm.results)
	perm.ranks = t(apply(perm.results, 1, rank))
	perm.pvals = (as.numeric(alternative == "greater") * (num.perms + 2 - perm.ranks[,1])/num.perms) +
		(as.numeric(alternative == "less") * (perm.ranks[,1])/num.perms)
	perm.pvals[which(perm.pvals > 1)] = 1
	perm.qvals = p.adjust(perm.pvals, method = "BH")

	#Create and return output
	output = cbind(gene.annot, signif(r.vec, digits), signif(r.vec^2, digits),
		signif(perm.pvals, digits), signif(perm.qvals, digits))
	colnames(output)[4] = "R"
	colnames(output)[5] = "R^2"
	colnames(output)[6] = "perm_pValue"
	colnames(output)[7] = "perm_qValue"
	rownames(output) = rownames(exp.mat)
	return(output)
	}

