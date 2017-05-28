#' @title A Function for Creating a List of Pearson Correlation Coefficients
#'
#' @description This function uses the corr.compute() function to compute gene-specific Pearson correlation coefficients in each
#'	group of samples defined in a sample annotation matrix.
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
#' @param sample.annot An optional two-column matrix of sample annotation data.  Column 1 = sample IDs, Column 2 = sample annotation 
#' (e.g. tumor vs. normal).  If NULL, sample annot will be created using the common sample IDs and a single group ('1').  Default = NULL.
#'
#' @return Returns a list whose length is the number of unique groups defined by sample.annot.  Each entry in the list is the
#'	output of corr.compute.
#'
#' @examples output.list = exp.mat = tcga.exp.convert(exp.mat)
#'
#'  cn.mat = tcga.cn.convert(cn.mat)
#'
#'  prepped.data = data.prep(exp.mat, cn.mat, gene.annot, sample.annot, log.exp = FALSE)
#'
#'  pd.exp = prepped.data[["exp"]]
#'
#'  pd.cn = prepped.data[["cn"]]
#'
#'  pd.ga = prepped.data[["gene.annotation"]]
#'
#'  pd.sa = prepped.data[["sample.annotation"]]
#'
#'  corr.list.compute(pd.exp, pd.cn, pd.ga, pd.sa)
#'
#' @export
corr.list.compute = function(
	exp.mat, 
	cn.mat, 
	gene.annot,
	sample.annot = NULL
	)
	{
	#Create sample.annot if no sample annotation is provided.  In this case there
	#is only a one group
	if (is.null(sample.annot))
		{
		sample.annot = matrix(NA, ncol(cn.mat), 2)
		sample.annot[,1] = colnames(cn.mat)
		sample.annot[,2] = 1
		}
		
	#Separate exp.mat and cn.mat based on the information in sample.annot, if appropriate
	sample.groups = sort(unique(sample.annot[,2]))
	exp.list = vector("list", length(sample.groups))
	cn.list = vector("list", length(sample.groups))
	names(exp.list) = sample.groups
	names(cn.list) = sample.groups
	for (i in (1:length(sample.groups)))
		{
		temp.names = sort(sample.annot[which(sample.annot[,2] == sample.groups[i]), 1])
		exp.list[[sample.groups[i]]] = exp.mat[,temp.names]
		cn.list[[sample.groups[i]]] = cn.mat[,temp.names]
		}
	
	#Compute R and R^2 values for each element of exp.list and cn.list
	corr.list = vector("list", length(exp.list))
	names(corr.list) = sample.groups
	for (i in c(1:length(corr.list)))
		{
		corr.list[[sample.groups[i]]] = corr.compute(exp.list[[sample.groups[i]]], 
			cn.list[[sample.groups[i]]], gene.annot)
		}

	return(corr.list)
	}