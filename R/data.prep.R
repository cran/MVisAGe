#' A Function for Preparing mRNAseq and Copy Number Data Matrices
#'
#' This function prepares mRNAseq and copy number data matrices for use in other mVisAGe functions.
#'
#' @param exp.mat A matrix of gene-level expression data (rows = genes, columns = samples).  
#' 	Both row names (gene names) and column names (sample IDs) must be given.
#'
#' @param cn.mat A matrix of gene-level DNA copy number data (rows = genes, columns = samples).  
#'	DNA methylation data can also be used.  Both row names (gene names) and column names (Sample IDs) must be given.
#'
#' @param gene.annot A three-column matrix containing gene position information.  Column 1 = chromosome number written in 
#'	the form 'chr1' (note that chrX and chrY should be written chr23 and chr24), Column 2 = position (in base pairs), Column 3 = cytoband.
#'
#' @param sample.annot An optional two-column matrix of sample annotation data.  Column 1 = sample IDs, Column 2 = categorical sample annotation 
#' 	(e.g. tumor vs. normal).  If NULL, sample annot will be created using the common sample IDs and a single group ('1').  Default = NULL.
#'
#' @param log.exp A logical value indicating whether or not the expression values have been log-transformed.  Default = FALSE.
#'
#' @param gene.list Used to restrict the output to a set of genes of interest, e.g. genes identified by GISTIC as having
#' 	recurrent copy number alterations.  Default = NULL, and in this case all genes are used.
#'
#' @return Returns a list with four components:  cn, exp, gene.annot, and sample.annot.  Each of cn, exp, and gene.annot have been restricted 
#'	to a common set of genes, and these appear in the same order.  Similarly, cn, exp, and sample.annot have been restricted to a common set
#' 	of subjects that appear in the same order.
#'
#' @examples exp.mat = tcga.exp.convert(exp.mat)
#'
#'  cn.mat = tcga.cn.convert(cn.mat)
#'
#'  data.prep(exp.mat, cn.mat, gene.annot, sample.annot, log.exp = FALSE)
#'
#' @export
data.prep = function(
	exp.mat,
	cn.mat,
	gene.annot,
	sample.annot = NULL,
	log.exp = FALSE,
	gene.list = NULL
	)
	{
	#Restrict exp.mat, cn.mat, and gene.annot to a common set of genes
	common.genes = intersect(intersect(rownames(exp.mat), rownames(cn.mat)), 
		rownames(gene.annot))
	if (!is.null(gene.list))
		{
		common.genes = intersect(common.genes, gene.list)
		}
	if (length(common.genes) == 0)
		{
		stop("Exiting because there are no common genes.  Please check the input files.")
		} else	{
				exp.mat = exp.mat[sort(common.genes),]
				cn.mat = cn.mat[sort(common.genes),]
				gene.annot = gene.annot[sort(common.genes),]
				print("Checking gene names")
				print("Gene names in expression vs. CN")
				print(all(rownames(exp.mat) == rownames(cn.mat)))
				print("Gene names in expression vs. gene annotation")
				print(all(rownames(exp.mat) == rownames(gene.annot)))
				}
				
	#Remove "chr" string from gene.annot
	gene.annot[,"chr"] = unlist(strsplit(gene.annot[,"chr"], split = "r", fixed = T))[seq(2, (2*nrow(gene.annot)), by = 2)]
	
	#Rewrite "X" and "Y" in gene.annot as 23 and 24
	gene.annot[which(gene.annot[,"chr"] == "X"), "chr"] = 23
	gene.annot[which(gene.annot[,"chr"] == "Y"), "chr"] = 24
	
	#Restrict exp.mat, cn.mat, and sample.annot to a common set of subjects
	if (!is.null(sample.annot))
		{
		common.subjects = intersect(intersect(colnames(exp.mat), colnames(cn.mat)), 
			sample.annot[,1])
		} else 	{
				common.subjects = intersect(colnames(exp.mat), colnames(cn.mat))
				sample.annot = cbind(colnames(exp.mat), rep(1, ncol(exp.mat)))
				}
		
	if (length(common.subjects) == 0)
		{
		stop("Exiting.  Check that sample names are the same in the input files.")
		} else	{
				exp.mat = exp.mat[,sort(common.subjects)]
				cn.mat = cn.mat[,sort(common.subjects)]
				print("Checking sample names")
				print("Sample names in expression vs. CN")
				print(all(colnames(exp.mat) == colnames(cn.mat)))
				
				if (!is.null(sample.annot))
					{
					sample.annot = sample.annot[which(sample.annot[,1] %in% common.subjects),]
					sample.annot = sample.annot[order(sample.annot[,1]),]
					print("Sample names in expression sample annotation")
					print(all(colnames(exp.mat) == sample.annot[,1]))
					}
				}
				
	#Remove genes with missing values in either exp.mat or cn.mat
	na.rows = union(which(is.na(rowSums(exp.mat))), which(is.na(rowSums(cn.mat))))
	if (length(na.rows) > 0)
		{
		exp.mat = exp.mat[-na.rows,]
		cn.mat = cn.mat[-na.rows,]
		gene.annot = gene.annot[-na.rows,]
		print("Missing values removed")
		print((length(which(is.na(rowSums(exp.mat)))) == 0) & 
			(length(which(is.na(rowSums(cn.mat)))) == 0))
		}
				
	#Log transform exp.mat, if necessary
	exp.mat = ((1 - log.exp) * log2(exp.mat + 1)) + (log.exp * exp.mat)
	
	#Make a list containing the output, then return
	output.list = vector("list", 4)
	names(output.list) = c("exp", "cn", "gene.annot", "sample.annot")
	output.list[["exp"]] = exp.mat
	output.list[["cn"]] = cn.mat
	output.list[["gene.annot"]] = gene.annot
	output.list[["sample.annot"]] = sample.annot
	return(output.list)
	}