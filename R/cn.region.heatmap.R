#' @title A Function for Creating a Heatmap of DNA Copy Number Data
#'
#' @description This function creates a heatmap of DNA copy number data for a given chromosomal region.
#'
#' @param cn.mat A matrix of gene-level DNA copy number data (rows = genes, columns = samples).  
#'	DNA methylation data can also be used.  Both row names (gene names) and column names (Sample IDs) must be given.
#'
#' @param gene.annot A three-column matrix containing gene position information.  Column 1 = chromosome number written in 
#'	the form 'chr1' (note that chrX and chrY should be written chr23 and chr24), Column 2 = position (in base pairs), Column 3 = cytoband.
#'
#' @param plot.chr The chromosome used to define the region of interest.
#'
#' @param plot.start The genomic position (in base pairs) where the region starts.
#'
#' @param plot.stop The genomic position (in base pairs) where the region stops.
#'
#' @param sample.annot An optional two-column matrix of sample annotation data.  Column 1 = sample IDs, Column 2 = sample annotation 
#' (e.g. tumor vs. normal).  If NULL, sample annot will be created using the common sample IDs and a single group ('1').  Default = NULL.
#'
#' @param sample.cluster Logical values indicating whether the samples should be clustered.  Default = FALSE.
#'
#' @param low.thresh Lower threshold for DNA copy number measurements.  All values less than low.thresh are set equal to low.thresh.  Default = -2.
#'
#' @param high.thresh Upper threshold for DNA copy number measurements.  All values greater than high.thresh are set equal to high.thresh.  Default = 2.
#'
#' @param num.cols Number of distinct colors in the heatmap.  Default = 50.
#'
#' @param collist Color scheme for displaying copy number values.  Default = ("blue", "white", "red").
#'
#' @param annot.colors Character vector used to define the color scheme for sample annotation.  Default = c("black", "red", "green", "blue", "cyan").
#'
#' @param plot.sample.annot Logical value used to specify whether the sample annotation information should be plotted.  Default = FALSE.
#'
#' @param plot.list A list produced by corr.list.compute().
#'
#' @param cytoband.colors Character vector of length two used to define the color scheme for annotating the cytoband.  Default = c("gray90", "gray60").
#'
#' @return NULL
#'
#' @examples exp.mat = tcga.exp.convert(exp.mat)
#'
#'  cn.mat = tcga.cn.convert(cn.mat)
#'
#'  prepped.data = data.prep(exp.mat, cn.mat, gene.annot, sample.annot, log.exp = FALSE)
#'
#'  pd.exp = prepped.data[["exp"]]
#'
#'  pd.cn = prepped.data[["cn"]]
#'
#'  pd.ga = prepped.data[["gene.annot"]]
#'
#'  pd.sa = prepped.data[["sample.annot"]]
#'
#'  output.list = corr.list.compute(pd.exp, pd.cn, pd.ga, pd.sa)
#'
#'  cn.region.heatmap(cn.mat = pd.cn, gene.annot = pd.ga, plot.chr = 11,
#'
#'	plot.start = 0e6, plot.stop = 135e6, sample.annot = pd.sa, plot.list = output.list)
#'
#' @export
cn.region.heatmap = function(
	cn.mat,
	gene.annot,
	plot.chr,
	plot.start,
	plot.stop,
	plot.list,
	sample.annot = NULL,
	sample.cluster = F,
	low.thresh = -2,
	high.thresh = 2,
	num.cols = 50,
	collist = c("blue", "white", "red"),
	annot.colors = c("black", "red", "green", "blue", "cyan"),
	plot.sample.annot = F,
	cytoband.colors = c("gray90", "gray60")
	)
	{
	#Restrict cn.mat and corr.list to the same set of genes
	common.list.genes = names(which(table(unlist(lapply(plot.list, rownames))) == length(plot.list)))
	common.list.genes = intersect(common.list.genes, rownames(cn.mat))
	
	if (length(common.list.genes) == 0)
		{
		stop("Exiting.  No common genes in cn.mat and plot.list.")
		} else if (length(common.list.genes) > 0)
			{
			for (i in (1:length(plot.list)))
				{
				plot.list[[i]] = plot.list[[i]][common.list.genes,]
				}
			cn.mat = cn.mat[common.list.genes,]
			gene.annot = gene.annot[common.list.genes,]
			}
	
	#Use gene.annot to order the rows of cn.mat and gene.annot by genomic position
	chr.ind = as.numeric(as.numeric(plot.list[[1]][,"chr"]) == plot.chr)
	start.ind = as.numeric(as.numeric(plot.list[[1]][,"pos"]) >= plot.start)
	stop.ind = as.numeric(as.numeric(plot.list[[1]][,"pos"]) <= plot.stop)
	plot.rows = which(chr.ind * start.ind * stop.ind == 1)
	cn.mat = cn.mat[plot.rows,]
	gene.annot = gene.annot[plot.rows,]
	for (i in (1:length(plot.list)))
		{
		plot.list[[i]] = plot.list[[i]][plot.rows,]
		}
	
	chr = as.numeric(gene.annot[,"chr"])
	pos = as.numeric(gene.annot[,"pos"])
	chr.pos.perm = order(chr, pos)
	cn.mat = cn.mat[chr.pos.perm,]
	gene.annot = gene.annot[chr.pos.perm,]
	for (i in (1:length(plot.list)))
		{
		plot.list[[i]] = plot.list[[i]][chr.pos.perm,]
		}
	
	#Define colors
	cn.mat[cn.mat < low.thresh] = low.thresh
	cn.mat[cn.mat > high.thresh] = high.thresh
	ColorRamp = colorRampPalette(collist)(num.cols)
	ColorLevels = seq(from = low.thresh, to = high.thresh, length = num.cols)
	ColorRamp_ex = ColorRamp[round(1 + (min(cn.mat) - low.thresh) * num.cols/(high.thresh - low.thresh)):round((max(x = cn.mat) - low.thresh) * num.cols/(high.thresh - low.thresh))]
	
	#Define layout
	lmat = matrix(c(0, 1, 3, 2), 2, 2, byrow = T)
	lhei = c(.25, 5)
	lwid = c(.25, 5)
	layout(lmat, lhei, lwid)
	
	#Cytoband annotation
	cytoband = sort(unique(plot.list[[1]][,"cytoband"]))
	cytoband.matrix = matrix(NA, nrow(cn.mat), length(cytoband))
	for (i in c(1:length(cytoband)))
		{
		cytoband.matrix[,i] = i * as.numeric(plot.list[[1]][,"cytoband"] == cytoband[i])
		}
	p.arm.entries = grep("p", plot.list[[1]][,"cytoband"])
	q.arm.entries = grep("q", plot.list[[1]][,"cytoband"])
	cytoband.vec = c(rep(1, length(p.arm.entries)), rep(2, length(q.arm.entries)))
	par(mar = c(.25, .25, .25, 2))
	image(as.matrix(cytoband.vec), axes = F, col = cytoband.colors, bty = "n")
	
	#Copy number heatmap
	par(mar = c(2, .25, .25, 2))
	if (!is.null(sample.annot) & sample.cluster) {
		#Cluster within each group defined by sample.annot
		clustered.samples = c()
		for (j in unique(sample.annot[,2]))
			{
			temp.subjects = sample.annot[which(sample.annot[,2] == j), 1]
			temp.cn = cn.mat[,temp.subjects] + rnorm(length(cn.mat[,temp.subjects]), 0, 1e-4)
			temp.dissim = (1 - cor(temp.cn))/2
			temp.clust = hclust(as.dist(temp.dissim), method = "average")
			clustered.samples = c(clustered.samples,
				temp.clust[["labels"]][temp.clust[["order"]]])
			}
		image(cn.mat[,clustered.samples], col=ColorRamp_ex, axes = F, bty = "n")
		} else if (sample.cluster) 
			{
			temp.dissim = (1 - cor(cn.mat))/2
			temp.clust = hclust(as.dist(temp.dissim), method = "average")
			clustered.samples = temp.clust[["labels"]][temp.clust[["order"]]]
			image(cn.mat[,clustered.samples], col = ColorRamp_ex, axes = F, bty = "n")
			} else 
				{
				image(cn.mat, col = ColorRamp_ex, axes = F, bty = "n")
				}
	
	#Sample annotation
	if (!is.null(sample.annot) & plot.sample.annot)
		{
		sample.groups = sort(unique(sample.annot[,2]))
		annot.matrix = matrix(NA, ncol(cn.mat), length(sample.groups))
		for (i in (1:length(sample.groups)))
			{
			annot.matrix[,i] = i * as.numeric(sample.annot[,2] == sample.groups[i])
			}
		annot.vec = rowSums(annot.matrix)
		names(annot.vec) = sample.annot[,1]
		#table(annot.vec, sample.annot[,2])
		if (sample.cluster)
			{
			#annot.vec = annot.vec[clustered.samples]
			par(mar = c(2, .25, .25, .25))
			image(t(as.matrix(annot.vec[clustered.samples])), axes = F, 
				col = annot.colors[c(1:length(sample.groups))], bty = "n")
			}
		}
	}
