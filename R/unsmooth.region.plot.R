#' A Function for Plotting Pearson Correlation Coefficients in a Given Genomic Region
#'
#' This function plots unsmoothed R or R^2 values produced by corr.list.compute() in a specified genomic region.
#'
#' @param plot.list A list produced by corr.list.compute().
#'
#' @param plot.chr The chromosome for which gene-level R or R^2 values will be plotted.
#'
#' @param plot.start The genomic position (in base pairs) where the plot will start.
#'
#' @param plot.stop The genomic position (in base pairs) where the plot will stop.
#'
#' @param plot.column "R" or "R^2" depending on whether Pearson correlation coefficients or squared Pearson
#'	correlation coefficients will be plotted.  Default = "R^2".
#'
#' @param plot.points Logical value specifying whether points should be included in the plot.  Default = TRUE.
#'
#' @param plot.lines Logical values specifying whether points should be connected by lines.  Default = FALSE.
#'
#' @param gene.names Logical value specifying whether gene names should appear on the plot.  Default = FALSE.
#'
#' @param annot.colors A vector of colors used for plotting values in different entries of plot.list.  Default = c("black", "red", "green", "blue", "cyan").
#'
#' @param vert.pad Amount of vertical white space in the plot.  Default = 0.05.
#'
#' @param num.ticks Number of ticks on the x-axis.  Default = 5.
#'
#' @param ylim.low Smallest value on the y-axis (used to control the range of values on the y-axis).  Default = NULL.
#'
#' @param ylim.high Largest value on the y-axis (used to control the range of values on the y-axis).  Default = NULL.
#'
#' @param pch.vec Vector specifying point characters for plotting values in different entries of plot.list.  Default = NULL.  See \code{\link{par}}.
#'
#' @param lty.vec Vector specifying line types for plotting values in different entries of plot.list.  Default = NULL.  See \code{\link{par}}.
#'
#' @param lwd.vec Vector specifying line widths for plotting values in different entries of  plot.list.  Default = NULL.  See \code{\link{par}}.
#'
#' @param plot.legend Logical value specifying whether a legend should be included.  Default = FALSE.
#'
#' @param legend.loc Character value specifing the location of the legend.  Default = "topright".  See See \code{\link{legend}}
#'
#' @return Creates a plot of gene-level R or R^2 values produced by corr.list.compute().
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
#'  unsmooth.region.plot(plot.list = output.list, plot.chr = 11, plot.start = 69e6, plot.stop = 70.5e6)
#'
#' @export
unsmooth.region.plot = function(
	plot.list, 
	plot.chr, 
	plot.start, 
	plot.stop,
	plot.column = "R", 
	plot.points = TRUE, 
	plot.lines = TRUE, 
	gene.names = TRUE,
	annot.colors = c("black", "red", "green", "blue", "cyan"), 
	vert.pad = .05, 
	num.ticks = 5,
	ylim.low = NULL,
	ylim.high = NULL,
	pch.vec = NULL,
	lty.vec = NULL,
	lwd.vec = NULL,
	plot.legend = TRUE,
	legend.loc = "topright"
	)
	{
	#Restrict plot.list to a common set of genes.  This may be necessary if 
	#plot.list has length greater than 1.
	common.list.genes = names(which(table(unlist(lapply(plot.list, rownames))) == length(plot.list)))
	for (i in (1:length(plot.list)))
		{
		plot.list[[i]] = plot.list[[i]][common.list.genes,]
		}
	
	#Order entries in plot.list by genomic position
	chr = as.numeric(plot.list[[1]][,"chr"])
	pos = as.numeric(plot.list[[1]][,"pos"])
	chr.pos.perm = order(chr, pos)
	for (i in (1:length(plot.list)))
		{
		plot.list[[i]] = plot.list[[i]][chr.pos.perm,]
		}

	#Restrict corr.output.list to a specific genomic window
	chr.ind = as.numeric(as.numeric(plot.list[[1]][,"chr"]) == plot.chr)
	start.ind = as.numeric(as.numeric(plot.list[[1]][,"pos"]) >= plot.start)
	stop.ind = as.numeric(as.numeric(plot.list[[1]][,"pos"]) <= plot.stop)
	plot.rows = which(chr.ind * start.ind * stop.ind == 1)
	for (i in (1:length(plot.list)))
		{
		plot.list[[i]] = plot.list[[i]][plot.rows,]
		}

	#Create a vector of y-coordinates, then define the limits on the y-axis
	plot.vals = c()
	for (i in (1:length(plot.list)))
		{
		plot.vals = c(plot.vals, as.numeric(plot.list[[i]][,plot.column]))
		}	
	if (!is.null(ylim.low))
		{
		ylim.low = min(ylim.low, min(plot.vals)) - vert.pad
		} else ylim.low = min(plot.vals) - vert.pad
	if (!is.null(ylim.high))
		{
		ylim.high = max(ylim.high, max(plot.vals)) + vert.pad
		} else ylim.high = max(plot.vals) + vert.pad

	#Create the basic plot
	if (is.null(pch.vec))
		{
		pch.vec = rep(1, length(plot.list))
		}
	plot(as.numeric(plot.list[[1]][,"pos"])/1e6, 
		as.numeric(plot.list[[1]][,plot.column]), axes = F,
		xlab = "", ylab = "", main = "", xlim = c(plot.start, plot.stop)/1e6, 
		ylim = c(ylim.low, ylim.high), pch = pch.vec[1], col = annot.colors[1])

	#Add points to the plot, if appropriate
	plot.matrix = matrix(NA, nrow(plot.list[[i]]), length(plot.list))
	plot.matrix[,1] = as.numeric(plot.list[[1]][,plot.column])
	if (plot.points && (length(plot.list) >= 2))
		{
		for (i in (2:length(plot.list)))
			{
			points(as.numeric(plot.list[[i]][,"pos"])/1e6, 
				as.numeric(plot.list[[i]][,plot.column]),
				pch = pch.vec[i], col = annot.colors[i])
			plot.matrix[,i] = as.numeric(plot.list[[i]][,plot.column])
			}
		}
		
	#Add lines to the plot, if appropriate
	if (is.null(lty.vec))
		{
		lty.vec = rep(1, length(plot.list))
		}
	if (is.null(lwd.vec))
		{
		lwd.vec = rep(1, length(plot.list))
		}
	if (plot.lines)
		{
		for (i in (1:length(plot.list)))
			{
			lines(as.numeric(plot.list[[i]][,"pos"])/1e6, 
				as.numeric(plot.list[[i]][,plot.column]),
				col = annot.colors[i], lty = lty.vec[i], 
				lwd = lwd.vec[i])
			}
		}
		
	#Add gene names to the plot, if appropriate
	if (gene.names)
		{
		text(rownames(plot.list[[1]]), x = as.numeric(plot.list[[1]][,"pos"])/1e6,
			y = apply(plot.matrix, 1, max, na.rm = T), pos = 3)
		}
		
	#Put ticks on the x-axis
	plot.ticks = signif(seq(plot.start/1e6, plot.stop/1e6, by = (plot.stop/1e6 - plot.start/1e6)/(num.ticks - 1)), 3)
	
	#Add axes
	axis(side = 1, at = plot.ticks, labels = plot.ticks, cex = 1)
	axis(side = 2, cex = 1)
	if (plot.column == "R^2")
		{
		yaxis.label = expression(rho^2)
		}	else yaxis.label = expression(rho)
	mtext(side = 1, text = paste(c("chr", plot.chr, " Position (Mb)"), collapse = ""), 
		line = 2.5)
	mtext(side = 2, text = yaxis.label, line = 2.5)
	
	#Add legend, if appropriate
	if (plot.legend)
		{
		legend(legend.loc, lwd = rep(3, length(plot.list)), 
			col = annot.colors[1:length(plot.list)],
			names(plot.list), bty = "n", inset = .05, pch = pch.vec, lty = lty.vec)
		}
	}