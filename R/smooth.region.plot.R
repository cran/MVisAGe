#' A Function for Plotting Smoothed Pearson Correlation Coefficients Genomewide
#'
#' This function plots smoothed R or R^2 values produced by corr.list.compute() genomewide.
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
#' @param annot.colors A vector of colors used for plotting values in different entries of plot.list.  Default = c("black", "red", "green", "blue", "cyan").
#'
#' @param vert.pad Amount of vertical white space in the plot.  Default = 0.
#'
#' @param ylim.low Smallest value on the y-axis (used to control the range of values on the y-axis).  Default = NULL.
#'
#' @param ylim.high Largest value on the y-axis (used to control the range of values on the y-axis).  Default = NULL.
#'
#' @param lty.vec Vector specifying line types for plotting values in different entries of plot.list.  Default = NULL.  See \code{\link{par}}.
#'
#' @param lwd.vec Vector specifying line widths for plotting values in different entries of  plot.list.  Default = NULL.  See \code{\link{par}}.
#'
#' @param plot.legend Logical value specifying whether a legend should be included.  Default = FALSE.
#'
#' @param legend.loc Character value specifing the location of the legend.  Default = "topright".  See See \code{\link{legend}}.
#'
#' @param loess.span A numerical value used to control the level of smoothing.  Smoothing is performed separately for each chromosome,
#'	and loess.span effectively defines the number of genes in each smoothing window.  Default = 250.
#'
#' @param expand.size A numerical value used to control smoothing at the ends of the region of interest.  Both ends of the region are artificially
#'	extended by expand.size genes, smoothing is performed on the expanded region, and then the smoothed values are restricted to the size
#' 	of the original region.  Default = 50.
#'
#' @param xaxis.label Text used to label the x-axis of the plot.  Default = "Chromosome".  See \code{\link{plot}}.
#'
#' @param yaxis.label Text used to label the y-axis of the plot.  Default = NULL.  See \code{\link{plot}}.
#'
#' @param main.label Text used to label the plot header.  Default = NULL.  See \code{\link{plot}}.
#'
#' @param axis.cex Numerical value used to specify the font size on the axes.  Default = 1.  See \code{\link{par}}.
#'
#' @param label.cex Numerical value used to specify the font size for the axis labels.  Default = 1.  See \code{\link{par}}.
#'
#' @param xaxis.line Numerical value used to specify location of xaxis.label.  Default = 0.  See \code{\link{mtext}}.
#'
#' @param yaxis.line Numerical value used to specify location of yaxis.label.  Default = 0.  See \code{\link{mtext}}.
#'
#' @param main.line Numerical value used to specify location of main.label.  Default = 0.  See \code{\link{mtext}}.
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
#'  smooth.region.plot(plot.list = output.list, plot.chr = 11, plot.start = 0e6, plot.stop = 135e6)
#'
#' @export
smooth.region.plot = function(
	plot.list,
	plot.chr,
	plot.start,
	plot.stop,
	plot.column = "R^2", 
	annot.colors = c("black", "red", "green", "blue", "cyan"), 
	vert.pad = 0.05, 
	ylim.low = NULL,
	ylim.high = NULL,
	plot.legend = TRUE,
	legend.loc = "topleft",
	lty.vec = NULL,
	lwd.vec = NULL,
	loess.span = 50,
	expand.size = 50,
	xaxis.label = "Position (Mb)",
	yaxis.label = NULL,
	main.label = NULL,
	axis.cex = 1,
	label.cex = 1,
	xaxis.line = 1.5,
	yaxis.line = 2.5,
	main.line = 0
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
	na.list.rows = c()
	for (i in (1:length(plot.list)))
		{
		plot.list[[i]] = plot.list[[i]][plot.rows,]
		na.list.rows = c(na.list.rows, which(is.na(plot.list[[i]][,plot.column])))
		}
		
	na.list.rows = unique(na.list.rows)
	if (length(na.list.rows) > 0)
		{
		for (i in (1:length(plot.list)))
			{
			plot.list[[i]] = plot.list[[i]][-na.list.rows,]
			}
		}
		
	#Define other terms for plotting
	m = length(plot.rows)
	dists = as.numeric(plot.list[[1]][,"pos"])[c(2:m)] - 
		as.numeric(plot.list[[1]][,"pos"])[c(1:(m - 1))]
	dists = c(0, dists)
	cumdists = cumsum(dists)
	
	#Loess smoothing for plotting purposes
	x.list = vector("list", length(plot.list))
	smoothed.list = vector("list", length(plot.list))
	names(x.list) = names(plot.list)
	names(smoothed.list) = names(plot.list)
	
	for (i in c(1:length(smoothed.list)))
		{
		loess.xvals = cumdists
		loess.start = min(loess.xvals)
		loess.xvals = loess.xvals - loess.start
		loess.n = length(loess.xvals)
		expand.loess.xvals = c(-loess.xvals[(expand.size + 1):2],
			loess.xvals, ((2 * max(loess.xvals)) - loess.xvals[(loess.n - 1):(loess.n - expand.size - 1)]))
		loess.yvals = as.numeric(plot.list[[i]][,plot.column])
		expand.loess.yvals = c(loess.yvals[(expand.size + 1):2],
			loess.yvals, loess.yvals[(loess.n - 1):(loess.n - expand.size - 1)])
		
		temp.loess2 = loess(expand.loess.yvals ~ expand.loess.xvals,
			span = loess.span/length(expand.loess.xvals), degree = 1, family = "gaussian")
		
		temp.xvals = (temp.loess2$x[(expand.size + 1):(expand.size + loess.n)] + loess.start)
		temp.yvals = temp.loess2$fitted[(expand.size + 1):(expand.size + loess.n)]
	
		x.list[[i]] = (temp.xvals + as.numeric(plot.list[[i]][1, "pos"]))/1e6
		smoothed.list[[i]] = temp.yvals
		}
	x.vals = sort(unlist(x.list))
			
	#Create a vector of y-coordinates, then define the limits on the y-axis
	plot.vals = c()
	for (i in (1:length(smoothed.list)))
		{
		plot.vals = c(plot.vals, as.numeric(smoothed.list[[i]]))
		}	
	if (!is.null(ylim.low))
		{
		ylim.low = min(ylim.low, min(plot.vals)) - vert.pad
		} else ylim.low = min(plot.vals) - vert.pad
	if (!is.null(ylim.high))
		{
		ylim.high = max(ylim.high, max(plot.vals)) + vert.pad
		} else ylim.high = max(plot.vals) + vert.pad
		
	plot(x.vals, 
		rep(0, length(x.vals)), 
		xlim = range(x.vals),
		ylim = c(ylim.low, ylim.high), 
		axes = FALSE, 
		ylab = "", 
		xlab = "", 
		type = "n")

	if (is.null(lty.vec))
		{
		lty.vec = rep(1, length(smoothed.list))
		}
	if (is.null(lwd.vec))
		{
		lwd.vec = rep(1, length(smoothed.list))
		}
	
	for (i in c(1:length(smoothed.list)))
		{
		lines(x.list[[i]], smoothed.list[[i]], lwd = lwd.vec[i], 
			col = annot.colors[i], lty = lty.vec[i])
		}

	axis(side = 1, cex.axis = axis.cex)
	axis(side = 2, cex.axis = axis.cex)
	mtext(xaxis.label, side = 1, cex = label.cex, line = xaxis.line)
	mtext(yaxis.label, side = 2, cex = label.cex, line = yaxis.line)
	mtext(main.label, side = 3, cex = label.cex, line = main.line)
	if (plot.legend)
		{
		legend(legend.loc, lwd = rep(3, length(smoothed.list)), 
			col = annot.colors[c(1:length(smoothed.list))], names(smoothed.list), 
			bty = "n", inset = .05, lty = lty.vec)
		}
	}
