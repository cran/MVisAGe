#' A Function for Plotting Smoothed Pearson Correlation Coefficients Across Multiple Chromosomes
#'
#' This function plots smoothed R or R^2 values produced by corr.list.compute() across multiple chromosomes or genomewide.
#'
#' @param plot.list A list produced by corr.list.compute().
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
#' @param plot.legend Logical value specifying whether a legend should be included.  Default = FALSE.
#'
#' @param legend.loc Character value specifing the location of the legend.  Default = "topright".  See See \code{\link{legend}}.
#'
#' @param lty.vec Vector specifying line types for plotting values in different entries of plot.list.  Default = NULL.  See \code{\link{par}}.
#'
#' @param lwd.vec Vector specifying line widths for plotting values in different entries of  plot.list.  Default = NULL.  See \code{\link{par}}.
#'
#' @param loess.span A numerical value used to control the level of smoothing.  Smoothing is performed separately for each chromosome,
#'	and loess.span effectively defines the number of genes in each smoothing window.  Default = 250.
#'
#' @param expand.size A numerical value used to control smoothing at the ends of chromosomes.  Both ends of each chromosome are artificially
#'	extended by expand.size genes, smoothing is performed on the expanded chromosome, and then the smoothed values are restricted to the size
#' 	of the original chromosome.  Default = 50.
#'
#' @param rect.colors A character vector of length two that controls the background color for each alternating chromosome.  Default = c("light gray", "gray").
#'
#' @param chr.label Logical value specifying whether chromosome numbers should appear on the plot.  Default = FALSE.
#'
#' @param xaxis.label Text used to label the x-axis of the plot.  Default = "Chromosome".  See \code{\link{plot}}.
#'
#' @param yaxis.label Text used to label the y-axis of the plot.  Default = NULL.  See \code{\link{plot}}.
#'
#' @param main.label Text used to label the plot header.  Default = NULL.  See \code{\link{par}}.
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
#' @param margin.vec Numerical vector specifying margin sizes.  Default = rep(1, 4).  See \code{\link{par}}.
#'
#' @return Creates a plot of gene-level R or R^2 values produced by corr.list.compute().  Values of R
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
#'  smooth.genome.plot(plot.list = output.list, lwd.vec = c(3, 3), lty.vec = c(1, 2))
#'
#' @export
smooth.genome.plot = function(
	plot.list, 
	plot.column = "R^2", 
	annot.colors = c("black", "red", "green", "blue", "cyan"), 
	vert.pad = 0.05, 
	ylim.low = NULL,
	ylim.high = NULL,
	plot.legend = TRUE,
	legend.loc = "bottomright",
	lty.vec = NULL,
	lwd.vec = NULL,
	loess.span = 250,
	expand.size = 50,
	rect.colors = c("light gray", "gray"),
	chr.label = TRUE,
	xaxis.label = "Chromosome",
	yaxis.label = NULL,
	main.label = NULL,
	axis.cex = 1,
	label.cex = 1,
	xaxis.line = 1.5,
	yaxis.line = 2.5,
	main.line = 0,
	margin.vec = rep(1, 4)
	)
	{
	#Restrict plot.list to a common set of genes.  This may be necessary if 
	#plot.list has length greater than 1.
	common.list.genes = names(which(table(unlist(lapply(plot.list, rownames))) == length(plot.list)))
	for (i in (1:length(plot.list)))
		{
		plot.list[[i]] = plot.list[[i]][common.list.genes,]
		}
		
	#Define terms for plotting
	chr = as.numeric(plot.list[[1]][,"chr"])
	pos = as.numeric(plot.list[[1]][,"pos"])
	chr.pos.perm = order(chr, pos)
	m = length(chr.pos.perm)
	dists = pos[chr.pos.perm][2:m] - pos[chr.pos.perm][1:(m - 1)]
	dists[dists < 0] = 0
	cumdists = cumsum(dists)
	cumdists = c(0, cumdists)/1e6
	for (i in (1:length(plot.list)))
		{
		plot.list[[i]] = plot.list[[i]][chr.pos.perm,]
		}
		
	#if (expand.size > length(plot.rows))
	#	{
	#	print("The parameter 'expand.size' is larger than the number of genes")
	#	print("in the region of interest.  Adjusting 'expand.size' accordingly.")
	#	expand.size = length(plot.rows) - 1
	#	}
		
	#Loess smoothing for plotting purposes
	x.list = vector("list", length(plot.list))
	smoothed.list = vector("list", length(plot.list))
	chr.list = vector("list", length(plot.list))
	names(x.list) = names(plot.list)
	names(smoothed.list) = names(plot.list)
	names(chr.list) = names(plot.list)

	for (i in c(1:length(smoothed.list)))
		{
		temp.xvals = c()
		temp.yvals = c()
		temp.chr = c()
		for (j in sort(unique(as.numeric(plot.list[[i]][,"chr"]))))
			{
			temp.entries = which(as.numeric(plot.list[[i]][,"chr"]) == j)
			temp.matrix = plot.list[[i]][temp.entries,]
			temp.start = min(cumdists[temp.entries])
			loess.xvals = cumdists[temp.entries]
			loess.start = min(loess.xvals)
			loess.xvals = loess.xvals - loess.start
			loess.n = length(loess.xvals)
			expand.loess.xvals = c(-loess.xvals[(expand.size + 1):2],
				loess.xvals, ((2 * max(loess.xvals)) - loess.xvals[(loess.n - 1):(loess.n - expand.size - 1)]))
			loess.yvals = as.numeric(temp.matrix[,plot.column])
			expand.loess.yvals = c(loess.yvals[(expand.size + 1):2],
				loess.yvals, loess.yvals[(loess.n - 1):(loess.n - expand.size - 1)])
			
			#temp.loess1 = loess(loess.yvals ~ loess.xvals,
			#	span = loess.span/length(loess.xvals), degree = 1, family = "gaussian")
			temp.loess2 = loess(expand.loess.yvals ~ expand.loess.xvals,
				span = loess.span/length(expand.loess.xvals), degree = 1, family = "gaussian")
			
			temp.xvals = c(temp.xvals, (temp.loess2$x[(expand.size + 1):(expand.size + loess.n)] + loess.start))
			temp.yvals = c(temp.yvals, temp.loess2$fitted[(expand.size + 1):(expand.size + loess.n)])
			temp.chr = c(temp.chr, rep(j, length(temp.entries)))
			
			j = j + 1
			}
		x.list[[i]] = temp.xvals
		smoothed.list[[i]] = temp.yvals
		chr.list[[i]] = temp.chr
		}
	x.vals = sort(unlist(x.list))
	chr.vals = sort(unlist(chr.list))
			
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

	par(mar = margin.vec)
	plot(x.vals, 
		rep(0, length(x.vals)), 
		xlim = range(cumdists), 
		ylim = c(ylim.low, ylim.high), 
		axes = FALSE, 
		ylab = "", 
		xlab = "", 
		type = "n")

	for (i in unique(chr.vals))
		{
 		rect(min(x.vals[which(chr.vals == i)]),
 			ylim.low,
			max(x.vals[which(chr.vals == i)]),
   			ylim.high,
			col = rect.colors[1 + i%%2], 
			border = rect.colors[1 + i%%2])
					
		if (chr.label == T)
			{
			if ((i %% 2) == 1)
				{
				segments(min(x.vals[which(chr.vals == i)]),
 					ylim.low[1],
					max(x.vals[which(chr.vals == i)]),
   					ylim.low[1],
					lwd = 2)

				mtext(i,
					side = 1,
					line = -.5,
	 				at = .5 *(min(x.vals[which(chr.vals == i)]) +
	      				max(x.vals[which(chr.vals == i)])), 
					adj = .5)
				}
			}
		}

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

	axis(2, cex.axis = axis.cex)
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
