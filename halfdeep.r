# plotting half-deep intervals and depth along an assembly

read_scaffold_lengths <- function(lengthsFilename,scaffoldsOfInterest=NULL)
	#
	# Read a file containing scaffold names and lengths. Result is either
	# (a) reduced to scaffolds of interest, in order or (b) sorted by
	# decreasing length. A column is added, giving offsets for which the
	# scaffolds can be arranged along a number line.
	#
	# typical input:
	#	scaffold_100_arrow_ctg1 29829
	#	scaffold_101_arrow_ctg1 29808
	#	scaffold_102_arrow_ctg1 28782
	#	scaffold_103_arrow_ctg1 27584
	#	scaffold_104_arrow_ctg1 27050
	#    ...
	#
	# returns, e.g.
	#	                 name    length    offset
	#	     super_scaffold_1 215872496         0
	#	     super_scaffold_2 165051286 215872496
	#	scaffold_2_arrow_ctg1 125510139 380923782
	#	     super_scaffold_Z  84526827 506433921
	#	scaffold_5_arrow_ctg1  82438829 590960748
	#    ...
	#
	{
	scaffolds = read.table(lengthsFilename,header=F,colClasses=c("character","integer"))
	colnames(scaffolds) <- c("name","length")

	if (!is.null(scaffoldsOfInterest))
		{
		# pick ordered subset
		scaffolds = scaffolds[scaffolds$name %in% scaffoldsOfInterest,]
		scaffoldsToNumber = 1:length(scaffoldsOfInterest)
		names(scaffoldsToNumber) = scaffoldsOfInterest
		scaffolds = scaffolds[order(scaffoldsToNumber[scaffolds$name]),]
		}
	else
		{
		# sort by decreasing length
		scaffolds = scaffolds[order(-scaffolds$length),]
		}

	scaffolds[,"offset"] = rev(sum(as.numeric(scaffolds$length))-cumsum(rev(as.numeric(scaffolds$length))))

	scaffolds;
	}


linearized_scaffolds <- function(scaffolds)
	#
	# Create scaffold-to-number-line mapping; input is of the form returned by
	# read_scaffold_lengths().
	#
	{
	scaffoldToOffset <- scaffolds$offset
	names(scaffoldToOffset) <- scaffolds$name

	scaffoldToOffset
	}


read_depth <- function(depthFilename,scaffoldToOffset=NULL)
	#
	# Read a file containing coverage depth. Typically the depth has been
	# averaged over non-overlapping windows, with each window represented by
	# a single location. Genomic locations are origin-1, and intervals are
	# closed.
	#
	# A scaffold-to-offset can be provided, as would be produced by
	# linearized_scaffolds(). Two columns are added, converting each interval
	# to positions along a number line.
	#
	# typical input:
	#	scaffold_100_arrow_ctg1 1    1    16.727
	#	scaffold_100_arrow_ctg1 1001 1001 19.365
	#	scaffold_100_arrow_ctg1 2001 2001 20.701
	#	scaffold_100_arrow_ctg1 3001 3001 22.764
	#	scaffold_100_arrow_ctg1 4001 4001 20.971
	#    ...
	#
	# returns, e.g.
	#	               scaffold start  end  depth          s          e
	#	scaffold_100_arrow_ctg1     1    1 16.727 1177990031 1177990031
	#	scaffold_100_arrow_ctg1  1001 1001 19.365 1177991031 1177991031
	#	scaffold_100_arrow_ctg1  2001 2001 20.701 1177992031 1177992031
	#	scaffold_100_arrow_ctg1  3001 3001 22.764 1177993031 1177993031
	#	scaffold_100_arrow_ctg1  4001 4001 20.971 1177994031 1177994031
	#    ...
	#
	{
	depth = read.table(depthFilename,header=F,colClasses=c("character","integer","integer","numeric"))
	colnames(depth) <- c("scaffold","start","end","depth")

	if (!is.null(scaffoldToOffset))
		{
		depth = depth[depth$scaffold %in% names(scaffoldToOffset),]
		depth[,"s"] = scaffoldToOffset[depth$scaffold] + depth$start
		depth[,"e"] = scaffoldToOffset[depth$scaffold] + depth$end
		}

	depth
	}


read_halfdeep <- function(halfDeepFilename,scaffoldToOffset=NULL)
	#
	# Read a file containing a list of genomic intervals. The intervals are
	# origin-1 and closed.
	#
	# A scaffold-to-offset can be provided, as would be produced by
	# linearized_scaffolds(). Two columns are added, converting each interval
	# to positions along a number line.
	#
	# typical input:
	#	scaffold_100_arrow_ctg1 9001  10000
	#	scaffold_101_arrow_ctg1 12001 14000
	#	scaffold_103_arrow_ctg1 2001  14000
	#	scaffold_106_arrow_ctg1 1     3000
	#	scaffold_106_arrow_ctg1 19001 20000
	#    ...
	#
	# returns, e.g.
	#	               scaffold start   end          s          e
	#	scaffold_100_arrow_ctg1  9001 10000 1177999031 1178000030
	#	scaffold_101_arrow_ctg1 12001 14000 1178031860 1178033859
	#	scaffold_103_arrow_ctg1  2001 14000 1178110244 1178122243
	#	scaffold_106_arrow_ctg1     1  3000 1178189381 1178192380
	#	scaffold_106_arrow_ctg1 19001 20000 1178208381 1178209380
	#    ...
	#
	{
	halfDeep = read.table(halfDeepFilename,header=F,colClasses=c("character","integer","integer"))
	colnames(halfDeep) <- c("scaffold","start","end")
	if (!is.null(scaffoldToOffset))
		{
		halfDeep = halfDeep[halfDeep$scaffold %in% names(scaffoldToOffset),]
		halfDeep[,"s"] = scaffoldToOffset[halfDeep$scaffold] + halfDeep$start
		halfDeep[,"e"] = scaffoldToOffset[halfDeep$scaffold] + halfDeep$end
		}
	
	halfDeep
	}


read_percentiles <- function(percentilesFilename)
	#
	# Extract percentile values from a shell script.
	#
	# typical input:
	#	export percentile40=52.698
	#	export percentile50=54.818
	#	export percentile60=56.938
	#	export halfPercentile40=26.349
	#	export halfPercentile50=27.409
	#	export halfPercentile60=28.469
	#
	# returns, e.g.
	#	    percentile40     percentile50     percentile60
	#	          52.698           54.818           56.938
	#	 halfPercentile40 halfPercentile50 halfPercentile60 
	#	           26.349           27.409           28.469 
	#
	{
	percentilesCommand = paste("cat",percentilesFilename,'| sed "s/^export //" | tr "=" " "')
	percentiles = read.table(pipe(percentilesCommand),header=F,colClasses=c("character","numeric"))
	colnames(percentiles) <- c("name","value")
	percentileToValue <- percentiles$value
	names(percentileToValue) <- percentiles$name

	percentileToValue
	}


halfdeep_plot <- function(scaffolds,depth,halfDeep,percentileToValue,
                          assemblyName="",
                          scaffoldsToPlot=NULL,
                          plotFilename=NULL,
                          tickSpacing=10000000,
                          width=17,height=7,pointsize=18,
                          yLabelSpace=7,
                          maxDepth=NA,
                          scaffoldInterval=NULL)
	#
	# Plot half-deep intervals and depth along an assembly
	#
	# scaffolds:              As returned by read_scaffold_lengths().
	# depth:                  As returned by read_depth().
	# halfDeep:               As returned by read_halfdeep(). If this is NULL,
	#                         no half-deep information is displayed.
	# percentileToValue:      As returned by read_percentiles()
	# assemblyName:           Name of the assembly. This contributes to the
	#                         plain title, and can contribute to the plot
	#                         file name.
	#                         Example: "bAlcTor1.pri.cur.20190613"
	# scaffoldsToPlot=NULL:   A vector of scaffolds to restrict the plot to.
	#                         By default, everything in scaffolds[] is
	#                         plotted.
	# plotFilename=NULL:      File to plot the data to. If this is NULL, the
	#                         data is plotted on the screen.
	# tickSpacing:            Spacing of evenly-spaced ticks along the
	#                         horizontal axis. If this is zero, these ticks
	#                         are inhibited. The default is 10Mbp.
	# width,height,pointsize: These are passed through to whatever function
	#                         creates the plot window.
	# yLabelSpace:            Space below the plot. This can be increased to
	#                         accommodate longer scaffold names.
	# maxDepth:               Depth greater than this is not shown in the plot;
	#                         i.e. the vertical axis stops at this value. By
	#                         default this is 1.5*median depth.
	# scaffoldInterval        subinterval to restrict the plot to. Typically
	#                         this would only be used when only one scaffold
	#                         is to be plotted. This is a (start,end) pair,
	#                         origin-zero, half-open.
	#
	{
	# if we have a scaffold subset, reduce our copy of the data to that subset

	if (!is.null(scaffoldInterval))
		stop("scaffoldInterval is not implemented yet")

	showHalfDeep = !is.null(halfDeep)

	if (!is.null(scaffoldsToPlot))
		{
		# validate the names in the subset

		badNames = scaffoldsToPlot[!(scaffoldsToPlot %in% scaffolds$name)]
		if (length(badNames) > 0)
			stop(paste("bad scaffold name(s):",paste(scaffoldsToPlot,collapse=", ")))

		# pick ordered subset
		scaffolds = scaffolds[scaffolds$name %in% scaffoldsToPlot,]
		scaffoldsToNumber = 1:length(scaffoldsToPlot)
		names(scaffoldsToNumber) = scaffoldsToPlot
		scaffolds = scaffolds[order(scaffoldsToNumber[scaffolds$name]),]
		scaffolds[,"offset"] = rev(sum(as.numeric(scaffolds$length))-cumsum(rev(as.numeric(scaffolds$length))))
		scaffoldToOffset = linearized_scaffolds(scaffolds)

		# reduce to ordered subset
		depth = depth[depth$scaffold %in% names(scaffoldToOffset),]
		depth[,"s"] = scaffoldToOffset[depth$scaffold] + depth$start
		depth[,"e"] = scaffoldToOffset[depth$scaffold] + depth$end

		if (showHalfDeep)
			{
			halfDeep = halfDeep[halfDeep$scaffold %in% names(scaffoldToOffset),]
			halfDeep[,"s"] = scaffoldToOffset[halfDeep$scaffold] + halfDeep$start
			halfDeep[,"e"] = scaffoldToOffset[halfDeep$scaffold] + halfDeep$end
			}
		}

	# fetch percentile values (used only for drawing and labeling)

	depth50        = percentileToValue["percentile50"]
	halfDepthLo    = percentileToValue["halfPercentile40"]
	halfDepthHi    = percentileToValue["halfPercentile60"]
	depthClip      = ifelse(is.na(maxDepth),1.5*depth50,maxDepth)
	depth50Str     = sprintf("%.1f",depth50)
	halfDepthLoStr = sprintf("%.1f",halfDepthLo)
	halfDepthHiStr = sprintf("%.1f",halfDepthHi)
	depthClipStr   = sprintf("%.1f",depthClip)

	# (housekeeping)

	scaffoldTicks = 1 + c(scaffolds$offset,sum(as.numeric(scaffolds$length)))
	scaffoldCenters = (as.numeric(scaffoldTicks[1:nrow(scaffolds)])+as.numeric(scaffoldTicks[2:(nrow(scaffolds)+1)])) / 2
	if (showHalfDeep)
		halfDeepCenters = (as.numeric(halfDeep$s)+as.numeric(halfDeep$e))/2

	xlim = c(1,max(scaffoldTicks))
	ylim = c(0,depthClip)

	depthColor          = if (showHalfDeep) rgb(.6,.6,.6) else "black"
	halfDepthColor      = "black"
	halfDeepMarkerColor = "red"
	depthLimitsColor    = "blue"

	guns = col2rgb(halfDeepMarkerColor) / 255
	halfDeepOverlayColor = rgb(guns[1],guns[2],guns[3],alpha=.3)

	# open plot window or file

	turnDeviceOff = F
	if (is.null(plotFilename))
		{
		quartz(width=width,height=height)
		}
	else
		{
		print(paste("drawing to",plotFilename))
		pdf(file=plotFilename,width=width,height=height,pointsize=pointsize)
		turnDeviceOff = T
		}

	# create empty plot

	if (showHalfDeep)
		{
		title = paste("half-deep intervals (red overlay) in ",assemblyName,
		              "\nmedian=",depth50Str,
		              "   40%ile/2=",halfDepthLoStr,
		              "   60%ile/2=",halfDepthHiStr,
		              sep="")
		ylab  = paste("aligned read depth in 1Kbp windows (gray, clipped at ",depthClipStr,")",sep="")
		}
	else
		{
		title = paste("coverage depth in ",assemblyName,"\nmedian=",depth50Str,sep="")
		ylab  = paste("aligned read depth in 1Kbp windows (clipped at ",depthClipStr,")",sep="")
		}

	par(mar=c(yLabelSpace,4,2.5,0.2)+0.1)     # BLTR
	options(scipen=10)
	plot(NA,xlim=xlim,ylim=ylim,main=title,xaxt="n",xlab="",ylab=ylab)

	# add horizontal axis

	if ((showHalfDeep) && (nrow(halfDeep) > 0))
		{
		axis(1,at=halfDeepCenters,labels=F,line=-0.5,col=halfDeepMarkerColor)                        # interval markers
		axis(1,at=c(-0.2*max(scaffoldTicks),1.2*max(scaffoldTicks)),labels=F,line=-0.5,col="white")  # erase unwanted horizonal 'axis'
		}

	if ((tickSpacing > 0) & (xlim[2]>=tickSpacing))          # equal-spaced ticks
		axis(1,at=seq(tickSpacing,xlim[2],by=tickSpacing),labels=F,col="gray")
	axis(1,at=scaffoldTicks,labels=F,tck=-0.04)              # scaffold ticks
	axis(1,at=scaffoldCenters,tick=F,labels=scaffolds$name,  # scaffold labels
		 las=2,cex.axis=0.7)

	# draw depth as gray, and depth in half-deep intervals as black

	points(depth$s,depth$depth,col=depthColor,pch=19,cex=0.3)

	if (showHalfDeep)
		{
		halfsies = (depth$depth>=halfDepthLo) & (depth$depth<=halfDepthHi)
		points(depth$s[halfsies],depth$depth[halfsies],col=halfDepthColor,pch=19,cex=0.1)
		}

	# draw half-deep intervals as red overlay rectangles; note that R often
	# does a poor job at filling narrow rectangles, so we also draw each
	# rectangle as a centered line; and, because the depth plot will overshoot
	# the upper limit, we multiply that limit by 1.2 here

	if ((showHalfDeep) && (nrow(halfDeep) > 0))
		{
		rect(halfDeep$s,ylim[1],halfDeep$e,ylim[2]*1.2,border=NA,col=halfDeepOverlayColor)  # LBRT

		halfDeepCentersX = matrix(nrow=3,ncol=length(halfDeepCenters))
		halfDeepCentersX[1,] = halfDeepCenters
		halfDeepCentersX[2,] = halfDeepCenters
		halfDeepCentersX[3,] = NA
		dim(halfDeepCentersX) = NULL
		halfDeepCentersY = matrix(nrow=3,ncol=length(halfDeepCenters))
		halfDeepCentersY[1,] = ylim[1]
		halfDeepCentersY[2,] = ylim[2]*1.2
		halfDeepCentersY[3,] = NA
		dim(halfDeepCentersY) = NULL
		lines(halfDeepCentersX,halfDeepCentersY,col=halfDeepOverlayColor)
		}

	# add horizontal lines to show median and half-deep limits

	lines(xlim,c(depth50,depth50),col=depthLimitsColor,lwd=2,lty=2)
	if (showHalfDeep)
		{
		lines(xlim,c(halfDepthHi,halfDepthHi),col=depthLimitsColor,lwd=2,lty=2)
		lines(xlim,c(halfDepthLo,halfDepthLo),col=depthLimitsColor,lwd=2,lty=2)
		}

	text(0,depth50,"median ",adj=1,cex=0.7,col=depthLimitsColor)
	if (showHalfDeep)
		{
		text(0,halfDepthLo,"half-40th ",adj=1,cex=0.7,col=depthLimitsColor)
		text(0,halfDepthHi,"half-60th ",adj=1,cex=0.7,col=depthLimitsColor)
		}

	# close the plot

	if (turnDeviceOff) dev.off()
	}


halfdeep_read_and_plot <- function(lengthsFilename,depthFilename,halfDeepFilename,
                                   percentilesFilename,
                                   plotFilenameTemplate=NULL,
                                   tickSpacing=10000000,
                                   width=17,height=7,pointsize=18,
                                   yLabelSpace=7)
	#
	# Plot half-deep intervals and depth, for several assemblies
	#
	{
	scaffolds = read_scaffold_lengths(lengthsFilename)
	scaffoldToOffset = linearized_scaffolds(scaffolds)
	depth = read_depth(depthFilename,scaffoldToOffset)
	halfDeep = read_halfdeep(halfDeepFilename,scaffoldToOffset)
	percentileToValue = read_percentiles(percentilesFilename)
	
	halfdeep_plot(scaffolds,depth,halfDeep,percentileToValue,assembly,
	              plotFilenameTemplate=plotFilenameTemplate,
	              tickSpacing=tickSpacing,
                  width=width,height=height,pointsize=pointsize,
                  yLabelSpace=yLabelSpace)
	}


read_control_freec <- function(controlFreecFilename,scaffoldToOffset=NULL)
	#
	# Read a file containing the copy number ouput from ControlFREEC.
	#
	# A scaffold-to-offset can be provided, as would be produced by
	# linearized_scaffolds(). One columns are added, converting each window's
	# on-scaffold position to a position along a number line.
	#
	# typical input:
	#	Chromosome Start Ratio   MedianRatio CopyNumber
	#	SUPER_2    1     42.2814 22.2608     45
	#	SUPER_2    1001  31.0994 22.2608     45
	#	SUPER_2    2001  30.5009 22.2608     45
	#	SUPER_2    3001  26.9936 22.2608     45
	#	SUPER_2    4001  33.8109 22.2608     45
	#	SUPER_2    5001  27.2106 22.2608     45
	#    ...
	#
	# returns, e.g.
	#	scaffold start Ratio   MedianRatio CopyNumber s
	#	SUPER_2     1  42.2814 22.2608     45         200529156
	#	SUPER_2  1001  31.0994 22.2608     45         200530156
	#	SUPER_2  2001  30.5009 22.2608     45         200531156
	#	SUPER_2  3001  26.9936 22.2608     45         200532156
	#	SUPER_2  4001  33.8109 22.2608     45         200533156
	#	SUPER_2  5001  27.2106 22.2608     45         200534156
	#    ...
	#
	{
	controlFreec = read.table(controlFreecFilename,header=T,colClasses=c("character","numeric","numeric","numeric","numeric"))
	colnames(controlFreec) <- c("scaffold","start","Ratio","MedianRatio","CopyNumber")
	if (!is.null(scaffoldToOffset))
		{
		controlFreec = controlFreec[controlFreec$scaffold %in% names(scaffoldToOffset),]
		controlFreec[,"s"] = scaffoldToOffset[controlFreec$scaffold] + controlFreec$start
		}

	controlFreec
	}


control_freec_plot <- function(scaffolds,depth,controlFreec,percentileToValue,
                               assemblyName="",
                               scaffoldsToPlot=NULL,
                               plotFilename=NULL,
                               tickSpacing=10000000,
                               tickLabels=F,
                               width=17,height=7,pointsize=18,
                               yLabelSpace=7,
                               maxDepth=NA,
                               scaffoldInterval=NULL)
	#
	# Plot half-deep intervals and depth along an assembly
	#
	# scaffolds:              As returned by read_scaffold_lengths().
	# depth:                  As returned by read_depth().
	# controlFreec:           As returned by read_control_freec(). If this is
	#                         NULL, no controlFreec copy number information is
	#                         displayed.
	# percentileToValue:      As returned by read_percentiles()
	# assemblyName:           Name of the assembly. This contributes to the
	#                         plain title, and can contribute to the plot
	#                         file name.
	#                         Example: "bAlcTor1.pri.cur.20190613"
	# scaffoldsToPlot=NULL:   A vector of scaffolds to restrict the plot to.
	#                         By default, everything in scaffolds[] is
	#                         plotted.
	# plotFilename=NULL:      File to plot the data to. If this is NULL, the
	#                         data is plotted on the screen.
	# tickSpacing:            Spacing of evenly-spaced ticks along the
	#                         horizontal axis. If this is zero, these ticks
	#                         are inhibited. The default is 10Mbp.
	# tickLabels:             If true, add numeric labels to the horizontal
	#                         axis.
	# width,height,pointsize: These are passed through to whatever function
	#                         creates the plot window.
	# yLabelSpace:            Space below the plot. This can be increased to
	#                         accommodate longer scaffold names.
	# maxDepth:               Depth greater than this is not shown in the plot;
	#                         i.e. the vertical axis stops at this value. By
	#                         default this is 1.5*median depth.
	# scaffoldInterval        subinterval to restrict the plot to. Typically
	#                         this would only be used when only one scaffold
	#                         is to be plotted. This is a (start,end) pair,
	#                         origin-zero, half-open.
	#
	{
	# if we have a scaffold subset, reduce our copy of the data to that subset

	if (!is.null(scaffoldInterval))
		{
		if (is.null(scaffoldsToPlot))
			{
			if (length(scaffoldsToPlot) > 1)
				stop("scaffoldInterval cannot be used with more than one scaffold")
			}
		else
			{
			if (length(unique(scaffolds$name)) > 1)
				stop("scaffoldInterval cannot be used with more than one scaffold")
			}
		}

	showControlFreec = !is.null(controlFreec)

	if (is.null(scaffoldsToPlot))
		{
		scaffoldLen = sum(scaffolds$length)
		}
	else
		{
		# validate the names in the subset

		badNames = scaffoldsToPlot[!(scaffoldsToPlot %in% scaffolds$name)]
		if (length(badNames) > 0)
			stop(paste("bad scaffold name(s):",paste(scaffoldsToPlot,collapse=", ")))

		# pick ordered subset
		scaffolds = scaffolds[scaffolds$name %in% scaffoldsToPlot,]
		scaffoldsToNumber = 1:length(scaffoldsToPlot)
		names(scaffoldsToNumber) = scaffoldsToPlot
		scaffolds = scaffolds[order(scaffoldsToNumber[scaffolds$name]),]
		scaffolds[,"offset"] = rev(sum(as.numeric(scaffolds$length))-cumsum(rev(as.numeric(scaffolds$length))))
		scaffoldToOffset = linearized_scaffolds(scaffolds)

		# reduce to ordered subset
		depth = depth[depth$scaffold %in% names(scaffoldToOffset),]
		depth[,"s"] = scaffoldToOffset[depth$scaffold] + depth$start
		depth[,"e"] = scaffoldToOffset[depth$scaffold] + depth$end

		if (showControlFreec)
			{
			controlFreec = controlFreec[controlFreec$scaffold %in% names(scaffoldToOffset),]
			controlFreec[,"s"] = scaffoldToOffset[controlFreec$scaffold] + controlFreec$start
			}

		scaffoldLen = sum(scaffolds$length[scaffolds$name==scaffoldsToPlot])
		}

	# fetch percentile values (used only for drawing and labeling)

	depth50        = percentileToValue["percentile50"]
	halfDepthLo    = percentileToValue["halfPercentile40"]
	halfDepthHi    = percentileToValue["halfPercentile60"]
	depthClip      = ifelse(is.na(maxDepth),1.5*depth50,maxDepth)
	depth50Str     = sprintf("%.1f",depth50)
	halfDepthLoStr = sprintf("%.1f",halfDepthLo)
	halfDepthHiStr = sprintf("%.1f",halfDepthHi)
	depthClipStr   = sprintf("%.1f",depthClip)

	# (housekeeping)

	scaffoldTicks = 1 + c(scaffolds$offset,sum(as.numeric(scaffolds$length)))
	scaffoldCenters = (as.numeric(scaffoldTicks[1:nrow(scaffolds)])+as.numeric(scaffoldTicks[2:(nrow(scaffolds)+1)])) / 2

	CNSpacing = depthClip / 32

	if (!is.null(scaffoldInterval))
		{
		xlim = scaffoldInterval
		depth = depth[(depth$s>=xlim[1])&(depth$s<=xlim[2]),]
		controlFreec = controlFreec[(controlFreec$s>=xlim[1])&(controlFreec$s<=xlim[2]),]
		scaffoldCenters = (scaffoldInterval[1]+scaffoldInterval[2])/2
		}
	else
		{
		xlim = c(1,max(scaffoldTicks))
		}

	ylim = if (showControlFreec) c(-5*CNSpacing,depthClip) else c(0,depthClip)
	ylimLong = ylim
	ylimLong[1] = ylimLong[1] - CNSpacing

	depthColor       = rgb(.6,.6,.6)
	halfDepthColor   = "black"
	depthLimitsColor = "blue"
	clippedColor     = "red"

	# clipping

	depth$clipped = ifelse(depth$depth<=depthClip,depth$depth,depthClip)
	depth$color   = ifelse(depth$depth<=depthClip,depthColor,clippedColor)

	controlFreec$CNclipped = ifelse(controlFreec$CopyNumber<=2,controlFreec$CopyNumber,3)

	# open plot window or file

	turnDeviceOff = F
	if (is.null(plotFilename))
		{
		quartz(width=width,height=height,pointsize=pointsize)
		}
	else
		{
		print(paste("drawing to",plotFilename))
		pdf(file=plotFilename,width=width,height=height,pointsize=pointsize)
		turnDeviceOff = T
		}

	# create empty plot

	if (assemblyName == "")
		title = paste("Copy Number\n(per ControlFREEC)",sep="")
	else
		title = paste("Copy Number on ",assemblyName,"\n(per ControlFREEC)",sep="")
	ylab = "depth (black/gray) and CN (blue/red)"

	par(mar=c(yLabelSpace,4,2.5,0.2)+0.1)     # BLTR
	options(scipen=10)
	plot(NA,xlim=xlim,ylim=ylim,main=title,xaxt="n",xlab="",ylab=ylab)

	# add horizontal axis

	if ((tickSpacing > 0) & (xlim[2]>=tickSpacing))              # equal-spaced ticks
		{
		leftTick = xlim[1] + (tickSpacing-1) - ((xlim[1] + (tickSpacing-1)) %% tickSpacing)
		ticks = seq(leftTick,xlim[2],by=tickSpacing)
		axis(1,at=ticks,labels=F,col="gray")
		if (tickLabels)
			{
			labeledTicks = ticks[abs(scaffoldCenters-ticks)>=tickSpacing/4]
			axis(1,at=labeledTicks,tick=F,                       # tick labels
			     labels=prettyNum(labeledTicks,big.mark=",",scientific=FALSE),
				 las=2,pos=ylim[1]-0.5,cex.axis=0.5)
			}
		}
	axis(1,at=scaffoldTicks,labels=F,tck=-0.08)                  # scaffold ticks
	axis(1,at=scaffoldCenters,tick=F,labels=scaffolds$name,      # scaffold labels
		 las=2,cex.axis=0.7)

	# draw depth as gray, and depth in half-deep intervals as black

	points(depth$s,depth$clipped,col=depth$color,pch=19,cex=0.3)

	halfsies = (depth$depth>=halfDepthLo) & (depth$depth<=halfDepthHi)
	points(depth$s[halfsies],depth$depth[halfsies],col=halfDepthColor,pch=19,cex=0.1)

	for (ix in 1:length(scaffoldTicks))
		lines(c(scaffoldTicks[ix],scaffoldTicks[ix]),ylimLong,col="black",lwd=1,lty=2)

	# add horizontal lines to show median and half-deep limits

	lines(xlim,c(depth50,depth50),col=depthLimitsColor,lwd=2,lty=2)
	lines(xlim,c(halfDepthHi,halfDepthHi),col=depthLimitsColor,lwd=2,lty=2)
	lines(xlim,c(halfDepthLo,halfDepthLo),col=depthLimitsColor,lwd=2,lty=2)

	text(xlim[1],depth50,    "median ",   adj=1,cex=0.5,col=depthLimitsColor)
	text(xlim[1],halfDepthLo,"half-40th ",adj=1,cex=0.5,col=depthLimitsColor)
	text(xlim[1],halfDepthHi,"half-60th ",adj=1,cex=0.5,col=depthLimitsColor)

	# add controlFreec copy number information
	# from bottom up, rows are copy number = 0, 1, 2, >2

	if (showControlFreec)
		{
		for (cn in -2:-5)
			{
			lines(xlim,c(cn*CNSpacing,cn*CNSpacing),col="red",lty=1)
			text(xlim[1]-0.005*(xlim[2]-xlim[1]),cn*CNSpacing,cex=0.4,adj=1,
			     ifelse(cn==-2,"CN>2",paste("CN=",cn+5,sep="")))
			}
		points(controlFreec$s,(controlFreec$CNclipped-5)*CNSpacing,pch=16,cex=0.5,
			   col=ifelse(controlFreec$CNclipped==1,"blue","red"))
		}

	# close the plot

	if (turnDeviceOff) dev.off()
	}

