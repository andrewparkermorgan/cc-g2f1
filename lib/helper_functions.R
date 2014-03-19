#############################################
## helper_functions.R
# purpose:	functions to support CC G2:F1 genetic map analyses 
# author:		GAC 12 Dec 2011
# udpated:	GAC 10 Jan 2012
#############################################


##############################
## make.density()
# function to convert cumulative map to densities on a fixed grid
# compute crossover density from cumulative map
# dscale is in kb
# function is designed for use with ddply
# e.g.  Den <- ddply(Map, .(Chr), 1000, 1/100)
make.density <- function(Map, dstep= 100, dscale=1000){
	  
  Den <- crossoverDensity(Map,dstep,dscale)

  #rescale event counts and tally events per Mb by Chromosome
  Den$Events <- 1000 * dscale * Den$Events
  Den$cM <- 1000 * dscale * Den$cM
  #
  #trim the first 3Mb off each Chr
  # Den <- Den[Den$Pos>3000000+1000*dscale/2, ]
  #
  # set Chr to factor and tidy up
  Den$Chr <- as.factor(Den$Chr)
  row.names(Den) <- NULL
  Den  
}

##############################
## crossoverDensity()
# Calculates the density from the cumulative map based on a sliding window.
# The returned value will be a dataframe with two columns: 
#   Pos is the window Center in bp, 
#   and columns in addition to Chr and bp, will be converted to a density 
#   over the window by taking the difference in cumulative values
#   between the startbp and endbp
#  note: it might be better to to return startbp and endbp instead of Pos
#
# map:
#   a cumulative distribution as created by organizeByChrom or
#   integrateCrossovers.  This should be a dataframe indexed by chromosome (Chr)
#   and a base-pair position (Pos)
#   
# stepSizeKb:
#   how far forward should we move the window for each step in kilobases.
#   This parameter allows for a trade-off between resolution (how smooth the
#   density curve is) and speed
#
# windowSizeKb:
#   how wide should the window be in kilobases
#
crossoverDensity <- function(Map, stepSizeKb = 100, windowSizeKb = 1000) {    
	stepSizeBp <- stepSizeKb * 1000
    windowSizeBp <- windowSizeKb * 1000
    
    densityMat <- NULL
    if(nrow(Map) >= 2) {
            startBp <- Map[1, "Pos"]
            endBp <- Map[nrow(Map), "Pos"]
            winStartsBp <- seq(startBp, endBp - windowSizeBp, stepSizeBp)
            winEndsBp <- winStartsBp + windowSizeBp
            winStartsCm <- .chrInterpolateFromBp(Map[,c("Pos","cM")], winStartsBp)
            winEndsCm <- .chrInterpolateFromBp(Map[,c("Pos","cM")], winEndsBp)
            
            winStartsEvents <- .chrInterpolateFromBp(Map[,c("Pos","cumEvents")], winStartsBp)
            winEndsEvents <- .chrInterpolateFromBp(Map[,c("Pos","cumEvents")], winEndsBp)
            
            
            Pos <- (winStartsBp + winEndsBp) / 2
            
            #Density <- ddply(CoxMap.avg, .(Chr), numcolwise(max))[,-1])

            cM <- (winEndsCm - winStartsCm) / windowSizeBp
            Events <- (winEndsEvents - winStartsEvents) / windowSizeBp
            #densityMat <- as.data.frame(cbind(Pos, Density[,-c(1,2)]))
            densityMat <- as.data.frame(cbind(Map$Chr[1], Pos, Events, cM))
            names(densityMat)[1] <- "Chr"
        }
    densityMat
}
#
.chrInterpolateFromBp <- function(Map, bpPositions) {
    foundIndices <- findInterval(bpPositions, Map[, 1])
    
    interpPositions <- c()
    mapRowCount <- nrow(Map)
    for(posIndex in 1 : length(bpPositions))
    {
        pos <- bpPositions[posIndex]
        foundIndex <- foundIndices[posIndex]
        if(foundIndex == 0) {
            interpPositions[posIndex] <- NA
        } else {
            foundBpPos <- Map[foundIndex, 1]
            foundCmPos <- Map[foundIndex, 2]
            if(pos == foundBpPos) {
                interpPositions[posIndex] <- foundCmPos
            } else {
                nextBpPos <- Map[foundIndex + 1, 1]
                nextCmPos <- Map[foundIndex + 1, 2]
                bpGap <- nextBpPos - foundBpPos
                cmGap <- nextCmPos - foundCmPos
                
                gapWeight <- (pos - foundBpPos) / bpGap
                interpPositions[posIndex] <- foundCmPos + (cmGap * gapWeight)
            }
        }
    }
    interpPositions
}

###################################################
# excursion()
# recursion function to operate on the precomputed e.scores in a Den dataframe
# computes a 1-D Smith-Waterman forward pass on e.scores
# and returns a dataframe with new element "E"
#
excursion <- function(Den){
  E <- rep(0, length(Den$e.score))
  E[1] <- max(0,Den$e.score[1])
  for(i in 2:length(Den$e.score)){
		 E[i] <- max(0, E[i-1]+Den$e.score[i])
	}
  Den$E <- E
  Den
}

###################################################
#make.excursion()
  #inputs
  #Den is a dnsity map dataframe with elements Chr, Pos, Events, cM
  #  as created by make.density()
  #dscale is the window size for computing the density estimate
  #mult is the enrichment factor for hotspot (e.g., 100) or coldspot (e.g., 1/100)
  #
  #return value is a new dataframe with added columns
  # e.score is a score of enrichment for each window
  # E is the result of Smith-Waterman forward pass = excursion score
  # function is designed for use with ddply
  # e.g.  Den <- ddply(Map, .(Chr), 1000, 1/100)
#
make.excursion <- function(Den, mult=100){
  # compute the e.scores
  lambda <- mean(Den$Events)
  gamma <- lambda * mult
  e.score <- lambda - gamma + Den$Events * log(gamma/lambda)
  Den <- transform(Den, e.score=e.score)
  
  # compute the forward pass
  Den <- excursion(Den)
  Den
} 

#################################################
#get.peaks()
# extract best segment(s) from an excursion plot
# takes a list of scores and the output of the 1-D forward pass
# it extracts the best segment and resets the scores over that segment to zero
# following extraction of the best segment, the forward pass can be re-executed 
# and the next segment extracted...
#
get.peaks <- function(Den,Thresh=10){ 
  peaks <- NULL
  while(max(Den$E)>Thresh){
    # find index of the maximum excursion
    indx.end <- max(which(Den$E == max(Den$E)))
    #
    #find the last zero preceeding the max index
    zero.list <- which(Den[1:indx.end,]$E == 0)
    indx.start <- ifelse(length(zero.list)==0,1,max(zero.list))

    # report region
    tmp <- cbind(Den[indx.start,c("Chr","Pos")], Den[indx.end,c("Pos","E")])
    names(tmp) <- c("Chr","Pos.start","Pos.end","E")
    peaks <- rbind(peaks, tmp)
    #
    #reset the scores to zero
    Den[indx.start:indx.end,"e.score"] <- 0
    #Den$E <- 0
    #
    # recompute the forward pass
    Den <- excursion(Den)
  }
  if(length(peaks$Pos.start)>1){
    peaks <- peaks[order(peaks$Pos.start),]
  }
  peaks
}

########################################
#  integrateCrossovers()
# Given the Events Data (in interval format Chr, Start, End) 
# this function computes a cumulative distribution of events
# output is a dataframe:
#    Pos = base pair position 
#    Events = accumulated (integrated) number of events
# with one row for every event boundary
# input is a dataframe with three columns:
#    chromosomeID, crossoverStartBp, crossoverEndBp
# designed to be called usinng ddply e.g., 
# Map <- ddply(Events, .(Chr))
#
integrateCrossovers <- function(Events, Chr=NULL) {
#   

	#trim the Events data and compute interval lengths 
	Events <- Events[,c("Start","End")]
	Events <- transform(Events, Length=End-Start)
	
	#temporary fix
	Events <- subset(Events, Start < End)
	
	if(any(Events$Length <= 0)){
		stop("Zero or Negative Interval Length")
	}

	#find all unique event boundaries
	Pos <- sort(unique( c(Events$Start, Events$End) ))
	
	#initialize cumulative events and count	
 	cumEvents <- rep(0,length(Pos))
	for(j in 2:length(Pos)){
		cumEvents[j] <- cumEvents[j-1] + 
			sum(((Pos[j]-Pos[j-1])/subset(Events, Start<=Pos[j-1] & End>=Pos[j])$Length))	
	}

    as.data.frame(cbind(Pos, cumEvents))
}

#####################
# pad cumulative events 
# assume sorted cumulative map with fields Chr and Pos
padMap <- function(Map)
{
	Pos.end <- c(197111960,181687110,159598996,155477512,152295665,149504681,152450636,
			131570700,124027719,129936781,121754330,121256667,120253871,125150241,
			103458457,98276177,95264893,90711723,61332329,166365178)

	Pos.start <- c(3036178,3011992,3033486,3270810,3073698,3001551,3073581,3083611,3084105,3009197,
			3113512,3081693,3006383,8527441,3096128,3577741,3070591,3012495,3125547,5023070)
	
	newMap <- NULL
	
	for(i in unique(Map$Chr)){
		tmpMap <- subset(Map, Chr==i)
		firstrow <- tmpMap[1,]
		firstrow$Pos <- Pos.start[i]
		lastrow <- tmpMap[dim(tmpMap)[1],]
		lastrow$Pos <- Pos.end[i]
		#
		newMap <- rbind(newMap,firstrow,tmpMap,lastrow)
	}
	newMap	
}

##########################################################
# .chrNameValues()
# returns integer values which can be used to order the given chromosome
# strings (ie the returned vector is suitable as an input to the order(...)
# function)
.chrNameValues <- function(chrNames) {
    origNames <- chrNames
    
    chrNames <- toupper(chrNames)
    chrNames <- sub("(CHR|CHROMOSOME)", "", chrNames)
    chrInt <- suppressWarnings(as.integer(chrNames)) 
    chrInt[chrNames == "X"] <- .Machine$integer.max - as.integer(2)
    chrInt[chrNames == "Y"] <- .Machine$integer.max - as.integer(1)
    chrInt[chrNames == "M"] <- .Machine$integer.max
    
    if(any(is.na(chrInt))) {
        stop("unrecognized chromosome name(s): ",
                paste(unique(origNames[is.na(chrInt)]), collapse = ", "))
    }
    
    chrInt
}

##########################################################
# .make.bin.labels()
# helper function to make friendly bin labels (eg. for histograms)
make.bin.labels <- function(breaks, right.closed = TRUE) {
	breaks <- sort(breaks)
	if (length(breaks) < 2) {
		return(NA)
	}
	labels <- character()
	while(length(breaks) > 1) {
		if (!is.finite(breaks[1])) {
			breaks <- breaks[-1]
			leq <- ifelse(right.closed, "<=", "<")
			labels <- c(labels, paste(leq, breaks[1]))
		}
		else {
			if (!is.finite(breaks[2])) {
				geq <- ifelse(right.closed, ">", ">=")
				labels <- c(labels, paste(geq, breaks[1]))
				breaks <- breaks[ -length(breaks) ]
			}
			else {
				paren.open <- ifelse(right.closed, "(", "[")
				paren.close <- ifelse(right.closed, "]", ")")
				labels <- c(labels, paste(paren.open, breaks[1], ", ", breaks[2], paren.close, sep = ""))
			}
		}
		breaks <- breaks[-1]
	}
	return(labels)
}
