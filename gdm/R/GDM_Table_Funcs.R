
gdm <- 
function (data, geo=FALSE, splines=NULL, knots=NULL) 
{
	options(warn.FPU = FALSE)
##
##	sanity check on the data table
##
	if ( ncol(data) < 6 )
	    stop("Not enough columns in data. (At minimum need: Observed, weights, X0, Y0, X1, Y1) ")
	if ( nrow(data) < 1 )
	    stop("Not enough rows in data")


##
##	check that the response data is [0..1]
##
	rtmp <- data[,1]
	if (length(rtmp[rtmp<0]) > 0)
 	    stop("Response data has negative values. Must be between 0 - 1.")
	if (length(rtmp[rtmp>1]) > 0)
	    stop("Response data has values greater than 1. Must be between 0 - 1.")


##
##	current data format is response,weights,X0,Y0,X1,Y1 before any predictor data (thus 6 leading columns)
##
	LEADING_COLUMNS <- 6
	if ( geo ) 
        {
	    nPreds <- ( ncol(data) - LEADING_COLUMNS ) / 2 + 1		
	}
	else 
        {
	    nPreds <- ( ncol(data) - LEADING_COLUMNS ) / 2
	}

	if ( nPreds < 1 )
	    stop("Data has NO predictors")


##
## 	setup the predictor name list
##
	if ( geo ) 
        {
	    if ( nPreds > 1 )
		predlist <- c("Geographic", names(data)[(LEADING_COLUMNS+1):(LEADING_COLUMNS+nPreds-1)])
	    else
		predlist <- c("Geographic")
	}
	else 
        {
  	    predlist <- names(data)[(LEADING_COLUMNS+1):(LEADING_COLUMNS+nPreds)]
	}


##
##	deal with the splines and knots
##
	if (is.null(knots))
	{
	    ##
	    ## generate knots internally from the data
	    ##
	    if ( is.null(splines) ) 
            {
		nSplines <- 3
		quantvec <- rep(0, nPreds * nSplines)
		splinvec <- rep(nSplines, nPreds)

		if ( geo ) 
                {
      		    ## get knots for the geographic distance
            	    v <- sqrt((data[,3]-data[,5])^2 + (data[,4]-data[,6])^2)
		    quantvec[1] <- min(v)
      		    quantvec[2] <- median(v)
            	    quantvec[3] <- max(v)

		    if ( nPreds > 1 ) 
                    {
		        ## get knots for the environmental predictors
			for (i in seq(from = 1, to = nPreds-1, by = 1)) 
                        {
            		    v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds-1])                 
      	            	    index = i * nSplines
		            quantvec[index+1] <- min(v)
		      	    quantvec[index+2] <- median(v)
      	      	      	    quantvec[index+3] <- max(v)
	      	      	}
		    }
    		}
		
		else 
                {
		    ## get knots for the environmental predictors after skipping geographic preds
	      	    for (i in seq(from = 1, to = nPreds, by = 1)) 
                    {
			v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds])                   
            	      	index = (i-1) * nSplines
	            	quantvec[index+1] <- min(v)
	      	        quantvec[index+2] <- median(v)
            	      	quantvec[index+3] <- max(v)
		    }
		}
	    }
		
	    else 
            {
		##
		##	otherwise check that the supplied splines vector has enough data and minumum spline values of 3
		##
		if ( length(splines) != nPreds ) 
                {
		    stop(paste("Number of splines does not equal the number of predictors. 
                   Splines argument has", length(splines), "items but needs", nPreds, "items."))
		}

  		## count the total number of user defined splines to dimension the knots vector
		quantvec <- rep(0, sum(splines))
		splinvec <- splines

		if ( geo ) 
                {
		    if ( splines[1] < 3 )
			stop("Must have at least 3 splines per predictor")

                    ## get knots for the geographic distance
      	      	    v <- sqrt((data[,3]-data[,5])^2 + (data[,4]-data[,6])^2)
		    quantvec[1] <- min(v)		## 0% knot
      	      	    quantvec[splines[1]] <- max(v)	## 100% knot

		    quant_increment <- 1.0 / (splines[1]-1)
		    this_increment <- 1
		    for (i in seq(from = 2, to = (splines[1]-1), by = 1)) 
                    {  
		        ## mid % knots
			quantvec[i] <- quantile(v,quant_increment*this_increment)
			this_increment <- this_increment + 1
		    }

		    if ( nPreds > 1 ) 
                    {
			## get knots for the environmental predictors
			current_quant_index <- splines[1]
			for (i in seq(from = 1, to = nPreds-1, by = 1)) 
                        {
      	      		    num_splines <- splines[i+1]
			    if ( num_splines < 3 )
			        stop("Must have at least 3 splines per predictor")

			    v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds-1])                 
			    quantvec[current_quant_index+1] <- min(v)	            ## 0% knot
	      		    quantvec[current_quant_index+num_splines] <- max(v)	    ## 100% knot

			    quant_increment <- 1.0 / (num_splines-1)
			    this_increment <- 1
			    for (i in seq(from = 2, to = (num_splines-1), by = 1)) 
                            {  
				## mid % knots
				quantvec[current_quant_index+i] <- quantile(v,quant_increment*this_increment)
				this_increment <- this_increment + 1
			    }
			
			    current_quant_index <- current_quant_index + num_splines
			}
		    }
                }

		else 
                {
		    ## get knots for the environmental predictors
		    current_quant_index <- 0
		    for (i in seq(from = 1, to = nPreds, by = 1)) 
                    {
      	      	        num_splines <- splines[i]
			if ( num_splines < 3 )
			    stop("Must have at least 3 splines per predictor")

                        v <- c(data[,i+LEADING_COLUMNS], data[,i+LEADING_COLUMNS+nPreds])          
			quantvec[current_quant_index+1] <- min(v)	        ## 0% knot
	      		quantvec[current_quant_index+num_splines] <- max(v)	## 100% knot

			quant_increment <- 1.0 / (num_splines-1)
			this_increment <- 1
			for (i in seq(from = 2, to = (num_splines-1), by = 1)) 
                        {  
			    ## mid % knots
			    quantvec[current_quant_index+i] <- quantile(v,quant_increment*this_increment)
			    this_increment <- this_increment + 1
			}
			current_quant_index <- current_quant_index + num_splines
		    }
		}
	    }
	}

	else
	{
	    ##
	    ## user defined knots supplied as an argument
	    ##
	    if ( is.null(splines) ) 
            {
	        ## check that there are nPreds * 3 knots in the user defined vector
		if ( length(knots) != (nPreds * 3) ) 
                {
		    stop(paste("When knots are supplied by the user, there should be", (nPreds * 3), "items in the knots argument, not", length(knots), "items."))
		}

		## now check that each of the three knots for each predictor are in ascending order
		for (i in seq(from = 1, to = nPreds, by = 1)) 
                {
                    index = i * 3
		    if ((knots[index-1] < knots[index-2] ) || 
                        (knots[index] < knots[index-2]) || 
                        (knots[index] < knots[index-1])) 
                    {
			stop(paste("Knots for ", predlist[i], "are not in ascending order."))
		    }
		}

                nSplines <- 3
		quantvec <- knots
		splinvec <- rep(nSplines, nPreds)
	    }
		
	    else 
            {
		## check that there are sum(splines) knots in the user defined vector
		if ( length(knots) != sum(splines) ) 
                {
		    stop(paste("When knots are supplied by the user, there should be", sum(splines), "items in the knots argument, not", length(knots), "items."))
		}
				
		## now check that each of the knots for each predictor are in ascending order
		index = 0
		for (i in seq(from = 1, to = nPreds, by = 1)) 
                {
		    for (j in seq(from = 2, to = splines[i], by = 1)) 
                    {
                        if (knots[index+j] < knots[index+j-1]) 
                        {
			    stop(paste("Knots for ", predlist[i], "are not in ascending order."))
			}
		    }
		    index <- index + splines[i]
		}

		quantvec <- knots
		splinvec <- splines
	    }
	}

      	p1 <- 0
	p2 <- 0
    	p3 <- 0
    	p4 <- 0
    	p5 <- rep(0,times=length(quantvec))
      	p6 <- rep(0,times=nrow(data))
      	p7 <- rep(0,times=nrow(data))
      	p8 <- rep(0,times=nrow(data))
	
##
##  Call the dll function
##
	z <- .C( "GDM_FitFromTable",
                 paste(getwd()),
                 as.matrix(data),
                 as.integer(geo),
                 as.integer(nPreds), 
                 as.integer(nrow(data)), 
                 as.integer(ncol(data)),
                 as.integer(splinvec),
                 as.double(quantvec),                 
                 gdmdev = as.double(p1),
                 nulldev = as.double(p2),
                 expdev = as.double(p3),
                 intercept = as.double(p4),         
                 coeffs = as.double(p5),
                 response = as.double(p6),
                 preddata = as.double(p7),
                 ecodist = as.double(p8), 
                 PACKAGE = "gdm"
               )

##	call <- match.call()
        m <- match.call(expand.dots = F)
        
    	gdmModOb <- structure(list(dataname = m[[2]],
                       geo = geo,
                       sample = nrow(data),
                       gdmdeviance = z$gdmdev,
                       nulldeviance = z$nulldev,
                       explained = z$expdev,
                       intercept = z$intercept,
                       predictors = predlist,
                       coefficients = z$coeffs,
                       knots = quantvec,
                       splines = splinvec,
                       creationdate = date(),
                       observed = z$response,
                       predicted = z$preddata,
                       ecological = z$ecodist))

        class(gdmModOb) <- c("gdm", "list")
        return(gdmModOb)
}



plot.gdm <- 
function (model, plot.layout = c(2,2), plot.color = "blue", plot.linewidth=2.0) 
{
      	options(warn.FPU = FALSE)
      	PSAMPLE <- 200
      	preddata <- rep(0,times=PSAMPLE)
			
	##
	## establish what plot layout to use
	##  
	thisplot <- 0
	one_page_per_plot <- FALSE
	if ((plot.layout[1]==1) && (plot.layout[2]==1)) 
		one_page_per_plot <- TRUE
	else
	  par(mfrow=plot.layout)	


      	##
      	## apply the link function and plot.....
      	##
      	plot( model$ecological, 
              model$observed, 
              xlab="Predicted Ecological Distance", 
              ylab="Observed Compositional Dissimilarity", type="n" )

        points( model$ecological, model$observed, pch=20, cex=0.25, col=plot.color )

        overlayX <- seq( from=min(model$ecological), to=max(model$ecological), length=PSAMPLE )
        overlayY <- 1 - exp( - overlayX )
        lines( overlayX, overlayY, lwd=plot.linewidth ) 
	thisplot <- thisplot + 1


        ##
        ## use the raw data and plot.....
        ##
	if (one_page_per_plot)
        {
      		dev.new()
        	dev.next()
	}
        plot( model$predicted, 
              model$observed, 
              xlab="Predicted Compositional Dissimilarity", 
              ylab="Observed Compositional Dissimilarity", type="n" )

        points( model$predicted, model$observed, pch=20, cex=0.25, col=plot.color )

        overlayX <- overlayY <- seq( from=min(model$predicted), to=max(model$predicted), length=PSAMPLE )
        lines( overlayX, overlayY, lwd=plot.linewidth ) 
	thisplot <- thisplot + 1


        ##
        ## determine the max of all the predictor data
        ##
        preds <- length(model$predictors)
        predmax <- 0
        splineindex <- 1
        for ( i in 1:preds ) 
        {  
            ## only if the sum of the coefficients associated with this predictor is > 0.....
            numsplines <- model$splines[i]
            if ( sum(model$coefficients[splineindex:(splineindex+numsplines-1)]) > 0 ) 
            {
	        ## get predictor plot Y-data                            
                z <- .C( "GetPredictorPlotData", 
                         pdata = as.double(preddata),
                         as.integer(PSAMPLE),
                         as.double(model$coefficients[splineindex:(splineindex+numsplines-1)]),
                         as.double(model$knots[splineindex:(splineindex+numsplines-1)]),
                         as.integer(numsplines), 
                         PACKAGE = "gdm"
                       )
		                        
		v <- max(z$pdata)
                if (v > predmax ) predmax <- v
            }
            splineindex <- splineindex + numsplines
	}

	
        ##
        ## plot the predictors with non-zero sum of coefficients
        ##      
        splineindex <- 1
	for ( i in 1:preds ) 
        {  
            ## only if the sum of the coefficients associated with this predictor is > 0.....
            numsplines <- model$splines[i]
            if ( sum(model$coefficients[splineindex:(splineindex+numsplines-1)]) > 0 ) 
            {
		if (one_page_per_plot)
                {
      		    dev.new()
		    dev.next()
		}
		
   	        else 
                {
                    thisplot <- thisplot + 1
		    if (thisplot > (plot.layout[1] * plot.layout[2])) 
                    {	
			dev.new()                              				
                        dev.next()					
			thisplot <- 1
                        par(mfrow=plot.layout)	
		    }
		}
                  
		## get predictor plot Y-data    
                 z <- .C( "GetPredictorPlotData", 
                           pdata = as.double(preddata),
                           as.integer(PSAMPLE),
                           as.double(model$coefficients[splineindex:(splineindex+numsplines-1)]),
                           as.double(model$knots[splineindex:(splineindex+numsplines-1)]),
                           as.integer(numsplines),
                           PACKAGE = "gdm"
                        )

                        
                  plot( seq(from=model$knots[[(i*3)-2]],to=model$knots[[(i*3)]], length=PSAMPLE),
                  	z$pdata, 
                        xlab=model$predictors[i], 
                        ylab=paste("f(", model$predictors[i], ")", sep="" ), 
            	  	ylim=c(0,predmax), type="l" )
      	    }
	      	splineindex <- splineindex + numsplines
      	}       
}


predict.gdm <- 
function (model, data) 
{
    options(warn.FPU = FALSE)

#
#  Call the dll function
#
    predicted <- rep(0,times=nrow(data))
    z <- .C( "GDM_PredictFromTable",
             as.matrix(data),
             as.integer(model$geo),
             as.integer(length(model$predictors)), 
             as.integer(nrow(data)), 
             as.double(model$knots),
             as.integer(model$splines),
             as.double(c(model$intercept,model$coefficients)),
             preddata = as.double(predicted),
             PACKAGE = "gdm"
           )
    return(z$preddata)
}


gdm.transform <- 
function (model, data) 
{
  options(warn.FPU = FALSE)
  transformed <- matrix(0, nrow(data), ncol(data))
  colnames(transformed) <- names(data)
  z <- .C("GDM_TransformFromTable", as.integer(nrow(data)), 
          as.integer(ncol(data)), as.integer(model$geo), as.integer(length(model$predictors)), 
          as.integer(model$splines), as.double(model$knots), as.double(model$coefficients), 
          as.matrix(data), trandata = as.double(transformed), 
          PACKAGE = "gdm")
  nRows <- nrow(data)
  nCols <- ncol(data)
  myVec <- z$trandata
  pos <- 1
  for (i in seq(from = 1, to = nCols, by = 1)) {
    tmp <- myVec[seq(from = pos, to = pos + nRows - 1)]
    transformed[, i] <- tmp
    pos <- pos + nRows
  }
  return(transformed)
}



summary.gdm <- 
function (model) 
{
        print( "", quote=F )    
        print( "", quote=F )    
        print( "GDM Modelling Summary", quote=F );
        print( paste( "Creation Date: ", model$creationdate ), quote=F );
        print( "", quote=F )    
##        call <- match.call()
        m <- match.call(expand.dots = F)
        print( paste( "Name: ", m[[2]] ), quote=F )
        print( "", quote=F )    
        print( paste( "Data: ", model$dataname ), quote=F )
        print( "", quote=F )    
        print( paste( "Samples: ", model$sample ), quote=F )
        print( "", quote=F )    
        print( paste( "Geographical distance used in model fitting? ", model$geo ), quote=F )
        print( "", quote=F )    
        print( paste( "NULL Deviance: ", model$nulldeviance ), quote=F )
        print( paste( "GDM Deviance: ", model$gdmdeviance ), quote=F )  
        print( paste( "Deviance Explained: ", model$explained ), quote=F )
        print( "", quote=F )    
        print( paste( "Intercept: ", model$intercept ), quote=F )
        print( "", quote=F )    
        thiscoeff <- 1
        thisquant <- 1
        for ( i in 1:length(model$predictors) ) {
        print( paste( "Predictor ",i,": ",model$predictors[[i]], sep="" ), quote=F )            
        print( paste( "Splines: ",model$splines[[i]], sep="" ), quote=F )
        numsplines <- model$splines[[i]]
        for ( j in 1:numsplines ) {
                if ( j == 1 ) print( paste( "Min Knot: ",model$knots[[thisquant]], sep="" ), quote=F )          
                else if ( j == numsplines ) print( paste( "Max Knot: ",model$knots[[thisquant]], sep="" ), quote=F )
                else print( paste( round(100/(numsplines-1),digits=2),"% Knot: ",model$knots[[thisquant]], sep="" ), quote=F )
                thisquant <- thisquant + 1
        }
        for ( j in 1:numsplines ) {
            print( paste( "Coefficient[",j,"]: ",model$coefficients[[thiscoeff]], sep="" ), quote=F )
            thiscoeff <- thiscoeff + 1
        }
        print( "", quote=F )                
    }   
}



formatsitepair <- function(bioData, bioFormat, dist="bray", abundance=F, 
                           siteColumn=NULL, XColumn, YColumn, sppColumn=NULL, 
                           abundColumn=NULL, sppFilter=0, predData, distPreds=NULL, 
                           weightType="equal", custWeightVect=NULL, samples=NULL){
  ##Takes species data from a variety of commonly used formats and transforms it into the
  ##required site-pair format for GDM
  ##Input Variables:
  ##bioData = the input species data for gdm
  ##bioFormat = the format of the species data table
  ##dist = the dissimilarity distance metric, defined by the vegdist function in the vegan package
  ##abundance = whether the species data are (F) pres/abs or (T) abundance
  ##XColumn = the x coordinates of the sample sites
  ##YColumn = the y coordinates of the sample sites
  ##siteColumn = the name of the column in the species data which holds the site identification, needed for table type 1
  ##abundColumn = the name of the column which holds the species data, needed for table type 2
  ##sppColumn = species code column, like name, needed for table type 2
  ##sppFilter = number of spp at a site for that site to be included in model fitting 
  ##predData = environmental data to be used in model fitting
  ##distPreds = a list of user prepared distance matrices to be used in model fitting
  ##weightType = how the weights of the site-pair table are calculated
  ##custWeightVect = custom weights vector
  ##samples = the maximum number of sites pairs to select at randomly from the full table (to reduce computational demands) 
  ##Output Variables:
  ##outTable = the fully calculated site-pair table for GDM
  ###########################
  #bioData <- commM
  #bioFormat <- 1
  #dist <- "bray"
  #abundance <- F
  #siteColumn <- "ID"
  #XColumn <- "lon"
  #YColumn <- "lat"
  #sppColumn <- NULL
  #sppFilter <- 0
  #abundColumn <- NULL
  #predData <- clim
  #distPreds <- NULL
  #weightType <- "equal"
  #custWeightVect <- NULL
  #samples <- NULL
  #################
  #bioData <- testData1a
  #bioFormat <- 1
  #dist <- "bray"
  #abundance <- F
  #siteColumn <- "site"
  #XColumn <- "Long"
  #YColumn <- "Lat"
  #sppColumn <- NULL
  #sppFilter <- 0
  #abundColumn <- NULL
  #predData <- envTab
  #distPreds <- NULL
  #weightType <- "equal"
  #custWeightVect <- NULL
  #samples <- NULL
  ###########################
  ##required libraries
  #require(raster)
  #require(reshape2)
  #require(vegan)
  #print(class(predData))
  ##if rows of bioData and predData do not match, exit function
  ##if bioFormat is not a number, exit function
  if(bioFormat %in% c(1:4)){} else{
    stop("Acceptable values for the bioFormat argument are: 1, 2, 3, or 4")
  }
  ##if samples is not a number, then exit function
  if(is.numeric(samples)==FALSE & is.null(samples)==FALSE){
    stop("samples argument must be a number")
  }
  ##makes sure that sppFilter is a number, if not exit function
  if(is.numeric(sppFilter)==FALSE){
    stop("sppFilter argument must be a number")
  }
  
  if(weightType %in% c("equal", "richness", "custom")){} else{
    stop("Acceptable values for the weightType argument are: equal, richness, or custom")
  }
  
  if(bioFormat==2 & is.null(siteColumn)==T){
    if(!(class(predData)=="RasterStack" | class(predData)=="RasterLayer" | class(predData)=="RasterBrick")){
      stop("A siteColumn needs to be provided in either the bioData or predData inputs")
    }
  }
  
  if(is.null(siteColumn)==F){
    if(!(siteColumn %in% colnames(bioData)) & (bioFormat==1 | bioFormat==2)){
      stop("Cannot find a match for siteColumn in the columns of bioData.")
    }
    if(bioFormat==3 & siteColumn %in% colnames(bioData)){
      wSite <- which(colnames(bioData)==siteColumn)
      bioData <- bioData[-wSite]
    }
  }
  
  if(bioFormat==3){
    if(weightType=="richness"){
      stop("Cannot weight by site richness when supplying the biological data as a distance matrix.")
    }
    if(nrow(bioData)!=ncol(bioData)){
      stop("Biological dissimularity has differing number of rows to number of columns, and therefore is not a true dissimularity matrix")
    }
  }
  
  ##warns if distPreds are not matrices
  for(mat in distPreds){
    if(class(mat)!="matrix"){
      warning("One or more of the provided distance predictor matrices are not of class 'matrix'.")
    }
  }
  
  toRemove <- NULL
  removeRand <- NULL
  ##checks input data format
  ##species data as site-species matrix
  if(bioFormat==1){
    ##processes site-species matrix
    siteNum <- which(colnames(bioData)==siteColumn)
    
    ##point locations
    x <- which(colnames(bioData)==XColumn)
    ##checks to see if coordinates are found with bioData, if not extracts them from envData
    ##coordinates need to be in bioData if envData is raster
    if(is.na(x[1])==T){
      x <- which(colnames(predData)==XColumn)
      y <- which(colnames(predData)==YColumn)
      locs <- predData[c(x,y)]
      sppDat <- bioData[,-c(siteNum)]
    }else{
      y <- which(colnames(bioData)==YColumn)
      locs <- bioData[c(x,y)]
      sppDat <- bioData[,-c(siteNum, x, y)]
    }
    
    ##checks unique sites against rasters
    if(class(predData)=="RasterStack" | class(predData)=="RasterLayer" | class(predData)=="RasterBrick"){
      spSiteCol = bioData[siteColumn]
      siteRaster <- predData[[1]]
      ##raster based site
      x = which(colnames(bioData)==XColumn)
      y = which(colnames(bioData)==YColumn)
      locs = bioData[c(x,y)]
      cellID <- as.data.frame(cellFromXY(siteRaster, locs))
      colnames(cellID)[colnames(cellID)=="cellFromXY(siteRaster, locs)"] <- "cellName"
      
      if(nrow(cellID)==sum(is.na(cellID$cellName))){
        stop("None of the given points intersect with the given raster data. Double check that you geography is correct and that the given XColumn and YColumn values are correct.")
      }
      names(spSiteCol)[names(spSiteCol)==siteColumn] <- "siteIDn"
      
      ##checks size, and if needed correcting
      ##if cellID less than siteColumn
      if(length(unique(cellID$cellName))<length(unique(spSiteCol$siteIDn))){
        warning("One or more sites falls in the same raster cell. Using cells as sites instead")
      }else if(length(unique(cellID$cellName))>length(unique(spSiteCol$siteIDn))){
        ##if siteColumn less than cellID
        warning("More unique cells identified than unique sites in site column, defaulting sites to cells")
      }
      
      ##new to place extraction here
      rastEXDat <- as.data.frame(extract(predData, cellID$cellName))
      predData <- cbind(cellID, locs, rastEXDat)
      predData <- predData[order(predData$cellName),]
      
      ##aggregates bioData data into classes by raster cells
      newDataTable <- bioData[-c(siteNum, x, y)]
      modDataTable <- cbind(cellID, newDataTable)
      bioData <- cbind(cellID, bioData[-c(siteNum)])
      bioData <- bioData[order(bioData$cellName),]
      
      ##filters out sites with low species count
      ##prep for filtering
      headDat <- bioData[,c(siteNum)]
      sppDat[sppDat>=1]=1
      sppDat[sppDat==0]=0
      sppTotals <- cbind(as.data.frame(headDat), apply(sppDat, 1, function(m){sum(as.numeric(m))}))
      ##filters out data
      filterBioDat <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] >= sppFilter)
      spSiteCol <- filterBioDat["headDat"]
      colnames(spSiteCol) <- "cellName"
      ##reassembles data after filtering
      bioData <- merge(spSiteCol, bioData, by="cellName")
      ##identifies what to remove from distpred matrix data based on spp filter
      filterOut <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] < sppFilter)
      toRemove <- which(headDat %in% filterOut$headDat)
      
      ##removes random samples and sets up to 
      if(is.null(samples)==FALSE){
        if(nrow(bioData)<samples){
          warning("After species filter, fewer remaining records remaining than specified in samples, continuing without farther removals")
        }else{
          randRows <- sample(1:nrow(bioData), samples)
          fullLength <- 1:nrow(bioData)
          removeRand <- fullLength[-(randRows)] 
          ##actual selection of the random rows to keep
          bioData <- bioData[c(randRows),]
        } 
      }
      
      ##removes items by
      predData <- predData[which(bioData$cellName %in% predData$cellName),]
      colnames(bioData)[colnames(bioData)=="cellName"] <- siteColumn
      colnames(predData)[colnames(predData)=="cellName"] <- siteColumn
      #uniqueCells = unique(cellID$cellName)
      if(abundance==T){
        ##compress by cell
        byCell <- aggregate(modDataTable, modDataTable[1], FUN=mean)
        ##removes unneeded columns
        inDataTable <- byCell[-c(1, 2)]
      }else{
        ##compress by cell
        byCell <- aggregate(modDataTable, modDataTable[1], FUN=sum)
        ##transforms data to binary
        byCell <- byCell[-c(1, 2)]
        byCell[byCell>=1] <- 1
        inDataTable <- byCell
      }
      
    }else{
      ##Table data
      ##filters out sites with low species count
      ##prep for filtering
      headDat <- bioData[,c(siteNum)]
      sppDat[sppDat>=1]=1
      sppDat[sppDat==0]=0
      sppTotals <- cbind(as.data.frame(headDat), apply(sppDat, 1, function(m){sum(as.numeric(m))}))
      ##filters out data
      filterBioDat <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] >= sppFilter)
      spSiteCol <- filterBioDat["headDat"]
      colnames(spSiteCol) <- siteColumn
      ##reassembles data after filtering
      bioData <- merge(spSiteCol, bioData, by=siteColumn)
      ##identifies what to remove from distpred matrix data based on spp filter
      filterOut <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] < sppFilter)
      toRemove <- which(headDat %in% filterOut$headDat)
      
      ##removes random samples and sets up to 
      if(is.null(samples)==FALSE){
        if(nrow(bioData)<samples){
          warning("After species filter, fewer remaining records remaining than specified in samples, continuing without farther removals")
        }else{
          randRows <- sample(1:nrow(bioData), samples)
          fullLength <- 1:nrow(bioData)
          removeRand <- fullLength[-(randRows)] 
          ##actual selection of the random rows to keep
          bioData <- bioData[c(randRows),]
        } 
      }
      
      ##rename for comparisons
      colnames(bioData)[colnames(bioData)==siteColumn] <- "gettingCoolSiteColumn"
      colnames(predData)[colnames(predData)==siteColumn] <- "gettingCoolSiteColumn"
      ##extracts predData by only those in bioData
      predData <- unique(predData)
      predData <- predData[which(bioData$gettingCoolSiteColumn %in% predData$gettingCoolSiteColumn),]
      ##rename for results
      colnames(bioData)[colnames(bioData)=="gettingCoolSiteColumn"] <- siteColumn
      colnames(predData)[colnames(predData)=="gettingCoolSiteColumn"] <- siteColumn
      #predSite <- which(names(predData)==siteColumn)
      #predData <- predData[-predSite]
      bx <- which(names(bioData)==XColumn)
      by <- which(names(bioData)==YColumn)
      siteNum <- which(names(bioData)==siteColumn)
      inDataTable <- bioData[-c(siteNum, bx, by)]
    }
    ##creates distance matrix
    if(abundance==F){
      distData <- vegdist(inDataTable, dist, binary=T)
    }else{
      distData <- vegdist(inDataTable, dist, binary=F)
    }
    
    ########################################################################
    ##species data as x,y,species list
  }else if(bioFormat==2){
    ##insert data if not available
    if(is.null(abundColumn)){
      warning("No abundance column specified, assuming biological data are presence")
      bioData["reallysupercoolawesomedata"] <- 1
      abundColumn <- "reallysupercoolawesomedata"
    }
    if(is.null(sppColumn)){
      stop("Need to define sppColumn argument when bioFormat=2")
    }
    if(sppColumn %in% names(bioData)){} else{
      stop("Cannot find sppColumn in bioData - check name?")
    }
    ##point locations
    x <- which(colnames(bioData)==XColumn)
    y <- which(colnames(bioData)==YColumn)
    locs <- bioData[,c(x,y)]
    
    ##processes coords, species list
    ##if environmental data is raster, checks size of data to makes sure there is not problems with end table binding
    if(class(predData)=="RasterStack" | class(predData)=="RasterLayer" | class(predData)=="RasterBrick"){
      ##raster based site
      ##site raster
      siteRaster <- predData[[1]]
      cellID <- as.data.frame(cellFromXY(siteRaster, locs))
      colnames(cellID)[colnames(cellID)=="cellFromXY(siteRaster, locs)"] <- "cellName"
      if(nrow(cellID)==sum(is.na(cellID$cellName))){
        stop("None of the given points intersect with the given raster data. Double just that you geography is correct and that the given XColumn and YColumn values are correct.")
      }
      
      ##if siteColumn has been provided....... now should always be true
      if(!is.null(siteColumn)){
        ##species site
        spSiteCol <- bioData[siteColumn]
        names(spSiteCol)[names(spSiteCol)==siteColumn] <- "siteIDn"
        ##checks size, and if needed correcting
        ##if cellID less than siteColumn
        if(length(unique(cellID$cellName))<length(unique(spSiteCol$siteIDn))){
          warning("One or more sites falls in the same raster cell. Using cells as sites instead")
        }else if(length(unique(cellID$cellName))>length(unique(spSiteCol$siteIDn))){
          ##if siteColumn less than cellID
          warning("more unique cells identified than unique sites in site column, defaulting sites to cells")
        }
        ##removes unneeded columns
        siteNum <- which(colnames(bioData)==siteColumn)
        newDataTable <- bioData[-c(siteNum, x, y)]
        headDat <- bioData[,c(siteNum)]
      }else{
        ##removes unneeded columns
        newDataTable <- bioData[-c(x, y)]
      }
      
      ##new to place extraction here
      rastEXDat <- as.data.frame(extract(predData, cellID$cellName))
      predData <- unique(cbind(cellID, locs, rastEXDat))
      
      ##casts data into site-species matrix
      modDataTable <- cbind(cellID, newDataTable)
      names(modDataTable)[names(modDataTable)==sppColumn] <- "spcodeUltimateCoolness"
      castData <- dcast(modDataTable, cellName~spcodeUltimateCoolness, value.var=abundColumn)
      
      ##aggregate data by cellID, and filters data by low spp number
      if(abundance==T){
        ##compress by cell
        cellName <- which(names(castData)=="cellName")
        byCell <- aggregate(castData, castData[cellName], FUN=mean)
        byCell <- byCell[-2]
        
        ##filters species data
        sppDat <- byCell[-c(cellName)]
        headDat <- byCell[c(cellName)]
        sppDat[sppDat>=1] <- 1
        sppDat[is.na(sppDat)] <- 0
        sppTotals <- cbind(headDat, rowSums(sppDat))
        ##filters out data
        filterBioDat <- subset(sppTotals, sppTotals["rowSums(sppDat)"] >= sppFilter)
        spSiteCol <- filterBioDat[cellName]
        ##reassembles data after filtering
        cellXY <- xyFromCell(siteRaster, spSiteCol$cellName)
        spSiteCol <- cbind(spSiteCol, cellXY)
        bioData <- merge(spSiteCol, byCell, by=cellName)
        #XColumn <- "x"
        #YColumn <- "y"
        ##identifies what to remove from distpred matrix data based on spp filter
        filterOut <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] < sppFilter)
        toRemove <- which(headDat$cellName %in% filterOut$cellName)
        
        ##removes random samples and sets up to 
        if(is.null(samples)==FALSE){
          if(nrow(bioData)<samples){
            warning("After species filter, fewer remaining records remaining than specified in samples, continuing without farther removals")
          }else{
            randRows <- sample(1:nrow(bioData), samples)
            fullLength <- 1:nrow(bioData)
            removeRand <- fullLength[-(randRows)] 
            ##actual selection of the random rows to keep
            bioData <- bioData[c(randRows),]
          } 
        }
        
        ##removes unneeded columns
        inDataTable <- bioData[-c(1,2,3)]
      }else{
        ##compress by cell
        cellName <- which(names(castData)=="cellName")
        byCell <- aggregate(castData, castData[cellName], FUN=sum)
        ##remove redunancy
        byCell <- byCell[-2]
        
        ##filters species data
        sppDat <- byCell[-c(cellName)]
        headDat <- byCell[c(cellName)]
        sppDat[sppDat>=1] <- 1
        sppTotals <- cbind(headDat, rowSums(sppDat))
        ##filters out data
        filterBioDat <- subset(sppTotals, sppTotals["rowSums(sppDat)"] >= sppFilter)
        spSiteCol <- filterBioDat[cellName]
        ##reassembles data after filtering
        cellXY <- xyFromCell(siteRaster, spSiteCol$cellName)
        spSiteCol <- cbind(spSiteCol, cellXY)
        bioData <- merge(spSiteCol, byCell, by=cellName)
        #XColumn <- "x"
        #YColumn <- "y"
        ##identifies what to remove from distpred matrix data based on spp filter
        filterOut <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] < sppFilter)
        toRemove <- which(headDat$cellName %in% filterOut$cellName)
        
        ##removes random samples and sets up to 
        if(is.null(samples)==FALSE){
          if(nrow(bioData)<samples){
            warning("After species filter, fewer remaining records remaining than specified in samples, continuing without farther removals")
          }else{
            randRows <- sample(1:nrow(bioData), samples)
            fullLength <- 1:nrow(bioData)
            removeRand <- fullLength[-(randRows)] 
            ##actual selection of the random rows to keep
            bioData <- bioData[c(randRows),]
          } 
        }
        
        ##transforms data to binary
        transBio <- bioData[-c(1,2,3)]
        transBio[transBio>=1] <- 1
        inDataTable <- transBio
      }
      predData <- predData[which(bioData$cellName %in% predData$cellName),]
      bioData <- bioData[order(bioData$cellName),]
      predData <- predData[order(predData$cellName),]
      colnames(bioData)[colnames(bioData)=="cellName"] <- siteColumn
      colnames(predData)[colnames(predData)=="cellName"] <- siteColumn
    }else{
      ##if environmental data is a table
      ##remove coords if needed
      newDataTable <- bioData[-c(x, y)]
      
      ##renames required columns
      if(siteColumn %in% colnames(newDataTable)){} else{
        stop("Predictor data table is missing site column.")
      }      
      names(newDataTable)[names(newDataTable)==siteColumn] <- "siteUltimateCoolness"
      names(newDataTable)[names(newDataTable)==sppColumn] <- "spcodeUltimateCoolness"
      castData <- dcast(newDataTable, siteUltimateCoolness~spcodeUltimateCoolness, value.var=abundColumn)
      
      ##filters species data
      siteName <- which(names(castData)=="siteUltimateCoolness")
      sppDat <- castData[-c(siteName)]
      headDat <- castData[c(siteName)]
      sppDat[sppDat>=1] <- 1
      sppDat[is.na(sppDat)] <- 0
      sppTotals <- cbind(headDat, rowSums(sppDat))
      ##filters out data
      filterBioDat <- subset(sppTotals, sppTotals["rowSums(sppDat)"] >= sppFilter)
      spSiteCol <- filterBioDat[siteName]
      ##reassembles data after filtering
      bioSiteCol <- which(names(bioData)==siteColumn)
      holdData <- bioData[c(bioSiteCol,x,y)]
      holdData <- aggregate(holdData, holdData[1], FUN=mean)
      holdData <- holdData[-1]
      siteXY <- merge(spSiteCol, holdData, by.x="siteUltimateCoolness", by.y=siteColumn, all=F)
      bioData <- merge(siteXY, castData, by="siteUltimateCoolness")
      ##identifies what to remove from distpred matrix data based on spp filter
      filterOut <- subset(sppTotals, sppTotals[colnames(sppTotals)[2]] < sppFilter)
      toRemove <- which(headDat$siteUltimateCoolness %in% filterOut$siteUltimateCoolness)
      
      ##removes random samples and sets up to 
      if(is.null(samples)==FALSE){
        if(nrow(bioData)<samples){
          warning("After species filter, fewer remaining records remaining than specified in samples, continuing without farther removals")
        }else{
          randRows <- sample(1:nrow(bioData), samples)
          fullLength <- 1:nrow(bioData)
          removeRand <- fullLength[-(randRows)] 
          ##actual selection of the random rows to keep
          bioData <- bioData[c(randRows),]
        } 
      }
      
      ##rename for comparisons
      colnames(predData)[colnames(predData)==siteColumn] <- "siteUltimateCoolness"
      predData <- unique(merge(siteXY, predData, by="siteUltimateCoolness"))
      predData <- predData[order(predData$siteUltimateCoolness),]
      bioData <- bioData[order(bioData$siteUltimateCoolness),]
      
      ##extracts predData by only those in bioData
      predData <- predData[which(bioData$siteUltimateCoolness %in% predData$siteUltimateCoolness),]
      
      ##rename for results
      colnames(bioData)[colnames(bioData)=="siteUltimateCoolness"] <- siteColumn
      colnames(predData)[colnames(predData)=="siteUltimateCoolness"] <- siteColumn
      
      if(!XColumn %in% colnames(predData)){
        thy <- which(colnames(predData)==paste(XColumn, ".y", sep=""))
        the <- which(colnames(predData)==paste(YColumn, ".y", sep=""))
        predData <- predData[-c(thy, the)]
        
        XColumn <- paste(XColumn, ".x", sep="")
        YColumn <- paste(YColumn, ".x", sep="") 
      }
      
      inDataTable <- bioData[-c(1,2,3)]
    }
    
    ##creates distance matrix
    ##another line of filtering...... don't remember why it is here
    inDataTable[is.na(inDataTable)] <- 0
    #inDataTableB <- inDataTable[apply(inDataTable, 1, function(x){sum(x)>=sppFilter}),]
    
    if(abundance==F){
      distData <- vegdist(inDataTable, dist, binary=T)
    }else{
      distData <- vegdist(inDataTable, dist, binary=F)
    }
    ########################################################################
    
    ##species data as site-site distance matrix
  }else if(bioFormat==3){
    ##site-site distance already calculated
    #castData = bioData
    distData <- lower.tri(as.matrix(bioData), diag=FALSE)
    distData <- as.vector(bioData[distData])
    predData <- unique(predData)
    predData <- predData[order(predData[siteColumn]),]
    ########################################################################
    
    ##site pair table, already preped 
  }else if(bioFormat==4){
    ##site-pair distance value
    outTable <- bioData
    ########################################################################
    
  }else{
    ##return error, bioFormat argument out of bounds
    stop(paste("bioFormat argument of '", as.character(bioFormat), "' is not an accepted input value", sep=""))
  }
  ########################################################################
  
  ##With the dissim distance calculated, creates and fills the table in gdm format
  if(bioFormat!=4){
    ##creates base site-pair table
    outTable <- as.data.frame(createSitePair(dist=distData, spdata=bioData, envInfo=predData, dXCol=XColumn, dYCol=YColumn, siteCol=siteColumn, weightsType=weightType, custWeights=custWeightVect))
  }else{
    outTable <- bioData
  }
  
  ##first checks that the size of the dissimilarity matrices, if any were provided
  ##then pastes any dissimilarity matrices onto the created site-pair table
  if(length(distPreds)>0){
    baseMat <- distPreds[[1]]
    ##checks to size of each dissimilarity matrix, to make sure they are all the same
    lapply(distPreds, function(mat, mat1){
      if(length(mat1)!=length(mat)){
        stop("The dimensions of your predictions matrices are not the same size.")
      }
    }, mat1=baseMat)
    
    if(length(toRemove)>0){
      distPreds <- lapply(distPreds, function(mat, tR){mat[-c(tR), -c(tR)]}, tR=toRemove)
      baseMat <- distPreds[[1]]
    }
    if(length(removeRand)>0){
      distPreds <- lapply(distPreds, function(mat, tR){mat[-c(tR), -c(tR)]}, tR=removeRand)
      baseMat <- distPreds[[1]]
    }
    
    ##checks the size of the dissimilarity matrices against the size of distData
    baseMatDat <- lower.tri(as.matrix(baseMat),diag=FALSE)
    baseMatDat <- as.vector(baseMat[baseMatDat])  
    if(nrow(outTable)!=length(baseMatDat)){
      stop("The dimensions of the distance predictor matrices do not match the biological data")
    }
    
    for(num in 1:length(distPreds)){
      #num <- 2
      ##isolate matrix
      matrixDat <- lower.tri(as.matrix(distPreds[[num]], diag=FALSE))
      sweetFreakenPredMatrix <- as.vector(distPreds[[num]][matrixDat])
      ##add matrix to table
      if(ncol(outTable)>6){
        ##break table up to insert matrix data columns
        baseSitePair <- outTable[,1:6]
        otherSitePair <- outTable[,7:ncol(outTable)]
        otherNames <- colnames(otherSitePair)
        s1SitePair <- as.data.frame(otherSitePair[,1:(ncol(otherSitePair)/2)])
        colnames(s1SitePair) <- otherNames[1:(ncol(otherSitePair)/2)]
        s2SitePair <- as.data.frame(otherSitePair[,(ncol(otherSitePair)/2+1):ncol(otherSitePair)])
        colnames(s2SitePair) <- otherNames[(ncol(otherSitePair)/2+1):ncol(otherSitePair)]
        ##formats data from dissimilarity matrices 
        s1Zeros <- as.data.frame(rep(0,length(sweetFreakenPredMatrix)))
        colnames(s1Zeros) <- paste("s1.matrix_", num, sep="")
        s2Mat <- as.data.frame(sweetFreakenPredMatrix)
        colnames(s2Mat) <- paste("s2.matrix_", num, sep="")
        ##restructures the output table
        outTable <- cbind(baseSitePair, s1SitePair, s1Zeros, s2SitePair, s2Mat)
      }else{
        ##formats data from dissimilarity matrices 
        s1Zeros <- as.data.frame(rep(0,length(sweetFreakenPredMatrix)))
        colnames(s1Zeros) <- paste("s1.matrix_", num, sep="")
        s2Mat <- as.data.frame(sweetFreakenPredMatrix)
        colnames(s2Mat) <- paste("s2.matrix_", num, sep="")
        ##restructures the output table
        outTable <- cbind(outTable, s1Zeros, s2Mat)
      }
    }
  }
  
  ##return output table
  class(outTable) <- c("gdmData", "data.frame")
  return(outTable)
}
##########################################################################


##########################################################################
createSitePair <- function(dist, spdata, envInfo, dXCol, dYCol, siteCol, 
                           weightsType, custWeights){
  ##Used in formatGDMData to transform data from a site-site distance matrix into
  ##a site pair format
  ##Input Variables:
  ##dist = the distance object representing the pair-wise dissimilarity of the species data
  ##envInfo = environmental data
  ##spdata = needed as may have X,Y coords
  ##dXCol = x coordinate
  ##dYCol = y coordinate
  ##siteCol = site column, taken from either the species or environmental tables
  ##predMatrices = Prediction matrices to be added to the site-pair table
  ##weightsType = the method of determining the site-pair weights
  ##custWeights = custom wieghts, as a vector, if given
  ##Output Variables:
  ##gdmTableFill = the complete site-pair table
  ###########################
  #dist = distData
  #spdata = bioData
  #envInfo = predData
  #dXCol = XColumn
  #dYCol = YColumn
  #siteCol = siteColumn
  #weightsType = weightType
  #custWeights = custWeightVect
  ###########################
  ##required libraries
  #require(raster)
  
  ##Create gdm ready table
  weightsType <- as.character(weightsType)
  distance <- as.vector(dist)
  ##calculates richness total, the sums of the two most populus sites
  if(weightsType[1]=="richness"){
    sppOnly <- spdata[-c(1,2,3)]
    sppSums <- rowSums(sppOnly)
    sppSiteSums <- cbind(spdata[1], sppSums)
    orderedSums <- sppSiteSums[order(-sppSiteSums[,2]),]
    richTotal <- orderedSums[1,2]+orderedSums[2,2]
  }
  
  ##Builds index needed for output gdm table format
  xCoord.S1 <- yCoord.S1 <- xCoord.S2 <- yCoord.S2 <- NULL
  s1 <- s2 <- NULL
  
  if((siteCol %in% colnames(envInfo))==T){
    count <- seq(length(unique(envInfo[,siteCol]))-1,1,-1)
  }else{
    count <- seq(length(unique(envInfo[,"siteUltimateCoolness"]))-1,1,-1)
  }
  s1 <- unlist(sapply(seq(length(count),1), function(y){c(s1, rep((max(count)-y)+1, times=y))}))
  s2 <- unlist(sapply(seq(length(count),1), function(y){c(s2, (max(count)-y+2):(max(count)+1))}))
  
  if(length(s1)!=length(distance)){
    stop("The length of distance values are not the same as the expected number of rows of the site-pair table, unable to proceed.")
  }
  
  if(weightsType[1]=="equal"){
    weights <- rep(1, times=length(distance))
  }else if(weightsType[1]=="custom"){
    weights <- custWeights
  }else{
    weights <- (sppSiteSums[s1, "sppSums"] + sppSiteSums[s2, "sppSums"]) / richTotal
  }
  gdmTable <- cbind(distance, weights)
  
  ##from environmental or species table, copy coordinates for site-pair table
  if((dXCol %in% colnames(envInfo))==T){
    if((siteCol %in% colnames(envInfo))==T){
      checkTab <- table(envInfo[siteCol])
    }else{
      checkTab <- table(envInfo["siteUltimateCoolness"])
    }
    
    if(sum(checkTab>1)>0){
      stop("A site has two or more unique entries of data associated with it. Double check you data for incosistancies.")
    }
    xCoord.S1 <- envInfo[s1, dXCol]
    xCoord.S2 <- envInfo[s2, dXCol]
    yCoord.S1 <- envInfo[s1, dYCol]
    yCoord.S2 <- envInfo[s2, dYCol]  
  }else if((dXCol %in% colnames(spdata))==T){
    xCoord.S1 <- spdata[s1, dXCol]
    xCoord.S2 <- spdata[s2, dXCol]
    yCoord.S1 <- spdata[s1, dYCol]
    yCoord.S2 <- spdata[s2, dYCol]
  }else{
    stop("X,Y Coordinates not found with unique sites, unable to complete site-pair table")
  }
  
  ##sets up output table
  gdmForm <- cbind(gdmTable, xCoord.S1, yCoord.S1, xCoord.S2, yCoord.S2)
  xhold <- which(names(envInfo)==dXCol)
  yhold <- which(names(envInfo)==dYCol)
  sitehold <- which(names(envInfo)==siteCol)
  sitehold2 <- which(names(envInfo)=="siteUltimateCoolness")
  envInfo <- envInfo[-c(xhold, yhold, sitehold, sitehold2)]
  
  ##fills output table
  if(ncol(envInfo)>0){
    gdmTableFill <- cbind(gdmForm, envInfo[s1,1:ncol(envInfo)], envInfo[s2,1:ncol(envInfo)])
    names.s1 <- paste("s1.",names(envInfo[1:ncol(envInfo)]), sep="")
    names.s2 <- paste("s2.",names(envInfo[1:ncol(envInfo)]), sep="")
    colnames(gdmTableFill) <- c(colnames(gdmTableFill)[1:6], names.s1, names.s2)
  }else{
    gdmTableFill <- gdmForm
  }
  
  ##returns results
  return(gdmTableFill)
}
##########################################################################


##########################################################################
isplineExtract <- function (gdmModel){
  ##Extracts Ispline data from a gdm model
  ##Input Variable:
  ##gdmModel = a gdm model object
  ##Output Variable:
  ##outData = a list of two matrices, one for x and y spline coordinates (graph space)
  ###########################
  #gdmModel = gdmOb
  ###########################
  ##Collects or sets simple data
  options(warn.FPU = FALSE)
  PSAMPLE <- 200
  preddata <- rep(0, times = PSAMPLE)
  pn <- gdmModel$predictors
  nPreds <- length(pn)
  yDat <- xDat <- matrix(0,PSAMPLE,nPreds)
  colnames(yDat) <- colnames(xDat) <- pn
  pmin <- 1
  pmax <- PSAMPLE
  
  ##cycles through each prodictor and fills the spline matrices
  splineindex <- 1
  for (i in 1:nPreds){ 
    ######
    #i<-1
    ######
    numsplines <- gdmModel$splines[i]
    z <- .C("GetPredictorPlotData", 
            pdata = as.double(preddata), 
            as.integer(PSAMPLE), 
            as.double(gdmModel$coefficients[splineindex:(splineindex + numsplines - 1)]), 
            as.double(gdmModel$knots[splineindex:(splineindex + numsplines - 1)]), 
            as.integer(numsplines),
            PACKAGE = "gdm")
    yDat[,i] <- z$pdata
    pmin <- pmin + PSAMPLE
    pmax <- pmax + PSAMPLE
    xDat[,i] <-  seq(from=gdmModel$knots[[(i*3)-2]],to=gdmModel$knots[[(i*3)]], length=PSAMPLE)
    splineindex <- splineindex + numsplines
  }
  
  ##lists and returns matrices
  outData <- list(xDat,yDat)
  return(outData)
}
##########################################################################