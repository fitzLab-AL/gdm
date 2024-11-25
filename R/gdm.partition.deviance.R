#' @title Perform Deviance Partitioning of a Fitted GDM
#'
#' @description Partitions deviance explained from GDM into different
#' user specified components - most typically environment versus space.
#'
#' @usage gdm.partition.deviance(sitePairTable, varSets=list(), partSpace=TRUE)
#'
#' @param sitePairTable A correctly formatted site-pair table from
#' \code{\link[gdm]{formatsitepair}}.
#'
#' @param varSets A list in which each element is a vector of variable names
#' across which deviance partitioning is to be performed, excluding
#' geographic distance (which is set by the partSpace argument). Variable names
#' must match those used to build the site-pair table. See example.
#'
#' @param partSpace Whether or not to perform the partitioning using
#' geographic space. Default=TRUE.
#'
#' @return A dataframe summarizing deviance partitioning results.
#'
#' @author Matt Fitzpatrick and Karel Mokany
#'
#' @examples
#' # set up site-pair table using the southwest data set
#' sppData <- southwest[, c(1,2,13,14)]
#' envTab <- southwest[, c(2:ncol(southwest))]
#' sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat",
#' sppColumn="species", siteColumn="site", predData=envTab)
#'
#' # EXAMPLE - Partition two groups of variables
#' # Make list of variable sets for partitioning
#' varSet <- vector("list", 2)
#'
#' # now, name the variable groups for partitioning
#' # note you do not need to add "space" as this is only needed
#' # for environmental variables
#'
#' # two groups (soils & climate)
#' names(varSet) <- c("soil", "climate")
#'
#' # lastly, add variable names for
#' varSet$soil <- c("awcA", "phTotal", "sandA", "shcA", "solumDepth")
#' varSet$climate <- c("bio5", "bio6", "bio15", "bio18", "bio19")
#' varSet
#'
#' # run the function to partition soils, climate, and space (partSpace=TRUE)
#' scgPart <- gdm.partition.deviance(sitePairTab, varSet, partSpace=TRUE)
#'
#' # EXAMPLE - Partition three groups of variables
#' # Make list of variable sets for partitioning
#' varSet <- vector("list", 3)
#' names(varSet) <- c("soil", "temp", "precip")
#' varSet$soil <- c("awcA", "phTotal", "sandA", "shcA", "solumDepth")
#' varSet$temp <- c("bio5", "bio6")
#' varSet$precip <- c("bio15", "bio18", "bio19")
#'
#' # partition soils, temperature, and precip
#' # note we can't also partition space given the function's limit to a
#' # maximum of three variable sets, so we set partSpace=FALSE
#' scPart <- gdm.partition.deviance(sitePairTab, varSet, partSpace=FALSE)
#'
#' @export
gdm.partition.deviance <- function(sitePairTable, varSets=list(), partSpace=TRUE){

  #sitePairTable <- sitePairTab
  #varSets <- varSet
  #partSpace <- F

  if(is(varSets, "list")==FALSE){
    stop("Variables names must be given in a list.")
  }

  if(length(varSets)>=3 & partSpace==T){
    stop("Cannot partition more than three variables sets.")
  }

  if(length(varSets)>3){
    stop("Cannot partition more than three variables sets.")
  }

  if(!"gdmData" %in% class(sitePairTable)){
    stop("Site-pair table must be class=gdmData.")
  }

  # number of variables sets to partition across,
  # excluding geo
  partLength <- length(varSets)

  # make a third 'total' variable set by combining
  # the two provided (assuming two are provided, otherwise, skip)
  if(partLength==2){
    varSets[[length(varSets)+1]] <- as.character(unlist(varSets))
    names(varSets)[partLength+1] <- paste(names(varSets)[1:partLength],
                                          collapse = " & ")}

  if(partLength==3){
    varSets[[length(varSets)+1]] <- as.character(unlist(varSets[1:2]))
    names(varSets)[partLength+1] <- paste(names(varSets)[1:2],
                                          collapse = " & ")
    varSets[[length(varSets)+1]] <- as.character(unlist(varSets[c(1,3)]))
    names(varSets)[partLength+2] <- paste(names(varSets)[c(1,3)],
                                          collapse = " & ")
    varSets[[length(varSets)+1]] <- as.character(unlist(varSets[2:3]))
    names(varSets)[partLength+3] <- paste(names(varSets)[2:3],
                                          collapse = " & ")
    varSets[[length(varSets)+1]] <- as.character(unlist(varSets[1:3]))
    names(varSets)[partLength+4] <- paste(names(varSets[1:3]),
                                          collapse = " & ")
    }

  # if partitioning using geo, setup to fit models w/ w/o geo
  if(partSpace){
    partGeo <- c(FALSE, TRUE)
  } else {partGeo <- FALSE}

  # objects to hold results
  dev <- NULL
  variableSets <- NULL

  # Work through each partition
  # Start with geo=T, as needed
  # then select variable sets from site-pair table & fit GDM to each
  if(partSpace){
    #print("geo = ")
    dev <- gdm(sitePairTable[,1:6], geo=T)$explained
    variableSets <- ("geo")
    # loop through all combinations of geo & other variable sets
    for(g in partGeo){
      for(v in 1:length(varSets)){
        #print(paste("geo = ", g, " and ", names(varSets)[v], sep=""))
        # select variable sets
        greppers <- c(paste("s1.", varSets[[v]], sep=""),
                      paste("s2.", varSets[[v]], sep=""))
        modSPT <- sitePairTable[, c(1:6,which(names(sitePairTable) %in% greppers))]
        # fit GDM, get deviance
        dev <- c(dev, gdm(modSPT, geo=g)$explained)
        if(g){
          ggg <- "geo & "
        } else {
          ggg <- ""
        }
        variableSets <- c(variableSets, paste(ggg, names(varSets)[v], sep=""))
      }
    }

    # data frame with results
    devTable <- data.frame(VARIABLE_SET = variableSets, DEVIANCE = dev,
                           stringsAsFactors = F)
    devTable[7,1] <- paste("ALL VARIABLES (", devTable[7,1], ")", sep="")

    # partition deviance
    if(length(varSets)==1){
      devTable[4,1] <- paste(names(varSets), "alone", sep=" ")
      devTable[4,2] <- devTable[3,2]-devTable[1,2]
      devTable[5,1] <- "geo alone"
      devTable[5,2] <- devTable[3,2]-devTable[2,2]
      devTable[6,1] <- "UNEXPLAINED"
      devTable[6,2] <- 100-devTable[3,2]
    }

    if(length(varSets)==3){
      devTable <- devTable[c(2,3,1,4:7),]
      devTable[9,1] <- paste(names(varSets)[1], "alone", sep=" ")
      devTable[9,2] <- devTable[7,2]-devTable[6,2]
      devTable[10,1] <- paste(names(varSets)[2], "alone", sep=" ")
      devTable[10,2] <- devTable[7,2]-devTable[5,2]
      devTable[11,1] <- "geo alone"
      devTable[11,2] <- devTable[7,2]-devTable[4,2]

      #Var1
      A <- devTable[1,2]
      #Var2
      B <- devTable[2,2]
      #Var3
      C <- devTable[3,2]

      # ALL VARIABLES
      A_B_C <- devTable[7,2]

      #var1 alone
      a <- devTable[9,2]
      #var2 alone
      b <- devTable[10,2]
      #var3 alone
      c <- devTable[11,2]

      eqs <- matrix(c(1,1,0,1,
                      1,0,1,1,
                      0,1,1,1,
                      1,1,1,1),4,4,byrow=TRUE)
      colnames(eqs) <- letters[4:7]
      diffs <- c(A-a, B-b, C-c, A_B_C-a-b-c)
      parts <- solve(eqs, diffs, vars=colnames(eqs))

      devTable[12,1] <- paste(names(varSets)[1], " intersect ", names(varSets)[2],
                              ", exclude geo", sep="")
      devTable[12,2] <- parts[1]
      devTable[13,1] <- paste(names(varSets)[1], " intersect geo",
                              ", exclude ", names(varSets)[2], sep="")
      devTable[13,2] <- parts[2]
      devTable[14,1] <- paste("geo intersect ", names(varSets)[2],
                              ", exclude ", names(varSets)[1], sep="")
      devTable[14,2] <- parts[3]
      devTable[15,1] <- paste(names(varSets)[2], " intersect ", names(varSets)[1],
                              " intersect geo", sep="")
      devTable[15,2] <- parts[4]
      devTable[8,1] <- "UNEXPLAINED"
      devTable[8,2] <- 100-devTable[7,2]
    }
  }

  # same as above, but excluding geo if requested
  if(partSpace==F & length(varSets)<=3){
    for(g in partGeo){
      for(v in 1:length(varSets)){
        #print(paste("geo = ", g, " and ", names(varSets)[v], sep=""))
        greppers <- c(paste("s1.", varSets[[v]], sep=""),
                      paste("s2.", varSets[[v]], sep=""))
        modSPT <- sitePairTable[, c(1:6,which(names(sitePairTable) %in% greppers))]

        dev <- c(dev, gdm(modSPT, geo=g)$explained)

        if(g){
          ggg <- "geo and "
        } else {
          ggg <- ""
        }
        variableSets <- c(variableSets, paste(ggg, names(varSets)[v], sep=""))
      }
    }

    devTable <- data.frame(variableSet = variableSets, deviance = dev,
                           stringsAsFactors = F)
    devTable[3,1] <- paste("ALL VARIABLES (", devTable[3,1], ")", sep="")

    devTable[5,1] <- paste(names(varSets)[1], "alone", sep=" ")
    devTable[5,2] <- devTable[3,2]-devTable[2,2]
    devTable[6,1] <- paste(names(varSets)[2], "alone", sep=" ")
    devTable[6,2] <- devTable[3,2]-devTable[1,2]
    devTable[7,1] <- paste(names(varSets)[1], " intersect ", names(varSets)[2],sep="")
    devTable[7,2] <- devTable[3,2]-devTable[5,2]-devTable[6,2]
    devTable[4,1] <- "UNEXPLAINED"
    devTable[4,2] <- 100-devTable[3,2]

    # vennPlot <- draw.pairwise.venn(area1=round(devTable[1,2], 2),
    #                                    area2=round(devTable[2,2], 2),
    #                                    cross.area = round(devTable[6,2], 2),
    #                                    category=c(names(varSets)[1],
    #                                               names(varSets)[2]),
    #                                    euler.d=T, cat.col="darkred",
    #                                scaled=T)
    # grid.newpage()
    #     grid.draw(vennPlot)

  }

  if(length(varSets)==7){
    ggg <- ""
    for(v in 1:length(varSets)){
      #print(paste("geo = ", g, " and ", names(varSets)[v], sep=""))
      greppers <- c(paste("s1.", varSets[[v]], sep=""),
                    paste("s2.", varSets[[v]], sep=""))
      modSPT <- sitePairTable[, c(1:6,which(names(sitePairTable) %in% greppers))]

      dev <- c(dev, gdm(modSPT, geo=F)$explained)

      variableSets <- c(variableSets, paste(ggg, names(varSets)[v], sep=""))
    }

    # data frame with results
    devTable <- data.frame(VARIABLE_SET = variableSets, DEVIANCE = dev,
                           stringsAsFactors = F)
    devTable[7,1] <- paste("ALL VARIABLES (", devTable[7,1], ")", sep="")

    # individual sets
    devTable[9,1] <- paste(names(varSets)[1], "alone", sep=" ")
    devTable[9,2] <- devTable[7,2]-devTable[6,2]
    devTable[10,1] <- paste(names(varSets)[2], "alone", sep=" ")
    devTable[10,2] <- devTable[7,2]-devTable[5,2]
    devTable[11,1] <- paste(names(varSets)[3], "alone", sep=" ")
    devTable[11,2] <- devTable[7,2]-devTable[4,2]

    #Var1
    A <- devTable[1,2]
    #Var2
    B <- devTable[2,2]
    #Var3
    C <- devTable[3,2]

    # ALL VARIABLES
    A_B_C <- devTable[7,2]

    #var1 alone
    a <- devTable[9,2]
    #var2 alone
    b <- devTable[10,2]
    #var3 alone
    c <- devTable[11,2]

    eqs <- matrix(c(1,1,0,1,
                  1,0,1,1,
                  0,1,1,1,
                  1,1,1,1),4,4,byrow=TRUE)
    colnames(eqs) <- letters[4:7]
    diffs <- c(A-a, B-b, C-c, A_B_C-a-b-c)
    parts <- solve(eqs, diffs, vars=colnames(eqs))

    # intersections
    devTable[12,1] <- paste(names(varSets)[1], " intersect ", names(varSets)[2],
                            ", exclude ", names(varSets)[3], sep="")
    devTable[12,2] <- parts[1]
    devTable[13,1] <- paste(names(varSets)[1], " intersect ", names(varSets)[3],
                            ", exclude ", names(varSets)[2], sep="")
    devTable[13,2] <- parts[2]
    devTable[14,1] <- paste(names(varSets)[2], " intersect ", names(varSets)[3],
                            ", exclude ", names(varSets)[1], sep="")
    devTable[14,2] <- parts[3]
    devTable[15,1] <- paste(names(varSets)[1], " intersect ", names(varSets)[2],
                            " intersect ", names(varSets)[3], sep="")
    devTable[15,2] <- parts[4]
    devTable[8,1] <- "UNEXPLAINED"
    devTable[8,2] <- 100-devTable[7,2]

  }
  devTable[,2] <- round(devTable[,2], 2)
  devTable[,2] <- ifelse(devTable[,2]<0, 0, devTable[,2])
  return(devTable)
}
