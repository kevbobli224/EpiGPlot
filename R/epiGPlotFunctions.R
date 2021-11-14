#' Parses epigenetic gene's expression values from .csv files
#'
#' Requires data to be in inst/extdata folder, and file should contain 4 columns, separated by tabs.
#'
#' @param data A file name for the epigenetic gene expression values file including the .csv suffix
#'
#' @return Returns a data frame containing 4 columns of: Sample class, sample name, normalized expression value, and quantile range.
#'
#' @importFrom utils read.delim
#'
#' @export
parseEpigeneticData <- function(data){
    if(!nchar(data)){
        stop("Error: Empty rda file name provided!")
    }
    dataPath <- paste("inst/extdata/", data, sep="")
    if(!file.exists(dataPath)){
        stop(sprintf("Error: data not found in: %s", dataPath))
    }
    dF <- utils::read.delim(dataPath, sep="\t")
    if(ncol(dF) < 4){
        stop("Error: invalid epigenetic genes expression value file format.")
    }
    return(dF)
}

#' Loads epigenetic gene's data from .rda file from the data subdirectory
#'
#' @param rda A string or char array containing the file name in the data subdirectory including the .rda suffix
#'
#' @return Returns a data frame containing 4 columns of: Sample class, sample name, normalized expression value, and quantile range.
#'
#' @export
loadEpigeneticData <- function(rda){
    if(!nchar(rda)){
        stop("Error: Empty rda file name provided!")
    }
    dataPath <- paste("data/", rda, sep="")
    if(!file.exists(dataPath)){
        stop("Error: Missing rda file for expression values!")
    }
    pEnv <- new.env()
    tmpDat <- get(load(dataPath, envir=pEnv))
    return(tmpDat)
}

#' Plots the the expression values of a epigenetic gene.
#'
#' Plots the expression vs quantile data in a scatterplot, linear regression can be visualized for user-specified sample class amongst other user-specified parameters
#'
#' @param layout An n by 4 data frame containing necessary information of epigenetic genes information, usually loaded by the package's given function; see parseEpigeneticData and loadEpigeneticData
#' @param normalized Whether the dataset "layout" expression value is normalized.
#' @param colour.legend Whether legend on the plot should be dispalyed
#' @param colour.grey Whether if greyscale should be applied to all visual details of the plot
#' @param size.arrow The size of the line segment if user chooses to display top expression values
#' @param size.pred The size of linear regression model's line segment if user chooses to display predictions on sample class
#' @param labels.which User can specify which sample to be displayed on the plot through a list/vector
#' @param labels.top User can specify how many highest expression value samples to be displayed
#' @param labels.size The font size of sample names
#' @param sample.class User can specify which sample class to display
#' @param sample.pred User can specify which sample class to perform linear regression on, and display their respective line segments
#'
#' @return Returns a plot that user can assign from its return value.
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#'
#' @export
plotEpigeneticEV <- function(layout,
                             normalized=FALSE,
                             colour.legend=TRUE,
                             colour.grey=FALSE,
                             size.arrow=1,
                             size.pred=1,
                             labels.which=NULL,
                             labels.top=3,
                             labels.size=3,
                             sample.class=NULL,
                             sample.pred=NULL){

    # Validate layout parameter for plotting.
    if(missing(layout)){
        stop("Error: layout is missing for plotting!")
    }

    # Validate if there are any NA values in given layout
    if(any(is.na(layout))){
        stop("Error: layout's data contains missing values!")
    }

    # Validate if layout is a data frame
    if(!is.data.frame(layout)){
        stop("Error: layout is not a valid data frame!")
    }

    # Initialize lLab for legend labels for sample class
    lLab <- c("cell_line", "fractionation", "primary_cell", "timecourse", "tissue")

    # Validate specified sample.class parameter for plotting.
    if(!is.null(sample.class)){
        if(length(sample.class) > 0 && sample.class %in% lLab){
            lLab <- lLab[sample.class %in% lLab]
        } else {
            stop("Error: invalid specified sample.class!")
        }
    }
    # Initialize colour matrix for sample class colour reference
    cMat <- getColourMatrix(layout)

    # Initialize ggplot's aesthetic
    aesP <- ggplot2::aes(x=layout[,1], y=layout[,2], colour=layout[,3])

    # Plot using geom_point or scatterplot
    ePlot <- ggplot2::ggplot() + ggplot2::geom_point(aesP)

    # Initialize x and y axis labels
    xLab <- if(!normalized) "Expression Values" else "Normalized Expression Values"
    yLab <- "Quantile over all genes"

    # Modify plot labels
    ePlot <- ePlot + labs(x=xLab, y=yLab)

    # User-specified legend parameters
    if(!colour.legend){
        ePlot <- ePlot + theme(legend.position="none")
    } else {
        ePlot <- ePlot + scale_color_manual(values=cMat[1,],labels=lLab) + labs(colour="Sample Class")
    }
    if(colour.grey){
        ePlot <- ePlot + scale_color_grey(labels=lLab) + theme_classic()
    }

    # Labelling of user-specified sample genes, if given
    if(!is.null(labels.which) || labels.top > 0){
        parsedLabels <- c()
        if(!is.null(labels.which)){
            parsedLabels <- sapply(labels.which, function(x) strsplit(x,".CNhs")[[1]][1], USE.NAMES = FALSE)
            if(!any(parsedLabels %in% rownames(layout))){
                stop("Error: specified label(s) does not exists in given layout data")
            }
            parsedLabels <- parsedLabels[parsedLabels %in% rownames(layout)]
        }

        if(labels.top > 0){
            sortedEV <- sort(layout[,1], decreasing = TRUE)[1:labels.top]
            parsedLabels <- c(parsedLabels, rownames(layout[which(layout[,1] %in% sortedEV),]))
        }
        ePlot <- ePlot + ggrepel::geom_text_repel(data=layout[parsedLabels,1:2],
                                                  min.segment.length = 0,
                                                  size=labels.size,
                                                  segment.size=size.arrow,
                                                  mapping=aes(x=layout[parsedLabels,1],
                                                              y=layout[parsedLabels,2],
                                                              label=parsedLabels,
                                                              colour=layout[parsedLabels,3]),
                                                  max.iter = 10000,
                                                  box.padding = .5,
                                                  show.legend = FALSE
                                                  )
    }

    # Performs linear regression on user-specified sample class, if given
    if(!is.null(sample.pred)){
        if(length(sample.pred) > 0 && any(sample.pred %in% lLab)){
            parsedPred <- lLab[lLab %in% sample.pred]
            oColNames <- colnames(layout)
            colnames(layout) <- c("x","y","c","S")
            for(pred in parsedPred){
                modelData <- layout[which(layout[,4] %in% pred),]
                mdl <- stats::lm(y~x,modelData)$coefficients
                ePlot <- ePlot + ggplot2::geom_abline(intercept = mdl[1],
                                             slope = mdl[2],
                                             colour=cMat[1,cMat[2,] %in% pred],
                                             size=size.pred)
            }
            colnames(layout) <- oColNames
        } else {
            stop("Error: invalid specified sample.pred!")
        }
    }

    # Returns the plot
    return(ePlot)
}
#' Create a layout dataset format for plotting and better data representation
#'
#' Cleans up the parsed csv/rda file data from obtained dataset, and returns a layout that can be used for plotting
#'
#' @param dataFrame The n x 4 data frame loaded from csv/rda using the functions above
#' @param normalized Whether or not to perform normalization on expression values
#' @param sample.class Chooses and subset a specific sample class from the data frame and use it for layout
#' @param class.colour User can specify a list/vector custom colour palette of length equals to number of given sample classes
#'
#' @return Returns a layout containing 4 columns of: expression value, quantile over all genes, colour, sample class with their row names assigned as their sample genes
#'
#' @importFrom grDevices col2rgb
#'
#' @export
layoutEpigeneticEV <- function(dataFrame,
                               normalized=FALSE,
                               sample.class=NA,
                               class.colour=NA){
    refDf <- dataFrame
    sampleClass <- c("cell_line", "fractionation", "primary_cell", "timecourse", "tissue")
    if(!is.na(sample.class) && any(sample.class %in% sampleClass)){
        refDf <- refDf[refDf[,1]==sample.class,]
        if(is.na(class.colour)){
            pColours <- getGeneClassColour(refDf[,1])
        } else{
            if(length(sample.class) != length(class.colour)){
                stop("Error: mismatch class.colour size compared to sample.class")
            }
            if("try-error" %in% class(try(grDevices::col2rgb(class.colour), silent=TRUE))){
                stop("Error: class.colour contains invalid colour element.")
            }
            pColours <- getGeneClassColour(refDf[,1], class.colour)
        }
    } else {
        if(length(class.colour)!=5 && !is.na(class.colour)){
            stop("Error: length of class.colour must be of length 5 for unspecified sample classes")
        } else if(is.na(class.colour)){
            pColours <- getGeneClassColour(refDf[,1])
        } else {
            pColours <- getGeneClassColour(refDf[,1], colourList = class.colour)
        }
    }

    if(normalized){
        scaledEV <- refDf[,3]
        xLab <- "Normalized Expression Value"
        yDf <- scale(dataFrame[,4])
    } else {
        scaledEV <- pScaleRange(refDf[,3])
        xLab <- "Expression Value"
        yDf <- refDf[,4]
    }
    layout <- data.frame(x=scaledEV, y=yDf)
    layout$Colours <- pColours

    rownames(layout) <- sapply(refDf[,2], function(x) strsplit(x,".CNhs")[[1]][1], USE.NAMES = FALSE)
    layout <- cbind(layout, refDf[,1])
    colnames(layout) <- c(xLab, "Quantile over all genes", "Colour", "Sample Class")
    return(layout)
}

#' Obtains a colour scheme from its corresponding sample class for plotting purposes
#'
#' @param dFClassName A n x 1 dataframe, providing a list of sample classes
#' @param colourList A custom colour palette for plotting purposes
#'
#' @return A n x 1 data frame with mapped colour string to its corresponding sample class
getGeneClassColour <- function(dFClassName,
                               colourList=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")){
    sColourList <- sort(colourList)
    colours <- dFClassName
    for(gClass in unique(dFClassName)){
        switch (gClass,
            "cell_line" = i <- 1,
            "fractionation" = i <- 2,
            "primary_cell" = i <- 3,
            "timecourse" = i <- 4,
            "tissue" = i <- 5
        )
        colours[colours==gClass] <- sColourList[i]
    }
    return(colours)
}

#' Obtains a scale range from 0 to 1 for normalized value plotting purposes
#'
#' @param dFExp a single column of data frame that provides expression values of sample genes
#'
#' @return Returns the vectorized operated scaled values over the entire column
pScaleRange <- function(dFExp){
    return((dFExp - min(dFExp))/(max(dFExp)-min(dFExp)))
}

#' Obtains a 2x5 colour matrix for cross-referencing colour values and sample classes
#'
#' @param dataFrame a layout data frame, usually used in the package's plotting function
#'
#' @return Returns a 2x5 matrix, where the first row contains hex colour values, and the second row contains its corresponding sample class assigned colour
getColourMatrix <- function(dataFrame){
    cMat <- matrix(unique(dataFrame[,3]), nrow=1, ncol=5)
    cMat <- rbind(cMat, unique(dataFrame[,4]))
    return(cMat)
}

# td <- loadEpigeneticData("expValues_NO66_HUMAN.rda")
# nc <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3")
# tdl <- layoutEpigeneticEV(td)
# (tdp <- plotEpigeneticEV(tdl,sample.class="primary_cell",labels.top = 7, sample.pred=tp))
#
# tdm <- td[td[,1]=="primary_cell",]
# colnames(tdm) <- c("c","s","y","x")
# tlm <- lm(x~y, tdm)
# tdp + geom_abline(intercept=tlm$coefficients[1], slope=tlm$coefficients[2], colour="#7570b3", size=1)
# head(tdm)
