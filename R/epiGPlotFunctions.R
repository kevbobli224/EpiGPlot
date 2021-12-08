#' Parses epigenetic gene's expression values from .csv files
#'
#' .csv file must contain 4 columns, separated by tabs. The column must contain following data in a particular order: sample class, sample name,  expression value, and quantile range. .csv files regarding epigenetic file data can be obtained through EpiFactors website.
#'
#' @param data A file name for the epigenetic gene expression values file including the .csv suffix
#'
#' @return Returns a data frame containing 4 columns of: Sample class, sample name, expression value, and quantile range.
#'
#' @importFrom utils read.delim
#'
#' @examples
#' \dontrun{
#' data <- parseEpigeneticData("./inst/extdata/expressions.csv")
#' dim(data) # 889 4
#' }
#'
#' @references
#' Medvedeva, Y. A., Lennartsson, A., Ehsani, R., Kulakovskiy, I. V., Vorontsov, I. E., Panahandeh, P., Khimulya, G., Kasukawa, T., &amp; Drabløs, F. (2015). Epifactors: A comprehensive database of human epigenetic factors and complexes. Database, 2015.
#' \href{http://database.oxfordjournals.org/content/2015/bav067.full}{Link}
#'
#' @export
parseEpigeneticData <- function(data){
    # Check if the data's path provided is empty
    if(!nchar(data)){
        stop("Error: Empty csv file name provided!")
    }
    # Check if the data exists in the given data's path
    if(!file.exists(data)){
        stop(sprintf("Error: data not found in: %s", data))
    }

    # Read file into variable dF
    dF <- utils::read.delim(data, sep="\t")

    # Check if number of columns provided in the .csv is not equal to 4
    if(ncol(dF) != 4){
        stop("Error: invalid epigenetic expression file, ensure the columns selected for exportation is equal to 4.")
    }
    return(dF)
}

#' Loads epigenetic gene's data from an existing .rda file
#'
#' @param rda A string or char array containing the file including the .rda suffix
#'
#' @return Returns a data frame containing 4 columns of: Sample class, sample name, expression value, and quantile range.
#'
#' @examples
#' \dontrun{
#' data <- loadEpigeneticData("./data/NO66_HUMAN.rda")
#' dim(data) # 889 4
#' }
#'
#' @export
loadEpigeneticData <- function(rda){
    # Check if the data's path provided is empty
    if(!nchar(rda)){
        stop("Error: Empty rda file name provided!")
    }
    # Check if the data exists in the given data's path
    if(!file.exists(rda)){
        stop("Error: Missing rda file for expression values!")
    }
    # Set new private environment such that data set can be assigned when loaded
    pEnv <- new.env()
    tmpDat <- get(load(rda, envir=pEnv))
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
#' @examples
#' \dontrun{
#' rawData <- EpiGPlot::loadEpigeneticData("data/NO66_HUMAN.rda")
#' gLayout <- EpiGPlot::layoutEpigeneticEV(rawData,
#'                                         normalized=TRUE,
#'                                         sample.class=c("cell_line", "tissue"))
#' gPlot <- EpiGPlot::plotEpigeneticEV(gLayout,
#'                                     normalized=TRUE,
#'                                     sample.pred = c("tissue"),
#'                                     labels.top=5)
#' gPlot
#'
#' gPlot <-
#' EpiGPlot::plotEpigeneticEV(gLayout, sample.pred = c("tissue", "cell_line"),
#'                          labels.which=c("chronic lymphocytic leukemia (T-CLL) cell line:SKW-3"))
#'
#' gPlot
#' }
#'
#' @export
plotEpigeneticEV <- function(layout,
                             normalized=FALSE,
                             colour.legend=TRUE,
                             colour.grey=FALSE,
                             size.arrow=1,
                             size.pred=1,
                             labels.which=NULL,
                             labels.top=0,
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
        if(length(sample.class) > 0 && any(sample.class %in% lLab)){
            layout <- layout[layout[,4] %in% sample.class,]
            lLab <- lLab[lLab %in% sample.class]
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
    ePlot <- ePlot + ggplot2::labs(x=xLab, y=yLab)

    # User-specified legend parameters
    if(!colour.legend){
        ePlot <- ePlot + ggplot2::theme(legend.position="none")
    } else {
        ePlot <- ePlot + ggplot2::scale_color_manual(values=cMat[1,],labels=lLab) + ggplot2::labs(colour="Sample Class")
    }
    if(colour.grey){
        ePlot <- ePlot + ggplot2::scale_color_grey(labels=lLab) + ggplot2::theme_classic()
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
                                                  mapping=ggplot2::aes(x=layout[parsedLabels,1],
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
#' @examples
#' \dontrun{
#' rawData <- EpiGPlot::loadEpigeneticData("data/NO66_HUMAN.rda")
#' gLayout <- EpiGPlot::layoutEpigeneticEV(rawData,
#'                                        sample.class=c("cell_line", "tissue"),
#'                                        class.colour=c("red", "blue"))
#' nrow(gLayout) # 390
#' unique(gLayout$Colour) # "red", "blue"
#' }
#' @export
layoutEpigeneticEV <- function(dataFrame,
                               normalized=FALSE,
                               sample.class=c(),
                               class.colour=c()){
    # Initialize a reference data frame that will be modified to accommodate
    # for layout specification
    refDf <- dataFrame

    # Check if user has specify a group of class to be displayed
    if(length(sample.class)>0){
        # Subset the input data based on specified list of class
        refDf <- refDf[which(refDf[,1] %in% sample.class),]
        # Check if user has specify a group of colours to be displayed
        if(length(class.colour)==0){
            # Get a list of colours based on the associated sample class
            pColours <- getGeneClassColour(refDf[,1])
        } else{
            # Check if length of sample.class and class.colour matches for
            # equal assignments
            if(length(sample.class) != length(class.colour)){
                stop("Error: mismatch class.colour size compared to sample.class")
            }
            # Check if class.colour contains invalid colours.
            if("try-error" %in% class(try(grDevices::col2rgb(class.colour), silent=TRUE))){
                stop("Error: class.colour contains invalid colour element.")
            }
            # Get a list of colours based on the associated sample class with
            # specified colours from class.colour
            pColours <- getGeneClassColour(refDf[,1], class.colour)
        }
    } else {
        # Check if specified class.colour is less than the number of unique sample classes
        if(length(class.colour)<length(unique(refDf[,1])) && length(class.colour) > 0){
            stop("Error: length of specified class.colour is lesser than the number of sample classes for unspecified sample classes")
        } else if(length(class.colour) == 0){
            # For unspecfied class.colour, assign pColours with default colour palette
            pColours <- getGeneClassColour(refDf[,1])
        } else {
            # Check if class.colour contains invalid colours.
            if("try-error" %in% class(try(grDevices::col2rgb(class.colour), silent=TRUE))){
                stop("Error: class.colour contains invalid colour element.")
            }
            # Get a list of colours based on the associated sample class with
            # specified colours from class.colour
            pColours <- getGeneClassColour(refDf[,1], class.colour)
        }
    }

    if(normalized){
        scaledEV <- pScaleRange(refDf[,3])
        xLab <- "Normalized Expression Value"
        yDf <- scale(refDf[,4])
    } else {
        scaledEV <- refDf[,3]
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
    # Check if unique sample classes exceed the number of colours
    if(length(unique(dFClassName)) > length(colourList)){
        stop("Error: unique class names exceed number of available colours.")
    }

    # Initialize a vector of classes assigned with its respective colour
    availClass <- unique(dFClassName)
    names(availClass) <- colourList[1:length(availClass)]

    # Assign colours to each sample classes
    colours <- dFClassName
    for(gClass in availClass){
        colours[colours==gClass] <- names(availClass)[availClass == gClass]
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

#' Obtains a 2 x n colour matrix for cross-referencing colour values and sample classes
#'
#' @param dataFrame a layout data frame, usually used in the package's plotting function
#'
#' @return Returns a 2 x n matrix, where the first row contains hex colour values, and the second row contains its corresponding sample class assigned colour
getColourMatrix <- function(dataFrame){
    # Obtain list of unique colours
    uCol <- unique(dataFrame[,3])
    # Create a 1 x uCol matrix
    cMat <- matrix(uCol, nrow=1, ncol=length(uCol))
    # Associate sample class to colour
    cMat <- rbind(cMat, unique(dataFrame[,4]))
    return(cMat)
}
