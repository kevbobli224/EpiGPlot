# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


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
parseEpigeneticData <- function(data="expressions.csv"){
    dataPath <- paste("../inst/extdata/", data, sep="")
    if(!file.exists(dataPath)){
        stop(sprintf("Error: data not found in: %s", dataPath))
    }
    dF <- read.delim(dataPath, sep="\t")
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
#' @import ggplot2
#'
#' @export
loadEpigeneticData <- function(rda){
    if(!nchar(rda)){
        stop("Error: Empty rda file name provided!")
    }
    dataPath <- paste("../data/", rda, sep="")
    if(!file.exists(dataPath)){
        stop("Error: Missing rda file for expression values!")
    }
    pEnv <- new.env()
    tmpDat <- get(load("../data/expValues.rda", envir=pEnv))
    return(tmpDat)
}

#' Plots the the expression values of a epigenetic gene.
#'
#' Plots the expression and quantile data in a scatterplot
#'
#' @param dataset An n by 4 data frame containing necessary information of epigenetic genes information, usually loaded by the package's given function; see parseEpigeneticData and loadEpigeneticData
#'
#' @import ggplot2
#' @import ggforce
plotEpigeneticEV <- function(dataset=expValues,
                             color.color= 'black'){
    update_geom_defaults(ggforce::GeomCircle, list(colours=color.color))
    aesP <- ggplot2::aes(x=tmpDat[,3], y=tmpDat[,4])
    ggplot(tmpDat[,3:4], mapping=aesP) + geom_point()
    ggplot2::geom_point(expValues[,3:4], mapping=aes(x=expValues[,3]))
}

# head(pScaleRange(tmpDat[,3]))
# head(tmpDat[,4])
layoutEpigeneticEV <- function(dataFrame,
                               normalized=FALSE){
    if(normalized){
        scaledEV <- pScaleRange(dataFrame[,3])
        xLab <- "Normalized Expression Value"
    } else {
        scaledEV <- dataFrame[,3]
        xLab <- "Expression Value"
    }

    layout <- cbind(labs(title="", x=xLab, y="Quantile over all Genes"))


}

legendEpigeneticEV <- function(){

}

# head(tmpDat[,1])
# unique(tmpDat[,1])
getGeneClassColour <- function(dataFrame,
                               colourList=c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")){

}

pScaleRange <- function(dataFrame){
    return((dataFrame - min(dataFrame))/(max(dataFrame)-min(dataFrame)))
}


