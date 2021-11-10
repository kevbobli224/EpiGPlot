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


#' Parses epigenetic genes expression values from .csv files
#'
#' Requires data to be in data folder, and file should contain 4 columns, separated by tabs.
#'
#' @param data A file name for the epigenetic gene expression values file
#'
#' @return A data frame containing 4 columns of: Sample class, sample name, normalized expression value, and quantile range.
#'
#' @importFrom utils read.delim
#'
#' @export
parseEpigeneticData <- function(data="expressions.csv"){
    dataPath <- paste("../inst/extdata", data, sep="")
    if(!file.exists(dataPath)){
        stop(sprintf("Error: data not found in: %s", dataPath))
    }
    dF <- read.delim(dataPath, sep="\t")
    if(ncol(dF) < 4){
        stop("Error: invalid epigenetic genes expression value file format.")
    }
    return(dF)
}

loadEpigeneticData <- function(){
    if(!file.exists("../data/expValues.rda")){
        stop("Error: Missing rda file for expression values!")
    }
    tmpEnv <- new.env()
    tmpDat <- get(load("../data/expValues.rda", envir=tmpEnv))
    rm(tmpEnv)
    return(tmpDat)
}
