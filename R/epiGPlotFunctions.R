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
    dataPath <- paste("inst/extdata/", data, sep="")
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
#' Plots the expression and quantile data in a scatterplot
#'
#' @param dataset An n by 4 data frame containing necessary information of epigenetic genes information, usually loaded by the package's given function; see parseEpigeneticData and loadEpigeneticData
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
normalized <- FALSE
plotEpigeneticEV <- function(layout,
                             normalized=FALSE,
                             color.legend=TRUE,
                             color.grey=FALSE,
                             color.arrow="red",
                             size.arrow=0.5,
                             labels.which=NULL,
                             sample.class=NULL,
                             sample.pred=NULL){

    # Validate layout parameter for plotting.
    if(missing(layout)){
        stop("Error: layout is missing for plotting!")
    }
    if(any(is.na(layout))){
        stop("Error: layout's data contains missing values!")
    }
    if(!is.data.frame(layout)){
        stop("Error: layout is not a valid data frame!")
    }

    # Validate specified sample.class parameter for plotting.
    lLab <- c("cell_line", "fractionation", "primary_cell", "timecourse", "tissue")
    if(!is.null(sample.class)){
        if(length(sample.class) > 0 && sample.class %in% lLab){
            lLab <- lLab[sample.class %in% lLab]
        } else {
            stop("Error: invalid specified sample.class!")
        }
    }

    # Initialize ggplot's aesthetic
    aesP <- ggplot2::aes(x=layout[,1], y=layout[,2], color=layout[,3])
    ePlot <- ggplot2::ggplot(data=layout[,1:2]) + ggplot2::geom_point(aesP)

    xLab <- "Quantile over all genes"
    yLab <- if(!normalized) "Expression Values" else "Normalized Expression Values"

    ePlot <- ePlot + labs(x=xLab, y=yLab)


    if(!color.legend){
        ePlot <- ePlot + theme(legend.position="none")
    } else {
        ePlot <- ePlot + scale_color_manual(values=unique(layout[,3]),labels=lLab) + labs(color="Sample Class")
    }
    if(color.grey){
        ePlot <- ePlot + scale_color_grey(labels=lLab) + theme_classic()
    }
    if(!is.null(labels.which)){
        # testWhich <- c("serous cystadenocarcinoma cell line:HTOA", "chronic lymphocytic leukemia (T-CLL) cell line:SKW-3")
        # Todo: replace parameters with function parameters
        ePlot <- ePlot + ggrepel::geom_text_repel(data=layout[labels.which,],
                                                  min.segment.length = 0,
                                                  mapping=aes(x=layout[labels.which,1],
                                                              y=layout[labels.which,2],
                                                              label=rownames(layout[labels.which,]),
                                                              segment.size=1,
                                                              segment.color="red",
                                                              ),
                                                  nudge_x = .15,
                                                  nudge_y = .5,
                                                  arrow= arrow(length = unit(0.015, "npc")),
                                                  max.iter = 1000)
        ePlot
    }
    return(ePlot)
}

# head(pScaleRange(tmpDat[,3]))
# head(tmpDat[,4])
# coords <- data.frame(x=scaledEV, y=tmpDat[,4])
# colnames(coords) <- c("Expression Value","Quantile over all genes")
# layout <- cbind(coords,labs(title="", x=xLab, y="Quantile over all Genes"))
# sample.class <- "cell_line"
layoutEpigeneticEV <- function(dataFrame,
                               normalized=FALSE,
                               sample.class="none",
                               predict.which="none"){
    refDf <- dataFrame
    sampleClass <- c("cell_line", "fractionation", "primary_cell", "timecourse", "tissue")
    pColours <- getGeneClassColour(refDf[,1])
    if(sample.class!="none" && sample.class %in% sampleClass){
        refDf <- refDf[refDf[,1]==sample.class,]
        pColours <- getGeneClassColour(refDf[,1])
    }

    if(normalized){
        scaledEV <- pScaleRange(refDf[,3])
        xLab <- "Normalized Expression Value"
        yDf <- scale(dataFrame[,4])
    } else {
        scaledEV <- refDf[,3]
        xLab <- "Expression Value"
        yDf <- refDf[,4]
    }
    layout <- data.frame(x=scaledEV, y=yDf)
    layout$Colours <- pColours
    # tmpSplit <- sapply(tmpDat[,2], function(x) strsplit(x,"CNhs")[[1]][1], USE.NAMES = FALSE)

    rownames(layout) <- sapply(refDf[,2], function(x) strsplit(x,".CNhs")[[1]][1], USE.NAMES = FALSE)
    colnames(layout) <- c(xLab, "Quantile over all genes", "Colour")
    return(layout)
}

legendEpigeneticEV <- function(){

}

# head(tmpDat[,1])
# unique(tmpDat[,1])
# layout[layout[,3]=="#e66101",3] <- "#1b9e77"
getGeneClassColour <- function(dFClassName,
                               colourList=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")){
    if(length(colourList) != 5){
        stop("Error: custom colour list must be of length 5!")
    }
    colours <- dFClassName
    for(gClass in unique(dFClassName)){
        switch (gClass,
            "cell_line" = i <- 1,
            "fractionation" = i <- 2,
            "primary_cell" = i <- 3,
            "timecourse" = i <- 4,
            "tissue" = i <- 5
        )
        colours[colours==gClass] <- colourList[i]
    }
    return(colours)
}

pScaleRange <- function(dataFrame){
    return((dataFrame - min(dataFrame))/(max(dataFrame)-min(dataFrame)))
}

td <- loadEpigeneticData("expValues_NO66_HUMAN.rda")
tdl <- layoutEpigeneticEV(td, sample.class="primary_cell")
(tdp <- plotEpigeneticEV(tdl,sample.class="primary_cell"))

tdm <- td[td[,1]=="primary_cell",]
colnames(tdm) <- c("c","s","y","x")
tlm <- lm(x~y, tdm)
tdp + geom_abline(intercept=tlm$coefficients[1], slope=tlm$coefficients[2], colour="#7570b3", size=1)
head(tdm)
