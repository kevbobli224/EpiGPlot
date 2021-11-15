library(EpiGPlot)

testthat::test_that("Check parsing functions and data set",{
    # Test for data frame's dimensions from the provided data set
    expect_equal(dim(EpiGPlot::NO66_HUMAN), c(889, 4))

    # Set working directory to package's directory
    setwd("../../../EpiGPlot")

    # Test loading rda data for data in data/ subdirectory
    testData <- EpiGPlot::loadEpigeneticData("NO66_HUMAN.rda")
    # Test for dimensions
    expect_equal(dim(testData), c(889, 4))
    # Incorrect input types for loadEpigeneticData
    expect_error(EpiGPlot::loadEpigeneticData(""), "Error: Empty rda file name provided!")
    expect_error(EpiGPlot::loadEpigeneticData("dne.rda"), "Error: Missing rda file for expression values!")

    # Test parsing csv data for data in inst/extdata/ subdirectory
    testData <- EpiGPlot::parseEpigeneticData("expressions.csv")
    # Test for dimensions
    expect_equal(dim(testData), c(889, 4))
    # Incorrect input types for parseEpigeneticData
    expect_error(EpiGPlot::parseEpigeneticData(""), "Error: Empty csv file name provided!")
    expect_error(EpiGPlot::parseEpigeneticData("dne.csv"), "Error: data not found in: inst/extdata/dne.csv")
})

testthat::test_that("Checking layout function", {
    testData <- EpiGPlot::NO66_HUMAN
    layout <- EpiGPlot::layoutEpigeneticEV(testData)
    # Test layout data frame's dimensions
    expect_equal(dim(layout), c(889,4))
    # Test for properly assigned column names
    expect_equal(colnames(layout), c("Expression Value", "Quantile over all genes", "Colour", "Sample Class"))

    # Incorrect input types for layoutEpigeneticEV
    expect_error(EpiGPlot::layoutEpigeneticEV(testData,
                                              sample.class=c("tissue"),
                                              class.colour = c("#ffffff", "#aaaaaa")),
                 "Error: mismatch class.colour size compared to sample.class")

    expect_error(EpiGPlot::layoutEpigeneticEV(testData, class.colour = c("red", "black", "blue", "yellow", "not a colour")),
                 "Error: class.colour contains invalid colour element.")

    expect_error(EpiGPlot::layoutEpigeneticEV(testData, class.colour = c("red")),
                 "Error: length of class.colour must be of length 5 for unspecified sample classes")

    # Test different parameters of layoutEpigeneticEV
    layout <- EpiGPlot::layoutEpigeneticEV(testData, sample.class=c("tissue"), normalized = TRUE)
    expect_equal(colnames(layout), c("Normalized Expression Value", "Quantile over all genes", "Colour", "Sample Class"))
    expect_equal(dim(layout), c(135,4))
    expect_equal(unique(layout[,4]), c("tissue"))
})
