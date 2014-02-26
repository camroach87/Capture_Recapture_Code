############################### PREAMBLE ################################
# PROGRAM: main
# CODER: Cameron Roach
# DATE: 24/2/2014
# DESCRIPTION: Sets directories and calls programs


# Set directories
dataDir <- file.path("./data")
srcDir <- file.path("./src")
outputDir <- file.path("./output")


# Make capture matrix
source(file.path(srcDir,"makeCaptureMatrix.R"))
