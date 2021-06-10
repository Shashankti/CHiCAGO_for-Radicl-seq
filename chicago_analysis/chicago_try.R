library(Chicago)
library(PCHiCdata)
vignette("Chicago")

testDesignDir <- file.path("Data/Processed/DesFile/")
dir(testDesignDir)
testDataPath <- file.path("Data/Processed/InputFile/")
dir(testDataPath)
files <- file.path(testDataPath, "13.chinput")
settingsFile <- file.path(system.file("extdata", package="PCHiCdata"),
                          "sGM12878Settings", "sGM12878.settingsFile")
cd <- setExperiment(designDir = testDesignDir, settingsFile = settingsFile)
cd@settings$otherEndIDcol= "otherID"
cd <- readAndMerge(files=files, cd=cd)
cd <- normaliseBaits(cd)
cd <- normaliseOtherEnds(cd)
cd <- estimateTechnicalNoise(cd)
cd <- estimateDistFun(cd)
cd@x$s_j <- !is.infinite(cd@x$s_j)
cd <- estimateBrownianComponent(cd)
cd <- getPvals(cd)
cd <- getScores(cd)

outputDirectory <- testDataPath
exportResults(cd,file.path(outputDirectory,"vignetteOutput"))
plottedBaitIDs <- plotBaits(cd,n=1)
#######
cd@settings$adjBait2bait = FALSE
###

#Example dataset
dataPath <- system.file("extdata", package="PCHiCdata")
testDesignDir <- file.path(dataPath, "hg19TestDesign")
dir(testDesignDir)
testDataPath <- file.path(dataPath, "GMchinputFiles")
dir(testDataPath)
files <- c(
  file.path(testDataPath, "GM_rep1.chinput"),
  file.path(testDataPath, "GM_rep2.chinput"),
  file.path(testDataPath, "GM_rep3.chinput")
)
settingsFile <- file.path(system.file("extdata", package="PCHiCdata"),
                          "sGM12878Settings", "sGM12878.settingsFile")
example <- setExperiment(designDir = testDesignDir, settingsFile = settingsFile)
example <- readAndMerge(files=files,cd=example)
example <- chicagoPipeline(example)
exportResults(example,file.path(outputDirectory,"exampleOutput"))

