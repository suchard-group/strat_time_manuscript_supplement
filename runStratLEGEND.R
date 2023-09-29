library(Andromeda)
library(Cyclops)

gpuDevice <- listOpenCLDevices()[1]

source("functions/LEGEND.R")

date <- "0521"
numberStratas <- c(5, 10, 20)
variances <- c(0.00970277, 0.00921064, 0.0106282)
pathToData <- "/home/jianxiao/LEGEND/Data/"
treatmentVarId <- readRDS(paste0(pathToData,"treatmentVarId.rds"))

for (i in 1:length(numberStratas)) {
    # Load data
    covariateData <- loadAndromeda(paste0(pathToData,"cmd_strata", numberStratas[i], ".zip"))

    # Run L1 Cox (PS stratified)
    fit <- runLEGEND(covariateData = covariateData,
		     treatmentVarId = treatmentVarId,
		     gpuDevice = gpuDevice,
		     useCrossValidation = FALSE,
		     variance = variances[i])

    # Save results
    saveRDS(fit, paste0(date, "_fitLEGENDStrata", numberStratas[i], ".rds"))
}
