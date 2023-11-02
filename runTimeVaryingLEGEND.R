library(Andromeda)
library(Cyclops)

gpuDevice <- listGPUDevices()[1]

source("functions/LEGEND.R")

date <- "0921"
safetyOutcome <- "Cough"
cut <- c(10) # cutoff for time-varying effect
variance <- NULL # 0.00146317
bootstrapFileName <- paste0(date, "_bootstrap_", safetyOutcome, ".txt")

# Load data
pathToData <- "/home/jianxiao/LEGEND/SafetyOutcomes/"
treatmentVarId <- readRDS(paste0(pathToData,"treatmentVarId", safetyOutcome, ".rds"))
covariateData <- loadAndromeda(paste0(pathToData,"cmd", safetyOutcome, ".zip"))

# Run L1 Cox with time-varying coefficients
fit <- runLEGEND(covariateData = covariateData,
		 treatmentVarId = treatmentVarId,
		 cut = cut,
		 gpuDevice = gpuDevice,
		 useCrossValidation = TRUE,
		 variance = variance,
		 runBootstrap = TRUE,
		 bootstrapFileName = bootstrapFileName)

# Save results
saveRDS(fit, paste0(date, "_results_", safetyOutcome, ".rds"))

