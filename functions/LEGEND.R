library(Andromeda)
library(Cyclops)
library(testthat)

runLEGEND <- function(covariateData,
		      treatmentVarId,
		      doublePrecision = TRUE,
		      seed = 123,
		      modelType = "cox_time",
		      runGPU = TRUE,
		      runCPU = TRUE,
		      gpuDevice,
		      threads = 1,
		      useCrossValidation = TRUE,
		      variance = NULL,
		      cut = NULL,
		      runBootstrap = FALSE,
		      bootstrapFileName = "bootstrapOut.txt",
		      bootstrapReplicates = 100) {

    if (!runGPU && !runCPU) stop("Nothing to run")
    if (runGPU && is.null(gpuDevice)) stop("Need to specify GPU device")

    set.seed(seed)
    
    if (doublePrecision) {
        tolerance <- 1E-4
    } else {
        tolerance <- 1E-3
    }

    if (is.null(cut)) {
        # Exclude treatment covariate from regularization
        exclude <- treatmentVarId
    } else {
        # Exclude treatment covariates from regularization
        exclude <- c(treatmentVarId:(treatmentVarId+length(cut)))

        # Specify time-varying effect
        timeEffectMap <- data.frame(covariateId = c(treatmentVarId))

        # Split outcome table
        longOut <- splitTime(shortOut = covariateData$outcomes, cut = cut)

        # Long data in Andromeda
        covariateData <- andromeda(outcomes = longOut,
				   covariates = covariateData$covariates %>%
					   mutate(stratumId = as.integer(1)))
    }

    ######################################
    # Run Experiments
    ######################################

    # Create prior and control for cross-validation
    if (useCrossValidation) { # regularization with cv
        prior <- createPrior(priorType = "laplace",
			     variance = 1,
			     exclude = exclude,
			     useCrossValidation = useCrossValidation)
        control <- createControl(cvType = "auto",
				 seed = seed,
				 threads = threads,
				 startingVariance = 0.01,
				 tolerance = 2e-07,
				 cvRepetitions = 1,
				 noiseLevel = "quiet")
    } else if (!is.null(variance)) { # regularization with a given variance
        prior <- createPrior(priorType = "laplace",
			     variance = variance,
			     exclude = exclude)
        control <- createControl()
    } else { # no regularization
        prior <- createPrior(priorType = "none")
        control <- createControl()
    }

    # Fit model
    if (runGPU) {
        cyclopsDataG <- convertToCyclopsData(outcomes = covariateData$outcomes,
                                            covariates = covariateData$covariates,
                                            timeEffectMap = timeEffectMap,
                                            modelType = modelType,
                                            checkRowIds = FALSE,
                                            quiet = TRUE)
        start <- Sys.time()
        fitGPU <- fitCyclopsModel(cyclopsDataG,
                               prior = prior,
                               control = control,
                               computeDevice = gpuDevice)
        delta_g <- Sys.time() - start
        writeLines(paste("GPU sparse analysis took", signif(delta_g,3), attr(delta_g,"units"),
             "(",
             signif(fitGPU$timeFit,3), attr(fitGPU$timeFit,"units"),
             ")"))
             
        if (runBootstrap) {
            bs <- runBootstrap(fitGPU, bootstrapFileName, as.character(treatmentVarId), bootstrapReplicates)
	}
    }

    if (runCPU) {
        cyclopsDataC <- convertToCyclopsData(outcomes = covariateData$outcomes,
                                            covariates = covariateData$covariates,
                                            timeEffectMap = timeEffectMap,
                                            modelType = modelType,
                                            checkRowIds = FALSE,
                                            quiet = TRUE)
        start <- Sys.time()
        fitCPU <- fitCyclopsModel(cyclopsDataC,
                               prior = prior,
                               control = control)
        delta_c <- Sys.time() - start
        writeLines(paste("CPU sparse analysis took", signif(delta_c,3), attr(delta_c,"units"),
             "(",
             signif(fitCPU$timeFit,3), attr(fitCPU$timeFit,"units"),
             ")")) 
    }


    ######################################
    # Check results
    ######################################

    # compare GPU v.s. CPU
    if (runGPU && runCPU) {
        expect_equal(coef(fitGPU), coef(fitCPU), tolerance = tolerance)
        expect_equal(logLik(fitGPU), logLik(fitCPU), tolerance = tolerance)
	writeLines(paste0("Speedup of total time:                ", as.numeric(delta_c)/as.numeric(delta_g)))
	writeLines(paste0("Speedup of target of parallelization: ", as.numeric(fitCPU$timeFit)/as.numeric(fitGPU$timeFit)))
    }

    # return results
    if (runGPU) {
        fitGPU
        return(list(fit = fitGPU))
    } else {
        fitCPU
        return(list(fit = fitCPU))
    }
}

