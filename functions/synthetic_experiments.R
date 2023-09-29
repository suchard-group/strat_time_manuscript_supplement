library(Matrix)
library(Andromeda)
library(Cyclops)
library(testthat)

simCox <- function(nRows,
		   nCovars = 1000,
		   beta,
		   nStrata = 1,
		   sparseness = 0.95,
		   zeroEffectSizeProp = 0.8,
		   ind = TRUE,
		   useAndr = TRUE,
		   seed = 123) {

    set.seed(seed)
 
    effectSizes <- data.frame(covariateId=1:nCovars,rr=exp(beta))

    # generate covariate matrix
    covariates <- rsparsematrix(nRows, nCovars, nnz = nRows*nCovars*(1-sparseness), rand.x = rnorm)
    covariates <- summary(covariates)

    if (ind) {
        covariates$x <- 1 # indicator
    }
    colnames(covariates) <- c("rowId", "covariateId", "covariateValue")

    # calculate sum of exb for each observation
    outcomes <- data.frame(rowId = 1:nRows, stratumId = round(runif(nRows,min=1,max=nStrata)), y=0)
    covariates <- merge(covariates,outcomes[,c("rowId","stratumId")])
    rowId_to_rr <- aggregate(rr ~ rowId, data=merge(covariates,effectSizes), prod)
    outcomes <- merge(outcomes,rowId_to_rr,all.x=TRUE)
    outcomes$rr[is.na(outcomes$rr)] <- 1

    # survival
    strataBackgroundProb <- runif(nStrata,min=0.01,max=0.03)
    outcomes$rate <-  strataBackgroundProb[outcomes$stratumId] * outcomes$rr
    outcomes$timeToOutcome <- 1+round(rexp(n=nrow(outcomes),outcomes$rate))
    outcomes$timeToCensor <- 1+round(runif(n=nrow(outcomes),min=0,max=499))
    outcomes$time <- outcomes$timeToOutcome
    outcomes$time[outcomes$timeToCensor < outcomes$timeToOutcome] <- outcomes$timeToCensor[outcomes$timeToCensor < outcomes$timeToOutcome]
    outcomes$y <- as.integer(outcomes$timeToCensor > outcomes$timeToOutcome)

    # wrap up
    outcomes <- outcomes[order(outcomes$stratumId,outcomes$rowId),]
    covariates <- covariates[order(covariates$stratumId,covariates$rowId,covariates$covariateId),]

    # break ties
    n <- dim(outcomes)[1]
    outcomes$time <- outcomes$time + rnorm(n, mean = 0, sd = 1E-3)

    sparsenessReal <- 1-(nrow(covariates)/(nRows*nCovars))
    writeLines(paste("Simulated sparseness = ",sparsenessReal*100,"%"))

    if (useAndr) {
        andr <- andromeda(out = outcomes, cov = covariates)
        return(andr)
    } else {
        return(list(out = outcomes, cov = covariates))
    }
}

runFixedL1 <- function(nRows,
		       nCovars = 1000,
		       nStrata = 1,
		       sparseness = 0.95,
		       zeroEffectSizeProp = 0.8,
		       effectSizeSd = 1,
		       ind = TRUE,
		       doublePrecision = TRUE,
		       variance = NULL,
		       seed = 123,
		       modelType = "cox_time",
		       runGPU = TRUE,
		       runCPU = TRUE,
		       gpuDevice) {
    
    set.seed(seed)

    ######################################
    # Simulate Data
    ######################################

    # initialize beta
    sd <- rep(effectSizeSd, nCovars) * rbinom(nCovars, 1, 1 - zeroEffectSizeProp)
    beta <- rnorm(nCovars,mean=0,sd=sd)

    andr <- simCox(nRows = nRows, nCovars = nCovars, beta = beta, nStrata = nStrata, sparseness = sparseness, zeroEffectSizeProp = zeroEffectSizeProp, ind = ind, seed = seed)
    
    ######################################
    # Run Experiments
    ######################################

    if (doublePrecision) {
        sd = 0.0001
        fp = 64
        tolerance <- 1E-4
    } else {
        sd = 0.1
        fp = 32
        tolerance <- 1E-3
    }

    if (is.null(variance)){
        prior = createPrior("none")
    } else {
        prior = createPrior("laplace", variance = variance)
    }

    if (runGPU) {
        # sparse gpu
        sparse_gpu_ptr <- convertToCyclopsData(andr$out, andr$cov, modelType = modelType, floatingPoint = fp)
        start <- Sys.time()
        sparse_gpu_fit <- fitCyclopsModel(sparse_gpu_ptr,
				          prior = prior,
					  computeDevice = gpuDevice)
        delta_g <- Sys.time() - start
        writeLines(paste("GPU sparse analysis took", signif(delta_g,3), attr(delta_g,"units"),
             "(",
             signif(sparse_gpu_fit$timeFit,3), attr(sparse_gpu_fit$timeFit,"units"),
             ")"))
    }

    if (runCPU) {
        # sparse cpu
        sparse_cpu_ptr <- convertToCyclopsData(andr$out, andr$cov, modelType = modelType, floatingPoint = fp)
        start <- Sys.time()
        sparse_cpu_fit <- fitCyclopsModel(sparse_cpu_ptr,
					  prior = prior)
        delta_c <- Sys.time() - start
        writeLines(paste("CPU sparse analysis took", signif(delta_c,3), attr(delta_c,"units"),
             "(",
             signif(sparse_cpu_fit$timeFit,3), attr(sparse_cpu_fit$timeFit,"units"),
             ")"))
    }

    # close Andromeda
    close(andr)

    ######################################
    # Check results
    ######################################

    # compare GPU v.s. CPU
    if (runGPU && runCPU) {
        expect_equal(coef(sparse_gpu_fit), coef(sparse_cpu_fit), tolerance = tolerance)
        expect_equal(logLik(sparse_gpu_fit), logLik(sparse_cpu_fit), tolerance = tolerance)
    }

    # compare with true beta
    if (runGPU && runCPU) {
        mse <- mean((beta - coef(sparse_gpu_fit))^2)
    } else if (runGPU) {
        mse <- mean((beta - coef(sparse_gpu_fit))^2)
    } else if (runCPU) {
        mse <- mean((beta - coef(sparse_cpu_fit))^2)
    }
    writeLines(paste("MSE = ", mse))

    # return time
    if (runGPU && runCPU) {
        return(list(total = c(delta_g, delta_c),
		    gh =    c(sparse_gpu_fit$timeFit, sparse_cpu_fit$timeFit)))
    } else if (runGPU) {
        return(list(total = delta_g, gh = sparse_gpu_fit$timeFit))
    } else if (runCPU) {
        return(list(total = delta_c, gh = sparse_cpu_fit$timeFit))
    }
}


