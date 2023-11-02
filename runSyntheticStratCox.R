library(Matrix)
library(Andromeda)
library(Cyclops)
library(testthat)

GpuDevice <- listGPUDevices()[1]

source("functions/synthetic_experiments.R")

N <- c(100000, 1000000)
P <- c(1000, 2000)
strataPool <- c(1, 5, 50, 500, 5000, 50000, 500000)

# Cox
for (j in 1:length(P)) {
for (i in 1:length(N)) {
	writeLines(paste0("----------P = ", P[j], ", N = ", N[i], "----------"))
	strata <- strataPool[strataPool < N[i]]
	total_times <- rep(0, 2*length(strata))
	gh_times <- rep(0, 2*length(strata))
	for (k in 1:length(strata)) {
		t <- runFixedL1(nRows = N[i],
				nCovars = P[j],
				nStrata = strata[k],
				modelType = if(strata[k] == 1) "cox" else "cox_time",
				gpuDevice = gpuDevice)
		total_times[c(k, k+length(strata))] <- t$total
		gh_times[c(k, k+length(strata))] <- t$gh
	}
	writeLines("Total time")
	print(total_times[1:length(strata)])
	print(total_times[(length(strata)+1):(2*length(strata))])
	writeLines("GH time")
	print(gh_times[1:length(strata)])
	print(gh_times[(length(strata)+1):(2*length(strata))])
}
}

