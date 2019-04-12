library("BiocManager", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library("minfi", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library("qsmooth", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library("wateRmelon", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")

# after batch effect removal


# quantile normalisation necessary?
library("quantro")
quantro(object = , groupFactor = , useMedianNormalized = TRUE)

if(type == 'functional'){
  preprocessFunnorm(nPCs = 5, sex = , bgCorr = FALSE, dyeCorr = FALSE)
}else if (type == 'SmoothQN'){
  qsmooth(object = , groupFactor = )
}

