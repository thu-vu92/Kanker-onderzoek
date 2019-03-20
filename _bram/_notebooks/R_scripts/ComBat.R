# https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

# keep only primary tumor
pheno_data <- read.table("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_meta_small.csv",
                         sep="\t", header=TRUE)
pheno_data <- as.data.frame(pheno_data)
pheno_data <- pheno_data[!duplicated(pheno_data[[1]]),]
rownames(pheno_data) <- pheno_data[,1]

pheno_data$sample <- seq.int(nrow(pheno_data))
pheno_data$Sample.name<-NULL
pheno_data$sample<-NULL

batch=as.factor(pheno_data$Batch)

####
modcombat = model.matrix(~as.factor(cat_packs)+as.factor(cat_age)
                         +as.factor(as.character(cat_gender))
                         +as.factor(cat_source_site)+as.factor(cat_sample_type), data=pheno_data)
mod0 = model.matrix(~1, data=pheno_data)

full_data <- read.table("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_raw.csv",
                        sep="\t", header=TRUE)

full_data <- as.data.frame(full_data)
full_data <- full_data[!duplicated(full_data[[1]]),]
rownames(full_data) <- full_data[,1]
full_data[,1]<- NULL

full_matrix = as.matrix(full_data)
#dimnames(full_matrix) <- list(rownames(full_data), colnames(full_data))

####
library("sva", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
source("/home/koekiemonster/DEV/GIT/RexR/_hackathon2018/_notebooks/CComBat.R")
corrected_batch_only = CComBat(full_matrix, batch = batch, mod = NULL, par.prior = FALSE, prior.plots = TRUE)
corrected = CComBat(full_matrix, batch = batch, mod = modcombat)
