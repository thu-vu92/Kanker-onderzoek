# https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

# keep only primary tumor
pheno_data <- read.table("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/gene_meta_combat.csv",
                         sep="\t", header=TRUE)
pheno_data <- as.data.frame(pheno_data)
pheno_data <- pheno_data[!duplicated(pheno_data[[1]]),]
pheno_data[, 'Array.name'] <- gsub("-", ".", as.character(pheno_data[, 'Array.name']))
rownames(pheno_data) <- pheno_data[,1]

pheno_data$sample <- seq.int(nrow(pheno_data))
pheno_data$Sample.name<-NULL
pheno_data$sample<-NULL

diag1 = as.character(pheno_data[pheno_data['cat_diagnosis']=='Lung Squamous Cell Carcinoma', 'Array.name'])
diag2 = as.character(pheno_data[pheno_data['cat_diagnosis']=='Lung Adenocarcinoma', 'Array.name'])

pheno1 = pheno_data[pheno_data$Array.name %in% diag1, ]
pheno2 = pheno_data[pheno_data$Array.name %in% diag2, ]

array_order = as.character(pheno_data[, 'Array.name'])
array_order1 = as.character(pheno1[, 'Array.name'])
array_order2 = as.character(pheno2[, 'Array.name'])

batch = as.factor(pheno_data$Batch)
batch1= as.factor(pheno1$Batch)
batch2= as.factor(pheno2$Batch)

#diag1 <- gsub("-", ".", diag1)
#diag2 <- gsub("-", ".", diag2)
####
modcombat<-NULL
# +as.numeric(as.factor(cat_diagnosis))+as.factor(cat_gender)+as.numeric(as.factor(cat_source_site))

# impute age with median age..
for(i in 1:ncol(pheno_data)){
  pheno_data[is.na(pheno_data[,i]), i] <- mean(pheno_data[,i], na.rm = TRUE)
}
modcombat = model.matrix(~as.factor(cat_gender), data=pheno_data) # +as.numeric(cat_age)+as.numeric(cat_packs), cat_smoke_years
modcombat1 = model.matrix(~as.factor(cat_gender), data=pheno1)
modcombat2 = model.matrix(~as.factor(cat_gender), data=pheno2)

mod0 = model.matrix(~1, data=pheno_data)
mod01 = model.matrix(~1, data=pheno1)
mod02 = model.matrix(~1, data=pheno2)

full_data <- read.table("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/gene_raw.csv",
                        sep="\t", 
                        header=TRUE)

full_data <- as.data.frame(full_data)
full_data <- full_data[!duplicated(full_data[[1]]),]
rownames(full_data) <- full_data[,1]
full_data[,1]<- NULL 

full_rows = rownames(full_data)

full_data1 <- full_data[, diag1]
full_data2 <- full_data[, diag2] 

sample_batch_df= pheno_data[c("Array.name", "Batch")]
names(sample_batch_df) <- c('sample_id', 'batch_id')

#####
BEclear::correctBatchEffect(data = full_data1, 
                            samples = sample_batch_df, 
                            adjusted = 0.01,
                            method = TRUE, 
                            colBlockSize = 0, 
                            rowBlockSize = 0,
                            epochs = 5,
                            dir = paste(root_folder, "_prepped/RNAex_BEclear_LSCC.csv", sep=""),
                            outputFormat = 'txt')

BEclear::correctBatchEffect(data = full_data2, 
                            samples = sample_batch_df, 
                            adjusted = 0.01,
                            method = TRUE, 
                            colBlockSize = 0, 
                            rowBlockSize = 0,
                            epochs = 5,
                            dir = paste(root_folder, "_prepped/RNAex_BEclear_LA.csv", sep=""),
                            outputFormat = 'txt')
