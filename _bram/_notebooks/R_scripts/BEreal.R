##### 
# /media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/
library("BEclear")
library("snowfall")

library("doParallel")
#cores<-detectCores()
#cl <- makeCluster(3, outfile="")
#registerDoParallel(cl)



print("..loading type data")
root_folder ="~/DATA/LungCancerResearch/" 
type_data <- read.table(paste(root_folder, "MethylationProbeTypes.csv", sep=""),
                        sep="\t", header=TRUE)
TypeI = as.character(type_data[type_data['Infinium_Design_Type']=='I', 'IlmnID'])
TypeII = as.character(type_data[type_data['Infinium_Design_Type']=='II', 'IlmnID'])


#####
print("..loading phenodata")
pheno_data <- read.table(paste(root_folder, "_prepped/methylation_combat_meta_small_cont.csv", sep=""),
                         sep="\t", header=TRUE)
pheno_data <- as.data.frame(pheno_data)
pheno_data <- pheno_data[!duplicated(pheno_data[[1]]),]
pheno_data[, 'Array.name'] <- gsub("-", ".", as.character(pheno_data[, 'Array.name']))
rownames(pheno_data) <- pheno_data[,1]

pheno_data$sample <- seq.int(nrow(pheno_data))
pheno_data$Sample.name<-NULL
pheno_data$sample<-NULL
#####
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
#####
# create sample_id, batch_id
sample_batch_df= pheno_data[c("Array.name", "Batch")]
names(sample_batch_df) <- c('sample_id', 'batch_id')
#####
print("..loading full data")
full_data <- read.table(paste(root_folder, "_prepped/methylation_raw_0nanmax.csv", sep=""),
                        sep="\t", 
                        header=TRUE)

full_data <- as.data.frame(full_data)
full_data <- full_data[!duplicated(full_data[[1]]),]
rownames(full_data) <- full_data[,1]
full_data[,1]<- NULL 

full_rows = rownames(full_data)
final_TypeI = Reduce(intersect, list(full_rows, TypeI))
final_TypeII = Reduce(intersect, list(full_rows, TypeII))

full_data1 = full_data[final_TypeI,]
full_data2 = full_data[final_TypeII,]
#####
rm(full_data)
print("..applying BEclear")
#registerDoParallel(cl)
#correctBatchEffect(data = full_data1, 
#                            samples = sample_batch_df, 
#                            adjusted = TRUE,
#			    colBlockSize = 0, 
#                            rowBlockSize = 0,
#                            epochs = 5,
#                            dir = paste(root_folder, "_prepped/BEclearType1", sep=""),
#                            outputFormat = 'txt')

rm(full_data1)
correctBatchEffect(data = full_data2, 
                            samples = sample_batch_df, 
                            adjusted = TRUE,
                            colBlockSize = 0, 
                            rowBlockSize = 0,
                            epochs = 5,
                            dir = paste(root_folder, "_prepped/BEclearType2", sep=""),
                            outputFormat = 'txt')
