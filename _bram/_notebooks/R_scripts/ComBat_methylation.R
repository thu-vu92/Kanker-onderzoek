# https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

##### 
type_data <- read.table("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/MethylationProbeTypes.csv", 
                        sep="\t", header=TRUE)
TypeI = as.character(type_data[type_data['Infinium_Design_Type']=='I', 'IlmnID'])
TypeII = as.character(type_data[type_data['Infinium_Design_Type']=='II', 'IlmnID'])


# keep only primary tumor
pheno_data <- read.table("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_meta_small_cont.csv",
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

full_data <- read.table("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_raw_0nanmax.csv",
                        sep="\t", 
                        header=TRUE)

full_data <- as.data.frame(full_data)
full_data <- full_data[!duplicated(full_data[[1]]),]
rownames(full_data) <- full_data[,1]
full_data[,1]<- NULL 

full_rows = rownames(full_data)
final_TypeI = Reduce(intersect, list(full_rows, TypeI))
final_TypeII = Reduce(intersect, list(full_rows, TypeII))

full_data1 <- full_data[, diag1]
full_data2 <- full_data[, diag2] 

full_matrix = as.matrix(full_data)




# make sure that the order of the samples is the same !!
####
library("sva", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
source("/home/koekiemonster/DEV/GIT/RexR/_hackathon2018/_notebooks/R_scripts/CComBat.R")

#corrected_batch_only_nonpar = CComBat(full_matrix, batch = batch, mod = NULL, par.prior = FALSE, prior.plots = TRUE, par.plots = TRUE)
#write.table(corrected_batch_only_nonpar, 
#            file="/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_0nans_primaryonly_nonparam.csv",
#            row.names=TRUE)

nullmode<-FALSE
splitByTarget<-FALSE
splitByType<-TRUE
param<-FALSE

if((splitByTarget==TRUE) & (splitByType==FALSE)){
  print("Split by target")
  full_matrix1 = as.matrix(full_data1)
  full_matrix2 = as.matrix(full_data2)
  #rm(full_matrix)
  #rm(full_data)
}else if((splitByTarget==FALSE)  &  (splitByType==TRUE)){
  print("Split by type")
  full_matrix1 = as.matrix(full_data[final_TypeI,])
  full_matrix2 = as.matrix(full_data[final_TypeII,])
  #rm(full_matrix)
  #rm(full_data)
}else if((splitByTarget==TRUE)  &  (splitByType==TRUE)){
  print("Split by target and  type")
  full_matrix11 = as.matrix(full_data1[final_TypeI,])
  full_matrix21 = as.matrix(full_data1[final_TypeII,])  
  full_matrix12 = as.matrix(full_data2[final_TypeI,])
  full_matrix22 = as.matrix(full_data2[final_TypeII,]) 
  #rm(full_matrix)
  #rm(full_data)
  }

  

if(param){parstring='param'}else{parstring='nonparam'}

    if (splitByTarget==FALSE && splitByType==FALSE){
      print("Full matrix")
        if(nullmode==FALSE){mod = modcombat}else{mod = mod0}
        corrected_batch_only_par = CComBat(full_matrix, batch = batch, mod = mod, par.prior = param)
        write.table(corrected_batch_only_par, 
                    file=paste("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_0nans_primaryonly_",
                               parstring, ".csv", sep=""),
                    row.names=TRUE)
    }else if (splitByTarget==TRUE && splitByType==TRUE){
      print("Split by target and type")
      if(nullmode==FALSE){mod = modcombat1}else{mod = mod01}
      corrected_batch_only_par = CComBat(full_matrix11, batch = batch1, mod = mod, par.prior = param)
      write.table(corrected_batch_only_par, 
                  file=paste("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_0nans_primaryonly_",
                             parstring, "_target1_type1.csv", sep=""),
                  row.names=TRUE)
      
      if(nullmode==FALSE){mod = modcombat1}else{mod = mod01}
      corrected_batch_only_par = CComBat(full_matrix21, batch = batch1, mod = mod, par.prior = param)
      write.table(corrected_batch_only_par, 
                  file=paste("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_0nans_primaryonly_",
                             parstring, "_target1_type2.csv", sep=""),
                  row.names=TRUE)   
      
      if(nullmode==FALSE){mod = modcombat2}else{mod = mod02}
      corrected_batch_only_par = CComBat(full_matrix12, batch = batch2, mod = mod, par.prior = param)
      write.table(corrected_batch_only_par, 
                  file=paste("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_0nans_primaryonly_",
                             parstring, "_target2_type1.csv", sep=""),
                  row.names=TRUE)
      
      if(nullmode==FALSE){mod = modcombat2}else{mod = mod02}
      corrected_batch_only_par = CComBat(full_matrix22, batch = batch2, mod = mod, par.prior = param)
      write.table(corrected_batch_only_par, 
                  file=paste("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_0nans_primaryonly_",
                             parstring, "_target2_type2.csv", sep=""),
                  row.names=TRUE)       
    }else{
        if(splitByTarget==TRUE){tstring<-'target'}else{tstring<-'type'}
      
        print(paste("Split by", tstring, sep=" "))
        
        if(nullmode==FALSE){if(splitByTarget==TRUE){mod = modcombat1; tbatch=batch1}else{mod = modcombat; tbatch=batch}}else{mod = mod01}
        
        corrected_batch_only_par = CComBat(full_matrix1, batch = tbatch, mod = mod, par.prior = param)
        write.table(corrected_batch_only_par, 
                    file=paste("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_0nans_primaryonly_",
                               parstring, "_",tstring,"1.csv", sep=""),
                    row.names=TRUE)
        
        if(nullmode==FALSE){if(splitByTarget==TRUE){mod = modcombat2; tbatch=batch2}else{mod = modcombat; tbatch=batch}}else{mod = mod02}
        corrected_batch_only_par = CComBat(full_matrix2, batch = tbatch, mod = mod, par.prior = param)
        write.table(corrected_batch_only_par, 
                    file=paste("/media/koekiemonster/DATA-FAST/genetic_expression/hackathon_2/Lung/_prepped/methylation_combat_0nans_primaryonly_",
                               parstring, "_",tstring,"2.csv", sep=""),
                  row.names=TRUE)
      }
 