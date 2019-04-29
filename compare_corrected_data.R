library(data.table)
library(ggplot2)
require(gridExtra)
library(stringr)

setwd("C:/Users/tvu032/Desktop/Beyond_Banking/Data/corrected")

############################# Methyl ###############################
load("C:/Users/tvu032/Desktop/Beyond_Banking/Kanker-onderzoek/prep_methy.RData")

dt_new <- fread("methyl_prep_new.csv")
dt <- copy(DT_methy_prep)

dim(dt_new)
dim(dt)

# Kernel density plot
par(mfrow=c(1,2))

rand_probe <- sample(names(dt),1)
d <- density(dt[[rand_probe]]) # returns the density data 
plot(d, main=paste("Old data - Methyl. density probeID", rand_probe))

d_new <- density(dt_new[[rand_probe]]) # returns the density data 
plot(d_new, main=paste("Corrected data - Methylation density of probeID", rand_probe))

# Plot against cancer subtypes
pheno <- fread("C:/Users/tvu032/Desktop/Beyond_Banking/Hackathon/Lung_Phenotype_Metadata.txt")

dt <- merge(dt, pheno[, .(SampleID, Diagnosis)], by="SampleID")
dt_new <- merge(dt_new, pheno[, .(SampleID, Diagnosis)], by="SampleID")

rand_probe <- sample(names(dt),1)
rand_probe2 <- sample(names(dt),1)

plot1 <- ggplot(dt, aes_string(x=rand_probe, y=rand_probe2, color="Diagnosis")) +
  geom_point() +
  labs(title="Old methylation data")

plot2 <- ggplot(dt_new, aes_string(x=rand_probe, y=rand_probe2, color="Diagnosis")) +
  geom_point() + 
  labs(title="New methylation data")

grid.arrange(plot1, plot2, ncol=2)



############################# Genex ###############################

dt_new <- fread("RNAexpression_combat_primaryonly_gender_nonparam.csv")
dt <- fread("C:/Users/tvu032/Desktop/Beyond_Banking/Hackathon/gen_exp.csv")

# Transpose new dataste
# Concatenate gene and chromosome and Start position
dt_new[,chrome := str_sub(chr, 4)]
dt_new[,gen_chr := paste(gene, chrome, start, sep="_")]

# Remove non-numeric cols
dt_new_copy <- copy(dt_new)
sampleID <- dt_new_copy$gen_chr
dt_new_copy[, ((length(dt_new_copy)-7) : length(dt_new)) := NULL] 
dt_new_copy$gen_chr <- NULL

# Transpose data
cols <- names(dt_new_copy)
dt_new_tran <- transpose(dt_new_copy)
colnames(dt_new_tran) <- sampleID
dt_new_tran <- cbind(cols, dt_new_tran)

dupcol <- names(dt_new_tran)[duplicated(names(dt_new_tran))]
dt_new_tran <- dt_new_tran[,names(dt_new_tran) %!in% dupcol, with=F]

rem = find_na(dt_new_tran)
names(dt_new_tran)[1] <- "SampleID"

dt_new <- copy(dt_new_tran)

fwrite(dt_new, "genex_prep_new")

# Check dims
dim(dt)
dim(dt_new)

length(setdiff(names(dt), names(dt_new)))
length(setdiff(names(dt_new), names(dt)))

com <- intersect(names(dt), names(dt_new))[-1]

# Kernel density plot
par(mfrow=c(1,2))

rand_probe <- sample(com,1)
d <- density(dt[[rand_probe]]) # returns the density data 
plot(d, main=paste("Old data - gene exp. density probeID", rand_probe))

d_new <- density(dt_new[[rand_probe]]) # returns the density data 
plot(d_new, main=paste("New data - gene exp. density probeID", rand_probe))


# Plot against cancer subtypes
pheno <- fread("C:/Users/tvu032/Desktop/Beyond_Banking/Hackathon/Lung_Phenotype_Metadata.txt")

dt <- merge(dt, pheno[, .(SampleID, Diagnosis)], by.x = "cols", by.y = "SampleID")
dt_new <- merge(dt_new, pheno[, .(SampleID, Diagnosis)], by = "SampleID")

dt_copy <- copy(dt)
dt_new_copy <- copy(dt_new)

names(dt_copy) <- gsub( "\\.", "_", names(dt_copy))
names(dt_new_copy) <-  gsub( "\\.", "_", names(dt_new_copy))

names(dt_copy) <- gsub( "-", "_", names(dt_copy))
names(dt_new_copy) <-  gsub( "-", "_", names(dt_new_copy))

com <- intersect(names(dt_copy), names(dt_new_copy))[-1]
rand_probe <- sample(com,1)
rand_probe2 <- sample(com,1)

plot1 <- ggplot(dt_copy, aes_string(x=rand_probe, y=rand_probe2, color="Diagnosis")) +
  geom_point() +
  labs(title="Old gene expression data")

plot2 <- ggplot(dt_new_copy, aes_string(x=rand_probe, y=rand_probe2, color="Diagnosis")) +
  geom_point() + 
  labs(title="New gene expression data")

grid.arrange(plot1, plot2, ncol=2)

