library(data.table)
library(readr)
library(ggplot2)
library(Rmisc)
library(xlsx)

#### DATA INLEZEN

DT1 <- fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/splitC/Lung_Methylation_prepped_chunk1.csv",sep = "~")
DT2 <- fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/splitC/Lung_Methylation_prepped_chunk2.csv",sep = "~")
DT3 <- fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/splitC/Lung_Methylation_prepped_chunk3.csv",sep = "~")
DT4 <- fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/splitC/Lung_Methylation_prepped_chunk4.csv",sep = "~")
DF = do.call("rbind", list(DT1, DT2, DT3, DT4))
load('Methy.RData')

### PCA BEREKENEN

PCA = prcomp(DF[,-1],scale. = T)

### BEREKEN STD EN VARIANCE
std_dev <- PCA$sdev
met_pca_var <- std_dev^2
prop_expl <- met_pca_var/sum(met_pca_var)

plot(prop_expl, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")
plot(cumsum(prop_expl), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")

### PCA FEATURES MAKEN EN SAMENVOEGEN MET PHENOTYPEN

components = predict(PCA,DF[,-1])
PCA_components = cbind(DF[,1],as.data.frame(components))

#### ALLE PHENOTYPEN INLADEN

phtype = fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/Lung_Phenotype_Metadata.txt")
colnames(phtype)[11] = "Time"
setnames(phtype,"Tumor Stage", "Tumor")
phtype <- phtype[Tumor != "NA"]
phtype[,Tumor2 := ifelse(grepl("iv", phtype$Tumor),"stage 4",
                         ifelse(grepl("iii", phtype$Tumor),"stage 3",
                                ifelse(grepl("ii", phtype$Tumor),"stage 2",
                                       "stage 1")))]
phtype$Tumor <- NULL
colnames(phtype)[4] = "Vital"
Pheno_with_pca =  merge(phtype,PCA_components,by = "SampleID",all = F)
Pheno_with_pca = Pheno_with_pca[-which(Pheno_with_pca$PC2>1500),]
colnames(Pheno_with_pca) <-  gsub('([[:punct:]])|\\s+','_',colnames(Pheno_with_pca))

### GENERALISATION

n <- length(phtype)
colm <- c(1:n)
mincolm <- -c(1:n)

### SCATTERPLOT

p12 = ggplot(Pheno_with_pca , aes(x = PC1,y = PC2,colour = Diagnosis)) + geom_point() + theme(legend.position="none")
p34 = ggplot(Pheno_with_pca , aes(x = PC3,y = PC4,colour = Diagnosis)) + geom_point()+ theme(legend.position="none")
p56 = ggplot(Pheno_with_pca , aes(x = PC5,y = PC6,colour = Diagnosis)) + geom_point()+ theme(legend.position="none")
p78 = ggplot(Pheno_with_pca , aes(x = PC7,y = PC8,colour = Diagnosis)) + geom_point()+ theme(legend.position="none")

multiplot(p12,p34,p56,p78,cols = 2)

#### TEST BOXPLOT VOOR DIAGNOSIS

library(reshape)
sequence = seq(n+1,n+10)
boxplot_data = melt(Pheno_with_pca[,c(c(1,3),..sequence)],id = c("SampleID","Diagnosis"))
ggplot(boxplot_data) + geom_boxplot(aes(x = variable,y = value,fill = Diagnosis),size = .7) + ylim(c(-10,10)) +
  labs(y = "Value",x = "Principle component",fill = "Type")

#### FUNCTIE OM VOOR IEDER  PHENOTYPE EEN BOXPLOT TE MAKEN

make_box = function(pheno,no_PCA){
  index = which(colnames(Pheno_with_pca)==pheno)
  myplot = 0
  if(length(index)){
    type = sapply(Pheno_with_pca[,..index],typeof)
    # Scale no_PCA
    if(missing(no_PCA)){
      if(type == "double" || type == "integer"){
        no_PCA = 8
      }else{
        no_PCA = ceiling(20/nrow(unique(Pheno_with_pca[,c(..index)])))
      }
    }
    colm = c(1,index)
    sequence = seq(n+1,n+no_PCA)
    boxplot_data = melt(Pheno_with_pca[,c(..colm,..sequence)],id = c("SampleID",pheno))
    colnames(boxplot_data)[2] = "Label"
    if(type == "double" || type == "integer"){
      min = round(min(boxplot_data[,2], na.rm=T))
      int = (round(max(boxplot_data[,2], na.rm=T))-min)/3
      boxplot_data[,2] <- data.frame(cut(as.numeric(unlist(boxplot_data[,2])), breaks=c(min, min+int, min+int+int, Inf), labels=c("low","middle","high")))
    }
    myplot <- ggplot(boxplot_data) + geom_boxplot(aes(x = variable,y = value, fill = Label),size = .7)
      + ylim(c(round(min(boxplot_data[,4])),round(max(boxplot_data[,4])))) +
      labs(y = "Value",x = paste("PC - ",pheno),fill = "Type") + theme(legend.position="bottom")
  }
  return(myplot)
}

### ALL DIFFERENT BOX PLOTS, MAKE LIST OF DESIRED PHENOTYPES IN BOX_PHENO

Box_pheno <- c("Diagnosis","Vital_Status","Gender","Ethnicity","Pack_Years","Smoking_Status","Age_At_Diagnosis__Years_","Relapse_Status","Vital")
myplots <- list()
for(i in 1:length(Box_pheno)) {
  myplots[[i]] <- make_box(Box_pheno[i])
}
pdf("boxplots.pdf",height=ceiling(length(Box_pheno)/2)*5,width=15)
multiplot(plotlist = myplots, cols=2)
dev.off()

### LOADINGS (RELATIEVE IMPORTANTIE PER FEATURE) OF INDIVIDUAL METHYLATIONS PER PCA
### TAKE ONLY ABSOLUTE LOADINGS OF A FIXED NUMBER OF PCA'S TO REDUCE CALCULATION TIME

no_PCA = 10
loadings = as.data.frame(PCA$rotation[,1:no_PCA])
loadings = abs(loadings)
loadings_sorted = loadings[order(-loadings$PC1),]
plot(loadings$PC1,type = "b")
plot(loadings_sorted$PC1, type = "b")

### MOVING AVERAGE OVER METHYLATION DATA, LENGTH INTERVAL N
### Change PC1 in definition of MA to get other principle components

library(TTR)
N = 1000
MA = as.data.frame(SMA(loadings$PC1,N))    # CHANGE PC HERE
MA = as.data.frame(MA[complete.cases(MA), ])
rownames(MA) = paste(rownames(loadings)[1:(nrow(loadings)-N)],rownames(loadings)[N:nrow(loadings)])
colnames(MA)[1] = "PC"
MA_sorted = MA[order(-MA$PC),]
plot(MA_sorted, type = "b")
plot(MA$PC,xlab="Methylation data", ylab=paste("Loading PC 2 and MA interval",N))

### LAST M PART OF METHYLATION DATA, BOTH ABSOLUTE AS MOVING AVERAGE

M = 10000
loadings_last = loadings[(nrow(loadings)-M+1):(nrow(loadings)),]
plot(loadings_last$PC1)
loadings_last_sorted = loadings_last[order(-loadings_last$PC1),]
plot(loadings_last_sorted$PC1, type = "b")

N = 100
MA_last = as.data.frame(SMA(loadings_last$PC1,N))
MA_last = as.data.frame(MA_last[complete.cases(MA_last), ])
rownames(MA_last) = paste(rownames(loadings_last)[1:(nrow(loadings_last)-N)],rownames(loadings_last)[N:nrow(loadings_last)])
colnames(MA_last)[1] = "PC1"
MA__last_sorted = MA_last[order(-MA_last$PC1),]
plot(MA__last_sorted, type = "b")
plot(MA_last$PC1,xlab="Methylation data", ylab="Loading")

### Save environment

save.image("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/Methy.RData")

