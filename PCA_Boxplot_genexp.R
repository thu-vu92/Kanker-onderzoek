library(data.table)
library(readr)
library(ggplot2)
library(Rmisc)
library(xlsx)

#### DATA INLEZEN

DT <- fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/gen_exp.csv",sep = "~",dec=".")
DF = as.data.frame(DT)
load('genexp.RData')

### PCA BEREKENEN

DF_NUM = DF[,-1]
m <- which(apply(DF_NUM, 2, var)!=0)
DF_NUM <- subset(DF_NUM,select=names(m))
PCA = prcomp(DF_NUM,scale. = T)

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

components = predict(PCA,DF_NUM)
PCA_components = cbind("SampleID" = DF[,1],as.data.frame(components))

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
    myplot <- ggplot(boxplot_data) + geom_boxplot(aes(x = variable,y = value, fill = Label),size = .7) + ylim(c(round(min(boxplot_data[,4])),round(max(boxplot_data[,4])))) +
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

### GENE GROUPS INFLUENCE ON PRINCIPLE PCA

genegroup = fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/GeneGroup.csv",sep = "~")

### LOADINGS (RELATIEVE IMPORTANTIE PER FEATURE) OF INDIVIDUAL GENES PER PCA

loadings = abs(PCA$rotation)
var.contrib <- data.frame(sweep(loadings, 2, colSums(loadings), "/"))

### PLOT BARGRAPH FOR EACH GENE GROUP

groups = unique(genegroup[,3])
no_PCA = 10
myplots <- list()
for (j in 1:nrow(groups)) {
  myplots[[j]] <- make_bar(groups[[j,1]])
}
pdf("barplots_3.pdf",height=ceiling(nrow(groups)/2)*5,width=15)
multiplot(plotlist = myplots,cols=2)
dev.off()

### BAR FUNCTION TO RETRIEVE BAR PLOT OF 1 GENE GROUP PER PCA
make_bar = function(group){
  genes = genegroup[grep(group,genegroup$GROUP),1]
  genes_index = c(grep(paste(unlist(genes),collapse="|") ,rownames(var.contrib)))
  sum_array = colSums(var.contrib[genes_index,c(1:no_PCA)])/length(genes_index)
  sum = data.frame(PCA=names(sum_array), percentage=sum_array, row.names=NULL)
  myplot <- ggplot(stack(var.contrib[,1:no_PCA])) + geom_boxplot(aes(x = ind,y = values)) + geom_point(data=sum, aes(y=percentage, x=PCA), colour = "red", size = 2) + scale_y_continuous(labels = scales::percent) + labs(x=paste("Average influence of",group,"gene group on PCA"),y="Percentage")
  return(myplot)
}

### PRINT GENES

Given_genes = colnames(DF)
write.xlsx(Given_genes, file = "C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/test.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)
