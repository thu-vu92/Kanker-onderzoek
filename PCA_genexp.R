### 16 April 2019
### A. van Schetsen
### Re-run Figure 2, 5, 6 manusscript PCA gen_exp

library(data.table)
library(readr)
library(ggplot2)
library(Rmisc)
library(xlsx)

### SET DIRECTORY

setwd("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/R_script_github")

#### DATA INLEZEN, change input DT

DT <- fread("gen_exp.csv",sep = "~",dec=".")
DF = as.data.frame(DT)
List <- fread("COSTIM_AGPRES_COINHIB_CYTCHEM.csv",sep = "~")

List_names <- List[,Gen_exp]
List_costim <- List[Group=="COSTIM",Gen_exp]
List_agpres <- List[Group=="AGPRES",Gen_exp]

# MAKE DATA TABLE GENE GROUPS

DT_List <- DT[,c("SampleID",List_names),with=FALSE]
DT_List[,mean_ALL:=rowMeans(.SD,na.rm=TRUE),.SDcols = List_names]
DT_List[,mean_COSTIM:=rowMeans(.SD,na.rm=TRUE),.SDcols = List_costim]
DT_List[,mean_AGPRES:=rowMeans(.SD,na.rm=TRUE),.SDcols = List_agpres]

### PCA BEREKENEN

DF_NUM = DF[,-1]
m <- which(apply(DF_NUM, 2, var)!=0)
DF_NUM <- subset(DF_NUM,select=names(m))
PCA = prcomp(DF_NUM,scale. = T)

### PCA FEATURES MAKEN EN SAMENVOEGEN MET PHENOTYPEN

components = predict(PCA,DF_NUM)
PCA_components = cbind("SampleID" = DF[,1],as.data.frame(components))

#### ALLE PHENOTYPEN INLADEN

phtype = fread("Lung_Phenotype_Metadata.txt")
colnames(phtype)[4] = "Vital"
colnames(phtype)[11] = "Time"
colnames(phtype)[25] = "Tumor"
colnames(phtype) <-  gsub('([[:punct:]])|\\s+','_',colnames(phtype))
phtype <- phtype[Tumor != "NA"]
phtype[,Tumor2 := ifelse(grepl("iv", phtype$Tumor),"stage 4",
                         ifelse(grepl("iii", phtype$Tumor),"stage 3",
                                ifelse(grepl("ii", phtype$Tumor),"stage 2",
                                       "stage 1")))]
phtype$Tumor <- NULL
Pheno_with_pca =  merge(phtype,PCA_components,by = "SampleID",all = F)
#Pheno_with_pca = Pheno_with_pca[-which(Pheno_with_pca$PC2>1500),]
Pheno_with_pca[Vital%in%c("FFPE Scrolls","Recurrent Tumor"),Vital:="Primary Tumor"]

### SCATTERPLOT, Figure 2

ggplot(Pheno_with_pca , aes(x = PC2,y = PC3,colour = Vital)) + geom_point()
ggplot(Pheno_with_pca , aes(x = PC5,y = PC2,colour = Diagnosis)) + geom_point()

### LOADINGS (RELATIEVE IMPORTANTIE PER FEATURE) OF INDIVIDUAL GENES PER PCA

loadings = abs(PCA$rotation)
var.contrib <- data.frame(sweep(loadings, 2, colSums(loadings), "/"))
no_PCA = 10
var.contrib = var.contrib[,c(1:no_PCA)]

### PLOT BARGRAPH FOR EACH GENE GROUP

group = "COSTIM"
make_bar(group)

### BAR FUNCTION TO RETRIEVE BAR PLOT OF 1 GENE GROUP PER PCA
make_bar = function(group){
  genes = List[Group==group,Gen_exp]
  genes_index = c(grep(paste(unlist(genes),collapse="|") ,rownames(var.contrib)))
  data_group = stack(var.contrib[genes_index,])
  data_group$label <- group
  data_all = stack(var.contrib)
  data_all$label <- "_All genes"
  boxplotdata = rbind(data_all,data_group)
  myplot <- ggplot(boxplotdata) + geom_boxplot(aes(x = ind,y = values, fill=label)) + scale_y_continuous(labels = scales::percent) + labs(x="PCA",y="Loadings")
  return(myplot)
}

### MAKE PDF, Figure 5
myplots <- list()
myplots[[1]] <- make_bar("AGPRES")
myplots[[2]] <- make_bar("COSTIM")
myplots[[3]] <- make_bar("COINHIB")
myplots[[4]] <- make_bar("CYTCHEM")
pdf("C:/Users/tvu032/Desktop/Beyond_Banking/Results/test.pdf",height=8,width=15)
multiplot(plotlist = myplots, cols=2)
dev.off()

### TEST CONNECT PHENO WITH DATA TO PLOT DIAGNOSIS, Figure 6

plotdata <-  merge(DT_List[,.(SampleID,mean_COSTIM,mean_AGPRES)],phtype[,.(SampleID,Diagnosis)],by = "SampleID",all = F)
ggplot(plotdata , aes(x = mean_COSTIM,y = mean_AGPRES,colour=Diagnosis)) + geom_point()















