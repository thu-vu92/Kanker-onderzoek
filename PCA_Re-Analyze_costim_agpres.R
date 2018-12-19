library(data.table)
library(readr)
library(ggplot2)
library(Rmisc)
library(xlsx)

#### DATA INLEZEN

DT <- fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/DATA/gen_exp.csv",sep = "~",dec=".")
List <- fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/DATA/COSTIM_AGPRES.csv",sep = "~")
DF = as.data.frame(DT)

### PCA BEREKENEN

DF_NUM = DF[,-1]
m <- which(apply(DF_NUM, 2, var)!=0)
DF_NUM <- subset(DF_NUM,select=names(m))
PCA = prcomp(DF_NUM,scale. = T)

### LOADINGS (RELATIEVE IMPORTANTIE PER FEATURE) OF INDIVIDUAL GENES PER PCA

loadings = abs(PCA$rotation)
var.contrib <- data.frame(sweep(loadings, 2, colSums(loadings), "/"))
no_PCA = 10
var.contrib = var.contrib[,c(1:no_PCA)]

### PLOT BARGRAPH FOR EACH GENE GROUP

group = "AGPRES"
make_bar(group)

### BAR FUNCTION TO RETRIEVE BAR PLOT OF 1 GENE GROUP PER PCA
make_bar = function(group){
  genes = List[Group==group,Gen_exp]
  genes_index = c(grep(paste(unlist(genes),collapse="|") ,rownames(var.contrib)))
  sum_array = colSums(var.contrib[genes_index,])/length(genes_index)
  sum = data.frame(PCA=names(sum_array), percentage=sum_array, row.names=NULL)
  data_group = stack(var.contrib[genes_index,])
  myplot <- ggplot(stack(var.contrib[,1:no_PCA])) + geom_boxplot(aes(x = ind,y = values)) + geom_point(data=sum, aes(y=percentage, x=PCA), colour = "red", size = 2) + scale_y_continuous(labels = scales::percent) + labs(x=paste("Average influence of",group,"gene group on PCA"),y="Percentage")
  return(myplot)
}

### SEPERATE

genes_agpres = List[Group=="AGPRES",Gen_exp]
index_agpres =c(grep(paste(unlist(genes_agpres),collapse="|") ,rownames(var.contrib)))
data_agpres = stack(var.contrib[genes_index,])
data_agpres$label <- "AGPRES"

genes_costim = List[Group=="COSTIM",Gen_exp]
index_costim =c(grep(paste(unlist(genes_costim),collapse="|") ,rownames(var.contrib)))
data_costim = stack(var.contrib[index_costim,])
data_costim$label <- "COSTIM"

data_all = stack(var.contrib)
data_all$label <- "_all"

boxplotdata = rbind(data_all,data_costim,data_agpres)
ggplot(boxplotdata) + geom_boxplot(aes(x = ind,y = values, fill=label)) + scale_y_continuous(labels = scales::percent) + labs(x="PCA",y="Percentage")

### MAKE PDF
myplots <- list()
myplots[[1]] <- make_bar("AGPRES")
myplots[[2]] <- make_bar("COSTIM")
pdf("test.pdf",height=8,width=15)
multiplot(plotlist = myplots, cols=2)
dev.off()

