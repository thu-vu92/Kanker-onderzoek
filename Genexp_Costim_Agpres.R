library(data.table)
library(readr)
library(ggplot2)
library(Rmisc)
library(xlsx)

#### DATA INLEZEN

DT <- fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/DATA/gen_exp.csv",sep = "~",dec=".")
List <- fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/DATA/COSTIM_AGPRES.csv",sep = "~")
colnames(List)[3] = "Type"  #KUTSPATIE IN CSV BESTAND

List_names <- List[,Gen_exp]
List_costim <- List[Group=="COSTIM",Gen_exp]
List_agpres <- List[Group=="AGPRES",Gen_exp]

# MAKE DATA TABLE GENE GROUPS
DT_List <- DT[,c("SampleID",List_names),with=FALSE]
DT_List[,mean_ALL:=rowMeans(.SD,na.rm=TRUE),.SDcols = List_names]
DT_List[,mean_COSTIM:=rowMeans(.SD,na.rm=TRUE),.SDcols = List_costim]
DT_List[,mean_AGPRES:=rowMeans(.SD,na.rm=TRUE),.SDcols = List_agpres]

#### ALLE PHENOTYPEN INLADEN

phtype = fread("C://Users/aschetsen001/Documents/Projects@PwC/LungCancer/DATA/Lung_Phenotype_Metadata.txt")
colnames(phtype)[4] = "Vital"
colnames(phtype)[11] = "Time"
colnames(phtype)[25] = "Tumor"
colnames(phtype) <-  gsub('([[:punct:]])|\\s+','_',colnames(phtype))

### TEST CONNECT PHENO WITH DATA TO PLOT DIAGNOSIS

plotdata <-  merge(DT_List[,.(SampleID,mean_COSTIM,mean_AGPRES)],phtype[,.(SampleID,Diagnosis)],by = "SampleID",all = F)
ggplot(plotdata , aes(x = mean_COSTIM,y = mean_AGPRES,colour=Diagnosis)) + geom_point()

### We see a great correlation between the two gene groups
### FUNCTION TO MAKE SCATTERPLOT WITH COLOR CORREPSONDING PHENOTYPES

make_scatter = function(pheno){
  index = which(colnames(phtype)==pheno)
  myplot = 0
  if(length(index)){
    cols = c(1,index)
    plotdata <- merge(DT_List[,.(SampleID,mean_COSTIM,mean_AGPRES)],phtype[,..cols],by = "SampleID",all=FALSE)
    colnames(plotdata)[4] = "Type"
    myplot = ggplot(plotdata , aes(x = mean_COSTIM,y = mean_AGPRES,colour=Type)) + geom_point() + ggtitle(pheno) + theme(legend.position="none")
      #+ theme(legend.position="bottom")
     
  }
  return(myplot)
}

### MAKE LIST OF PHENOTYPES SCATTERPLOT

Scatter_pheno <- c("Diagnosis","Vital_Status","Gender","Pack_Years","Smoking_Status","Age_At_Diagnosis__Years_","Relapse_Status","Vital")
Scatter_pheno <- colnames(phtype[,2:41]) # DO ALL
myplots <- list()
for(i in Scatter_pheno) {
  myplots[[i]] <- make_scatter(i)
}
pdf("scatterplots.pdf",height=ceiling(length(Scatter_pheno)/2)*5,width=15)
multiplot(plotlist = myplots, cols=2)
dev.off()

### WHAT ABOUT THE TIME PHENOTYPE? SAYS NOTHING

plotdata <-  merge(DT_List[,.(SampleID,mean_COSTIM,mean_AGPRES)],phtype[,.(SampleID,Time)],by = "SampleID",all = F)
plotdata = plotdata[-which(plotdata$Time>1000),]
ggplot(plotdata , aes(x = mean_COSTIM,y = mean_AGPRES,colour=Time)) + geom_point()

### FUNCTION TO MAKE BOXPLOT OF COLOR

make_box = function(pheno){
  index = which(colnames(phtype)==pheno)
  myplot = 0
  if(length(index)){
    cols = c(1,index)
    plotdata <- merge(DT_List[,.(SampleID,mean_COSTIM,mean_AGPRES)],phtype[,..cols],by = "SampleID",all=FALSE)
    boxplot_data <- melt(plotdata,id = c("SampleID",pheno))
    colnames(boxplot_data)[2] = "Type"
    myplot = ggplot(boxplot_data) + geom_boxplot(aes(x = variable,y = value, fill = Type),size = .7) + theme(legend.position="bottom")
  }
  return(myplot)
}

### CHECK INFLUENCE OF COSTIM AND AGPRES IN BOXPLOT WITH A FILL OF THE PHENOTYPES

Boxplot_pheno <- c("Diagnosis","Vital_Status","Gender","Pack_Years","Smoking_Status","Age_At_Diagnosis__Years_","Relapse_Status","Vital")
myplots <- list()
for(i in Boxplot_pheno) {
  myplots[[i]] <- make_box(i)
}
pdf("boxplots.pdf",height=ceiling(length(Boxplot_pheno)/3)*5,width=15)
multiplot(plotlist = myplots, cols=2)
dev.off()


### CHECK OVERALL INFLUENCE DIFFERENT COMPONENTS INSIDE THE GENE GROUP

boxplot_data = stack(DT_List[,c(List_names),with=F])
setDT(boxplot_data)
boxplot_data[,group:="COSTIM"]
boxplot_data[ind %in% List_agpres,group:="AGPRES"]
ggplot(boxplot_data) + geom_boxplot(aes(x = ind,y = values,fill=group))  + theme(axis.text.x = element_text(angle = 45, hjust = 1))

###