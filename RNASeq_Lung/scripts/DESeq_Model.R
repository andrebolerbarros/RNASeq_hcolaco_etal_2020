rm(list=ls())

#Load the required libraries
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("org.Mm.eg.db")
library("ggplot2")
library("genefilter")
library("apeglm")

#Determine the ggplot2 predefined theme
theme_set(theme_bw())


##### Data Importation ############################################

#Determine the working directory
setwd("~/Colaco_etal_2019/RNASeq_Lung/")

#Create the required folders for the analysis
dir.create(path = "GeneralPlots/", showWarnings = FALSE)
dir.create(path = "SessionInfo/", showWarnings = FALSE)

#Save a file with the information of the sessionInfo
writeLines(capture.output(sessionInfo()),paste0("SessionInfo/SessionInfo_",Sys.Date(),".txt",sep=""))

#The ReadsPerGene files are saved in a folder called gene_counts. We need to save a list of files to import
gc_tab<-"gene_counts/"
file.names<-paste0(gc_tab,dir(gc_tab,pattern="tab"),sep="")

#Importing each file as data-frame and save them in a list
outfile<-list()
for (i in file.names) {
  outfile[[i]]<-read.table(i,row.names=1,comment.char="")
}

#Remove the first for lines, since they will not be used in this analysis
for (i in 1:length(outfile)) {
  outfile[[i]]<-outfile[[i]][-c(1:4),]
}

#Now, we have a data-frame for each sample. The column to be considered for each sample is based on the strandness of the sequencing protocol
tab<-data.frame(outfile[[1]][,3])
rownames(tab)<-rownames(outfile[[1]])

for (i in 2:length(outfile)) {
  tab[,i]<-outfile[[i]][,3]
  rownames(tab)<-rownames(outfile[[1]])
}

#Re-name the columns to get only the sample names and remove the path of each file & extension
colnames(tab)<-gsub("^.*(lg[5-8][a-f]+).*","\\1",names(outfile))

#Import the metadata file & add a new variable that combines treatment & infection
design<-read.csv2("sampleTable.csv",header=T,row.names = 1,sep=";")
design$group<-factor(paste(design$drugs,design$infection,sep="_")) 

#Check if the rownames of design match exactly w/ the column names. This is important since DESeq2 will not match metadata rownames with raw data colnames
all(rownames(design)==colnames(tab)) 

#Determine the baseline level for the group variable
design$group<-relevel(x = design$group,ref = "PBS_NonInf")

#Create a new variable for the LibrarySize of each sample
design$LibrarySize<-colSums(tab)

#Make a plot for the library sizes of each sample and save it as tiff
tiff("GeneralPlots/LibrarySize_Plot.tiff", width = 1500, height = 1000,res = 150)

ggplot(data=design, aes(x=rownames(design),y=(design$LibrarySize/1E6)))+
  geom_bar(stat="identity",fill="dodgerblue4")+
  ylab("Total Reads aligned to Genes (Million Reads)")+
  xlab("Samples")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()


##### DESeq Object & Normalization############################################

#Create the DESeq object
dm <- DESeqDataSetFromMatrix(countData = tab, colData = design, design = ~ group)

#Remove the rows with zero counts for all samples. This will help mainly for computational purposes
dm<-dm[rowSums(counts(dm)) > 0 , ]

#Estimate Size Factors for each sample
dm<-estimateSizeFactors(dm)

#Perform the regularized log transformation
dlog<- rlog(dm,blind=F)


##### Clustering Plots ############################################

#Make a heatmap of all the samples to assess clustering among conditions
tiff("GeneralPlots/Heatmap_AllSamples.tiff", width = 1500, height = 1000,res = 150)
ddist<-dist(t(assay(dlog)))
heatmap1<-pheatmap(ddist,cluster_rows=T,cluster_cols=T,annotation_col = design[,c("drugs","infection")])
print(heatmap1)

dev.off()

#Make a PCA to assess clustering of all samples, by plotting sample names with colors accordingly with the level in the category 'group'
pca1<-plotPCA(dlog, intgroup = c("group"),returnData=T)
percentVar<-round(100*attr(pca1,"percentVar"))

tiff("GeneralPlots/PCA_AllSamples_Names.tiff", width = 1500, height = 1000,res = 150)

data1<-ggplot(pca1,aes(PC1,PC2,colour=group))+
  geom_text(aes(colour=group),label=rownames(design))+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

data1
dev.off()

#Make a PCA to assess clustering of all samples, by plotting the samples with colors accordingly with the level in the category 'group'

tiff("GeneralPlots/PCA_AllSamples.tiff", width = 1500, height = 1000,res = 150)

data1<-ggplot(pca1,aes(PC1,PC2,colour=group))+
  geom_point(aes(colour=group),size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

data1
dev.off()

##### Running DESeq ###################################

#Confirm the variable of interest is indeed the one used for the DESeq design
design(dm)<- ~ group

#Run DESeq
dm<-DESeq(dm)

#Create the folders necessary to save the variables and environments produced in this script
dir.create(path = "variables/", showWarnings = FALSE)
dir.create(path = "environments/", showWarnings = FALSE)

##### Saving R Environment ###################################

#Save the R environment and the DESeq object
saveRDS(dm, paste0("variables/dm_",Sys.Date(),".rds",sep="")) 
save.image(file = paste0("environments/DESeq_Model_",Sys.Date(),".RData",sep=""))
