library("org.Mm.eg.db")
library("xlsx")

setwd("C:/Users/asbarros/Desktop/Bioinfo/RNASeq_July19/R_Analysis/")
options(java.parameters = "- Xmx2048m")
source("scripts/DESeq_Model_NoLi61.R",echo = T)
source("scripts/theme_geometry_function.R")

cbPalette <- c("#999999","#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

##########################################################################################################################################################
dir.create(path = "DispersionPlots/", showWarnings = FALSE)

res1<-results(dm,contrast = c("Group", "PBS_Inf","PBS_Ctrl"),alpha = 0.05)
res2<-results(dm,contrast = c("Group", "Doxy_Ctrl","PBS_Ctrl"),alpha=0.05)
res3<-results(dm,contrast = c("Group", "Doxy_Inf","PBS_Inf"),alpha = 0.05)
res4<-results(dm,contrast = c("Group", "Doxy_Inf","Doxy_Ctrl"),alpha = 0.05)

res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res4$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res4),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res4$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res4), column="GENENAME", keytype="ENSEMBL", multiVals="first")

summary(res1)
summary(res2)
summary(res3)
summary(res4)

ashr_de1<-lfcShrink(dm, contrast = c("Group", "PBS_Inf","PBS_Ctrl"),type="ashr",res = res1)
ashr_de2<-lfcShrink(dm, contrast = c("Group", "Doxy_Ctrl","PBS_Ctrl"),type="ashr",res = res2)
ashr_de3<-lfcShrink(dm, contrast = c("Group", "Doxy_Inf","PBS_Inf"),type="ashr",res = res3)
ashr_de4<-lfcShrink(dm, contrast = c("Group", "Doxy_Inf","Doxy_Ctrl"),type="ashr",res = res4)

normal_de1<-lfcShrink(dm, contrast = c("Group", "PBS_Inf","PBS_Ctrl"),type="normal",res = res1)
normal_de2<-lfcShrink(dm, contrast = c("Group", "Doxy_Ctrl","PBS_Ctrl"),type="normal",res = res2)
normal_de3<-lfcShrink(dm, contrast = c("Group", "Doxy_Inf","PBS_Inf"),type="normal",res = res3)
normal_de4<-lfcShrink(dm, contrast = c("Group", "Doxy_Inf","Doxy_Ctrl"),type="normal",res = res4)

res1_sub<-res1[,-4]
res2_sub<-res2[,-4]
res3_sub<-res3[,-4]
res4_sub<-res4[,-4]

check1_ashr<-vector()
check2_ashr<-vector()
check3_ashr<-vector()

check1_normal<-vector()
check2_normal<-vector()
check3_normal<-vector()

for (i in 1:ncol(normal_de1)) {
  check1_normal[i]<-all(as.data.frame(res1)[,i] == as.data.frame(normal_de1)[,i],na.rm=T)
  check2_normal[i]<-all(as.data.frame(res2)[,i] == as.data.frame(normal_de2)[,i],na.rm=T)
  check3_normal[i]<-all(as.data.frame(res3)[,i] == as.data.frame(normal_de3)[,i],na.rm=T)
}

for (i in 1:ncol(ashr_de1)) {
  check1_ashr[i]<-all(as.data.frame(res1_sub)[,i] == as.data.frame(ashr_de1)[,i],na.rm=T)
  check2_ashr[i]<-all(as.data.frame(res2_sub)[,i] == as.data.frame(ashr_de2)[,i],na.rm=T)
  check3_ashr[i]<-all(as.data.frame(res3_sub)[,i] == as.data.frame(ashr_de3)[,i],na.rm=T)
}

print(check1_ashr)
print(check2_ashr)
print(check3_ashr)
print(check1_normal)
print(check2_normal)
print(check3_normal)

###############################################################################

Joint<-data.frame(GeneSymbol=res2$symbol,GeneName=res2$geneName,
                  pAdj_NI=ashr_de2$padj, log2FC_NI=ashr_de2$log2FoldChange,
                  pAdj_Inf=ashr_de3$padj,log2FC_Inf=ashr_de3$log2FoldChange)


Joint$Sign<-"NS"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf < 0.05]<-"Significant for Both"
Joint$Sign[Joint$pAdj_NI < 0.05 & Joint$pAdj_Inf >= 0.05]<-"Only Non-Infected"
Joint$Sign[Joint$pAdj_NI >= 0.05 & Joint$pAdj_Inf < 0.05]<-"Only Infected"
Joint$line<-seq(-30,30,length.out = nrow(Joint))

tiff("DispersionPlots/DP_Doxy_vs_PBS.tiff", width = 1500, height = 1500,res = 200)

ggplot(Joint) +
  theme_geometry(Joint$log2FC_NI, Joint$log2FC_Inf,interval = 8)+
  annotate("text", x = -0.1, y = -0.1, label = "0", size = 3)+
  annotate("text", x = 38, y = 0, label = "Log2FC in Non-Infection", fontface = 'italic',size = 4,angle = 90,)+
  annotate("text", x = 0, y = 38, label = "Log2FC in Infection", size = 4,fontface = 'italic')+
  geom_segment(aes(x = -32, y = 0, xend=36,yend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  geom_segment(aes(x = 0, y = -32, yend=36,xend=0),color = "black", size = 1,
               arrow=arrow(ends="last",length = unit(0.3, "cm")))+
  xlim(-40, 40)+ 
  geom_line(aes(x=line,y=line),linetype = "dashed",color="darkgrey")+
  ylim(-40, 40)+
  ggtitle("Dispersion Plot")+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),legend.position = "bottom")+
  geom_point(aes(x = log2FC_NI, y = log2FC_Inf,color=Sign), alpha=0.7, size=4)+
  scale_colour_manual(values=cbPalette)

dev.off()
################################################################################
save.image(file = paste0("environments/DispersonPlots_",Sys.Date(),".RData",sep=""))