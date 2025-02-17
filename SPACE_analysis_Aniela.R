#Implemented a few optional data imputation approaches for masspec missing values but found little to no use
#in tested data sets when used conservatively (replace max 3 NAs = 50% of values per protein)
#Adjusted pVal for multiple testing topTable (limma) remove BH if no significant hits detected

setwd("~/Documents/ARGH/SPACE/")
library(stringr)
library(ggplot2)
library(ggrepel)
library(limma)
library(gplots)
library(plotly)
library(goeveg)
library(ggfortify)
library(factoextra)
library(cluster)
library(dplyr)
library(tidyr)

#If input is Sinas MaxQuant database search result file use below code chunk

#write.table(data,"data_clean.txt", sep = "\t", row.names = F)

#remove data that is only identified ny site and potential contaminants
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#data <- read.csv("proteinGroups_FDR_05.txt", header = T, fill = TRUE, sep="\t")
#data<- subset(data, data$Only.identified.by.site != "+" & data$Reverse != "+" & data$Potential.contaminant != "+")

#Annotation gymnastics
#Uniprot.ID = str_split_fixed(data$Protein.IDs, ";", 2)
#Uniprot.ID = as.data.frame(Uniprot.ID)
#Uniprot.ID = str_split_fixed(Uniprot.ID$V1, "-", 2)
#Uniprot.ID = as.data.frame(Uniprot.ID)
#data = cbind(data, Uniprot.ID[,1]) 
#colnames(data)[ncol(data)] <- "Uniprot.ID"
#Gene.names <- str_split_fixed(data$Gene.names, ";", 2)
#Gene.names <- as.data.frame(Gene.names)
#Gene.names[Gene.names==""] <- NA
#data = cbind(data, Gene.names[,1]) 
#colnames(data)[ncol(data)] <- "Gene.names.1st"

#calculate the number of NAs per treatment group
#data1 <- data[,c(1,7,108:121)]
#from here jump to calculate number of NAs per treatment group


#***If input is JBC proteome discoverer database search result file use below code chunk
#extract abundances or ibaqs from data file
filen = "PR758_Aniela_20231018_all_protein.csv"
data <- read.csv(filen, header = T, fill = TRUE, sep=",")
data1 = data[,c(4,5,93:111)]
colnames(data1)
data1$Description = sub(".*GN=", "", data1$Description)
data1$Description = sub(" PE=.*", "", data1$Description)
data1$Description = sub(" OS=.*", "", data1$Description)

head(data1$Description)

stra = c(rep("PP2_OSMI_",3),rep("PP2_TMG_",3),rep("global_ctrl",3),rep("global_OSMI",3),rep("PP2_ctrl_",3),"PP2_IgG",rep("global_TMG",3))
strb = c(rep(1:3,5),1,1:3)
colnames(data1)=c("Uniprot.ID","Gene.names",paste0(stra,strb))
head(data1)
data1 = data1[,c(1,2,15:17,3:8,9:14,19:21)]
#***

#calculate the number of NAs per treatment group
data1[,c(3:20)][data1[,c(3:20)]==0] <- NA
data1$NA.PP2_C = rowSums(is.na(data1[,c(3:5)]))
data1$NA.PP2_O = rowSums(is.na(data1[,c(6:8)]))
data1$NA.PP2_T = rowSums(is.na(data1[,c(9:11)]))
data1$NA.global_C = rowSums(is.na(data1[,c(12:14)]))
data1$NA.global_O = rowSums(is.na(data1[,c(15:17)]))
data1$NA.global_T = rowSums(is.na(data1[,c(18:20)]))

#log2 transform IBAQs/Abundances and do quantile normalisation within Pulldowns but accross treatment groups
data1[,c(3:20)] = log2(data1[,c(3:20)])
y = data1[,c(3:11)]
plotDensities(y)
data1[,c(3:11)] = normalizeQuantiles(as.matrix(data1[,c(3:11)]))
x = data1[,c(3:11)]
plotDensities(x)

y = data1[,c(12:20)]
plotDensities(y)
data1[,c(12:20)] = normalizeQuantiles(as.matrix(data1[,c(12:20)]))
x = data1[,c(12:20)]
plotDensities(x)

#For datasets where there are n=3 measurements in one treatment group and 3 NAs in the other, SR imputates the missing values
#NA by using the lowest measured peptide intensity within treatment group this is to prevent division/0
for (i in c(3:5)){
  data1[,i]= ifelse(data1$NA.PP2_C == 3 &  data1$NA.PP2_O == 0 , min(data1[,3:5], na.rm = T) , data1[,i])
}

for (i in c(6:8)){
  data1[,i]= ifelse((data1$NA.PP2_O == 3 &  data1$NA.PP2_C == 0), min(data1[,6:8], na.rm = T) , data1[,i])
}

for (i in c(9:11)){
  data1[,i]= ifelse((data1$NA.PP2_T == 3 &  data1$NA.PP2_C == 0), min(data1[,9:11], na.rm = T) , data1[,i])
}

for (i in c(3:5)){
  data1[,i]= ifelse(data1$NA.PP2_C == 3 &  data1$NA.PP2_T == 0 , min(data1[,3:5], na.rm = T) , data1[,i])
}

for (i in c(12:14)){
  data1[,i]= ifelse(data1$NA.global_C == 3 &  data1$NA.global_O == 0 , min(data1[,12:14], na.rm = T) , data1[,i])
}

for (i in c(15:17)){
  data1[,i]= ifelse((data1$NA.global_O == 3 &  data1$NA.global_C == 0) , min(data1[,15:17], na.rm = T) , data1[,i])
}

for (i in c(12:14)){
  data1[,i]= ifelse(data1$NA.global_C == 3 &  data1$NA.global_T == 0 , min(data1[,12:14], na.rm = T) , data1[,i])
}

for (i in c(18:20)){
  data1[,i]= ifelse(data1$NA.global_T == 3 &  data1$NA.global_C == 0 , min(data1[,18:20], na.rm = T) , data1[,i])
}

#(Test) normalisation to (IPed) protein: Use log and quantile normalised abundance/IBAQ and average multiple peptides, then substract from dataframe
pola = as.numeric(data1[data1$Gene.names=="POLR2A",3:20])
polb = as.numeric(data1[data1$Gene.names=="POLR2B",3:20])
pol = rowMeans(cbind(pola,polb))
tasty = sweep(data1[,3:20],2,pol,'-')
data1[,3:20] = tasty

#check if values are normally distributed
ggplot(gather(data1[,3:11]), aes(value)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~key, scales = 'free_x')

#Data imputation (see end of this script)

#run PCA and check for outliers
test = na.omit(data1[,c(2,3:11)])
#start here when using imputed data (end of script)
#test2 = t(test[,2:7]) 
test2 = t(test[,2:10]) 
colnames(test2)= test[,1]
test2 = data.frame(test2)
head(test2)
test2$treatment = as.factor(c(rep("PP2_Ctrl",3),rep("PP2_OSMI",3),rep("PP2_TMG",3)))
#test2$treatment = as.factor(c(rep("8WG16_dTAG",3),rep("8WG16_DMSO",3)))
pca.space <- prcomp(test2[,c(1:ncol(test2)-1)], center = T, scale. = T)
summary(pca.space)
autoplot(pca.space, data = test2, colour = 'treatment',label = TRUE, label.label = rownames(test2))


#Select replicates to run in treatment groups and build linear model
samps <- factor(rep(c("PP2_Ctrl","PP2_OSMI"), each = 3))
design = model.matrix(~0 + samps)
colnames(design) <- gsub("samps", "", colnames(design))
contrast <- makeContrasts(log2FC.PP2.OSMI = (PP2_OSMI - PP2_Ctrl),
                           levels = design)
fit <- lmFit(data1[,c(3:8)], design)#
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)#
volcanoplot(fit)
tt = topTable(fit,number=Inf, fit$genes, adjust.method ="BH",confint=T,sort.by='none') 
data2 = cbind(data1, tt[,c(1,6,7)])

colnames(data2)[(ncol(data2)-2):ncol(data2)] = paste0(colnames(data2)[(ncol(data2)-2):ncol(data2)], ".PP2.OSMI")

#Prepare data for volcano plot. Add gene names to hits and change colours. Set adj.p.value and FC cutoffs below

co_pVal <- 0.05
co_FC <- 1

data2$delabel <- NA
data2$colour <- "Black"

#for differentially changed proteins add gene.name entry
data2[which(data2$adj.P.Val.PP2.OSMI < co_pVal & data2$logFC.PP2.OSMI > co_FC) ,"delabel"] <- data2[which(data2$adj.P.Val.PP2.OSMI < co_pVal & data2$logFC.PP2.OSMI > co_FC) ,"Gene.names"]
data2[which(data2$adj.P.Val.PP2.OSMI < co_pVal & data2$logFC.PP2.OSMI < -co_FC) ,"delabel"] <- data2[which(data2$adj.P.Val.PP2.OSMI < co_pVal & data2$logFC.PP2.OSMI < -co_FC) ,"Gene.names"]
#for differentially changed proteins add colour label
data2[which(data2$adj.P.Val.PP2.OSMI < co_pVal & data2$logFC.PP2.OSMI > co_FC) ,"colour"] <- "Blue"
data2[which(data2$adj.P.Val.PP2.OSMI < co_pVal & data2$logFC.PP2.OSMI < -co_FC) ,"colour"] <- "Red"

data2 = data2[order(data2$adj.P.Val.PP2.OSMI),]
head(data2)

#Volcanoplot
png(paste0("Vplot_","PP2-SPACE",colnames(design)[2],"vs",colnames(design)[1],".png"), w=642, h=407, pointsize=20)
ggplot(data=data2, aes(x=logFC.PP2.OSMI, y=-log10(adj.P.Val.PP2.OSMI), label=delabel)) +
  geom_point(colour=data2$colour) + 
  theme_minimal() +
  geom_text_repel() +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red") +
  ggtitle(paste("ChIP-SPACE",colnames(design)[2],"vs",colnames(design)[1]))
dev.off()


samps <- factor(c(rep("PP2_Ctrl", each = 3), rep("PP2_TMG", each = 3)))
design = model.matrix(~0 + samps)
colnames(design) <- gsub("samps", "", colnames(design))
contrast <- makeContrasts( log2FC.PP2.TMG = (PP2_TMG - PP2_Ctrl),
                           levels = design)
fit <- lmFit(data2[,c(3:5,9:11)], design)#
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)#
volcanoplot(fit)
tt = topTable(fit, fit$genes, adjust="BH",number=Inf,confint=T,sort.by='none') 
data3 = cbind(data2, tt[,c(1,6,7)])

colnames(data3)[(ncol(data3)-2):ncol(data3)] = paste0(colnames(data3)[(ncol(data3)-2):ncol(data3)], ".PP2.TMG")

data3$delabel <- NA
data3$colour <- "Black"
data3[which(data3$adj.P.Val.PP2.TMG < co_pVal & data3$logFC.PP2.TMG > co_FC) ,"delabel"] <- data3[which(data3$adj.P.Val.PP2.TMG < co_pVal & data3$logFC.PP2.TMG > co_FC) ,"Gene.names"]
data3[which(data3$adj.P.Val.PP2.TMG < co_pVal & data3$logFC.PP2.TMG < -co_FC) ,"delabel"] <- data3[which(data3$adj.P.Val.PP2.TMG < co_pVal & data3$logFC.PP2.TMG < -co_FC) ,"Gene.names"]
#for differentially changed proteins add colour label
data3[which(data3$adj.P.Val.PP2.TMG < co_pVal & data3$logFC.PP2.TMG > co_FC) ,"colour"] <- "Blue"
data3[which(data3$adj.P.Val.PP2.TMG < co_pVal & data3$logFC.PP2.TMG < -co_FC) ,"colour"] <- "Red"


data3 = data3[order(data3$adj.P.Val.PP2.TMG),]
head(data3)

png(paste0("Vplot_","PP2-SPACE",colnames(design)[2],"vs",colnames(design)[1],".png"), w=642, h=407, pointsize=20)
ggplot(data=data3, aes(x=logFC.PP2.TMG, y=-log10(adj.P.Val.PP2.TMG), label=delabel)) +
  geom_point(colour=data3$colour) + 
  theme_minimal() +
  geom_text_repel() +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red") +
  ggtitle(paste("ChIP-SPACE",colnames(design)[2],"vs",colnames(design)[1]))
dev.off()

write.csv(data3, paste0(filen,"_proteinGroups_analzyed_ProtDisc_ouput_norm_POLR2A_and_POLR2B.csv"))


#Data imputation optional below using data1 df run seperately per different IP
colnames(data1)
dfSPACE_imputed <- data.frame(
  original = data1[,11:16],
  imputed_missForest = missForest(data1[,11:16])$ximp,
  imputed_pmm = complete(mice(data1[,11:16], method = "pmm")),
  imputed_cart = complete(mice(data1[,11:16], method = "cart")),
  imputed_lasso = complete(mice(data1[,11:16], method = "lasso.norm")),
  imputed_rf = complete(mice(data1[,11:16], method = "rf")))

colnames(dfSPACE_imputed)
#Imputation per IP across treatments (eg. quantile normalised log data)
dfSPACE_imputed_comb <- data.frame(
  original <- c(dfSPACE_imputed[,1],dfSPACE_imputed[,2],dfSPACE_imputed[,3],dfSPACE_imputed[,4],dfSPACE_imputed[,5],dfSPACE_imputed[,6]),
  imputed_missForest <- c(dfSPACE_imputed[,7],dfSPACE_imputed[,8],dfSPACE_imputed[,9],dfSPACE_imputed[,10],dfSPACE_imputed[,11],dfSPACE_imputed[,12]),
  imputed_pmm <- c(dfSPACE_imputed[,13],dfSPACE_imputed[,14],dfSPACE_imputed[,15],dfSPACE_imputed[,16],dfSPACE_imputed[,17],dfSPACE_imputed[,18]),
  imputed_cart <- c(dfSPACE_imputed[,19],dfSPACE_imputed[,20],dfSPACE_imputed[,21],dfSPACE_imputed[,22],dfSPACE_imputed[,23],dfSPACE_imputed[,24]),
  imputed_lasso <- c(dfSPACE_imputed[,25],dfSPACE_imputed[,26],dfSPACE_imputed[,27],dfSPACE_imputed[,28],dfSPACE_imputed[,29],dfSPACE_imputed[,30]),
  imputed_rf <- c(dfSPACE_imputed[,31],dfSPACE_imputed[,32],dfSPACE_imputed[,33],dfSPACE_imputed[,34],dfSPACE_imputed[,35],dfSPACE_imputed[,36])
)
#dTAG treatments
#dfSPACE_imputed_comb <- data.frame(
#  original <- c(dfSPACE_imputed[,1],dfSPACE_imputed[,2],dfSPACE_imputed[,3]),
#  imputed_missForest <- c(dfSPACE_imputed[,7],dfSPACE_imputed[,8],dfSPACE_imputed[,9]),
#  imputed_pmm <- c(dfSPACE_imputed[,13],dfSPACE_imputed[,14],dfSPACE_imputed[,15]),
#  imputed_cart <- c(dfSPACE_imputed[,19],dfSPACE_imputed[,20],dfSPACE_imputed[,21]),
#  imputed_lasso <- c(dfSPACE_imputed[,25],dfSPACE_imputed[,26],dfSPACE_imputed[,27]),
#  imputed_rf <- c(dfSPACE_imputed[,31],dfSPACE_imputed[,32],dfSPACE_imputed[,33])
#)
#DMSO treatments
#dfSPACE_imputed_comb <- data.frame(
#  original <- c(dfSPACE_imputed[,4],dfSPACE_imputed[,5],dfSPACE_imputed[,6]),
#  imputed_missForest <- c(dfSPACE_imputed[,10],dfSPACE_imputed[,11],dfSPACE_imputed[,12]),
#  imputed_pmm <- c(dfSPACE_imputed[,16],dfSPACE_imputed[,17],dfSPACE_imputed[,18]),
#  imputed_cart <- c(dfSPACE_imputed[,22],dfSPACE_imputed[,23],dfSPACE_imputed[,24]),
#  imputed_lasso <- c(dfSPACE_imputed[,28],dfSPACE_imputed[,29],dfSPACE_imputed[,30]),
#  imputed_rf <- c(dfSPACE_imputed[,34],dfSPACE_imputed[,35],dfSPACE_imputed[,36])
#)

#visualise imputed data distributions to pick the most appropriate method
h1 <- ggplot(dfSPACE_imputed_comb, aes(x = original)) +
  geom_histogram(fill = "#ad1538", color = "#000000", position = "identity") +
  ggtitle("Original distribution") +
  theme_classic()
h2 <- ggplot(dfSPACE_imputed_comb, aes(x = imputed_missForest)) +
  geom_histogram(fill = "#15ad4f", color = "#000000", position = "identity") +
  ggtitle("missForest-imputed distribution") +
  theme_classic()
h3 <- ggplot(dfSPACE_imputed_comb, aes(x = imputed_pmm)) +
  geom_histogram(fill = "#1543ad", color = "#000000", position = "identity") +
  ggtitle("pmm-imputed distribution") +
  theme_classic()
h4 <- ggplot(dfSPACE_imputed_comb, aes(x = imputed_cart)) +
  geom_histogram(fill = "red", color = "#000000", position = "identity") +
  ggtitle("cart-imputed distribution") +
  theme_classic()
h5 <- ggplot(dfSPACE_imputed_comb, aes(x = imputed_lasso)) +
  geom_histogram(fill = "green", color = "#000000", position = "identity") +
  ggtitle("lasso-imputed distribution") +
  theme_classic()
h6 <- ggplot(dfSPACE_imputed_comb, aes(x = imputed_rf)) +
  geom_histogram(fill = "blue", color = "#000000", position = "identity") +
  ggtitle("rf-imputed distribution") +
  theme_classic()
plot_grid(h1, h2, h3, h4, h5, h6, nrow = 2, ncol = 3)

#Pick your method:
#MissForest
dfSPACE_imputed_mF <- data.frame(
  data1[,1:2],missForest(data1[,3:8])$ximp,missForest(data1[,11:16])$ximp
)
head(dfSPACE_imputed_mF)

#RandomForest
dfSPACE_imputed_rf <- data.frame(
  data1[,1:2],
  complete(mice(data1[,3:8], method = "rf")),
  complete(mice(data1[,11:16], method = "rf"))
)
head(dfSPACE_imputed_rf)

#Lasso-imputed
dfSPACE_imputed_lasso <- data.frame(
  data1[,1:2],
  complete(mice(data1[,3:8], method = "lasso.norm")),
  complete(mice(data1[,11:16], method = "lasso.norm"))
)
head(dfSPACE_imputed_lasso)

#chose an imputation method below
imputed <- dfSPACE_imputed_mF
imputed <- dfSPACE_imputed_lasso
imputed <- dfSPACE_imputed_rf

#This replaces NAs with imputed data points if each treatment has max 1 NAs all other data wont be affected
for (i in c(3:8)){
  for(j in c(1:nrow(data1))){
   data1[j,i] = ifelse(is.na(data1[j,i]),ifelse(data1$NA.8WG16[j] <=1 & data1$NA.8WG16_D[j] <=1, dfSPACE_imputed[j,i],data1[j,i]),data1[j,i])
  }
}
for (i in c(11:16)){
  for(j in c(1:nrow(data1))){
    data1[j,i] = ifelse(is.na(data1[j,i]),ifelse(data1$NA.PP2[j] <=1 & data1$NA.PP2_D[j] <=1, dfSPACE_imputed[j,i],data1[j,i]),data1[j,i])
  }
}

#This replaces NAs with imputed data points if there is a max of 3 NAs per protein all other data wont be affected
for (i in c(3:8)){
  for(j in c(1:nrow(data1))){
    data1[j,i] = ifelse(is.na(data1[j,i]),ifelse(data1$NA.8WG16[j] + data1$NA.8WG16_D[j] <=3, dfSPACE_imputed[j,i],data1[j,i]),data1[j,i])
  }
}
for (i in c(11:16)){
  for(j in c(1:nrow(data1))){
    data1[j,i] = ifelse(is.na(data1[j,i]),ifelse(data1$NA.PP2[j] + data1$NA.PP2_D[j] <=3, dfSPACE_imputed[j,i],data1[j,i]),data1[j,i])
  }
}




#Use test d.f to check PCA
test <- imputed[,2:14]

