library(biomaRt)
library(Biostrings)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)


setwd("~/Downloads/")
filename = "3utr"
fasta = read.csv2(filename,sep = "")
colnames(fasta)
head(fasta)
str(fasta)
dim(fasta)
fasta[2,2]
sequences = fasta[,2]
fastaFile <- RNAStringSet(x=sequences, use.names=TRUE)


#load fastafile to be analysed
rm(fastaFile)

filepath = "hairpin.fa"
fastaFile <- readRNAStringSet(filepath)
#fastaFile <- DNAStringSet(fastaFile)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
gencode32 <- data.frame(seq_name, sequence)
dim(gencode32)

#use targetscan data.Remove ensemblID version from ID field
ts7.2 = read.delim("Predicted_Targets_Info.default_predictions.txt")
ts7.2 <-  ts7.2[ts7.2$Species.ID==9606,]
ts7.2$Gene.ID <- substr(ts7.2$Gene.ID,1,15)
ts7.2$Transcript.ID <- substr(ts7.2$Transcript.ID,1,15)

summary(seq_name %in% ts7.2$Gene.ID)
res <- ts7.2[ts7.2$Gene.ID %in% seq_name,]

#filter results for seed type eg "7mer-a1"   "7mer-m8"   "8mer"      "cent.4.14" "cent.4.15" "cent.5.15"
eightmer <- res[res$Seed.match=="8mer",]
sevenmer <- res[which(res$Seed.match == "7mer-a1" | res$Seed.match == "7mer-m8"),]

miRsites <- sevenmer %>%
  group_by(miR.Family) %>%
  summarise(count = n_distinct(Transcript.ID))

miRsites <- miRsites[order(-miRsites$count),]
head(miRsites,20)
hist(miRsites$count,main=filepath,xlab="miR target sites", prob = T)



#save different results to compare
un.8mer <- miRsites
un.7mer <- miRsites
up.8mer <- miRsites
up.7mer <- miRsites
down.8mer <- miRsites
down.7mer <- miRsites


colnames(un.8mer) <- c("mir","count")
colnames(up.8mer) <- c("mir","count")
colnames(down.8mer) <- c("mir","count")
colnames(un.7mer) <- c("mir","count")
colnames(up.7mer) <- c("mir","count")
colnames(down.7mer) <- c("mir","count")

#**********************WAYPOINT 1**************************
rm(miRsites,eightmer,sevenmer,res)


#merge datasets and replace empty points (NaN) with either 0 or 0.1 for normalisation
all <- merge(up.8mer,down.8mer,by = "mir",all=T)
all <- merge(all,un.8mer,by = "mir",all=T)
all[is.na(all$count),19] <- 0.1
all[is.na(all)] <- 0

all$up_norm <- all$count.x/all$count
all$down_norm <- all$count.y/all$count
all$FC <- ifelse(is.infinite(all$down_norm/all$up_norm),100,all$down_norm/all$up_norm)
all[is.na(all$FC),7] <- 1

#replace inf values with 2x the highest finite value as placeholder
#all[is.infinite(all$up_norm),5] <- 2*max(all[is.finite(all$up_norm),5])
#all[is.infinite(all$down_norm),6] <- 2*max(all[is.finite(all$down_norm),6])

all <- all[order(-all$FC),]

#optionally add further motif (information) to dataframe all. Eg is miR contains fbe. SAVE IN df.m fist
all$miR <- sub("-5p","",all$mir)
all$miR <- sub("-3p","",all$miR)
#all$miR <- paste0(substr(all$miR,1,4),gsub('-.$', '', substr(all$miR,5,10))) #this removes family members (-1,-2) after miR name

#filter miRNAs

all <- all[which(nchar(substr(all$miR,5,10))<4),]

summary(df.m$miR %in% all$miR)
all <- all[all$miR %in% df.m$miR,]
all <-  merge(df.m[, c("miR","motifs","motifs per kb")], all, by="miR",all=T)
all <-  all[!is.na(all$mir),]
all[is.na(all$motifs),2:3] <- 0
all$motifs <- as.integer(all$motifs)

#normalised to unchanged miRsite occurence use box-cox transformation
ggplot(all,aes(x=log1p(up_norm),y=log1p(down_norm),color=motifs))+
#  xlim(0,7) + ylim(0,7)+
  geom_point(size=1)+
  #stat_smooth(method = "lm",col = "#C42126",se = T,size = 1)+
  theme_minimal()+
  ggtitle("Enriched 8mer target sites sRNA-seq siFOX2") +
  xlab("log1p up regulated mRNAs (norm)") + ylab("log1p down regulated mRNAs (norm)") +
 geom_text(size=1,label=substr(all$mir,5,20), nudge_x = 0.05, nudge_y = 0.05, check_overlap = T) +
  stat_smooth(method = "lm",col = "#C42126",se = T,size = 1) +
#  geom_line(test, aes(x= up_norm,y=down_norm,color= "blue")) +
  scale_color_viridis(option = "C")

#not normalised 
ggplot(all,aes(x=MPKB.x/trans.x,y=MPKB.y/trans.y,color=motifs))+
  #xlim(0,2.5) + ylim(0,8)+
  #geom_point()+
  stat_smooth(method = "lm",col = "#C42126",se = T,size = 1)+
  theme_minimal()+
  ggtitle("Overrepresented miRNA sites") +
  xlab("Occurence in upregulated set") + ylab("Occurence in downregulated set") +
  geom_text(size=2,label=substr(all$mir,5,20), nudge_x = 0.1, nudge_y = 0.1, check_overlap = T) +
  scale_color_viridis(option = "C")

#get genesymbols from mart 
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("ensembl_gene_id","hgnc_symbol")
filters <- c("ensembl_gene_id")
values <- seq_name
lookup <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
#lookup$hgnc_symbol %in% ts7.2$Gene
#res <- ts7.2[ts7.2$Gene %in% lookup$hgnc_symbol,]


celf <- RNAString("UGUUU")
are <- RNAString("AUUUA")
fbe <- RNAString("gcaug")
fbe <- RNAString("GNNYG")
mbnl <- RNAString("YGCY")
mbnl <- RNAString("YYYGCYKY")

mir32_5p <- RNAString("auugcaca")
mir32_3p <- RNAString("AAUUUAGU")

mir20b_5p <- RNAString("GAGAACCA")
mir20b_3p <- RNAString("CAGAACAA")

dna.fbe <- DNAString("GNNYG")
dna.canon <- DNAString("tgcatgt")
dna.mbnl <- DNAString("YGCY")
dna.hur <- DNAString("ATTTA")

mir32_5p <- DNAString("attgcaca")
mir32_3p <- DNAString("AATTTAGT")

mir32_5p_target <- reverseComplement(mir32_5p)
mir32_3p_target <-  reverseComplement(mir32_3p)
mir20_5p_target <- reverseComplement(mir20b_5p)
mir20_3p_target <-  reverseComplement(mir20b_3p)

len5p <- mir20_5p_target@length
len3p <- mir20_3p_target@length


m0 <- vmatchPattern(mir20_5p_target,fastaFile,fixed="subject")
m2 <- vcountPattern(mir20_5p_target,fastaFile,fixed="subject")
m1 <- fastaFile@ranges@width
m3 <- 1000*(m2/m1)
m4 <- fastaFile@ranges@NAMES
df.m <- as.data.frame(cbind(m1,m2,m3))
colnames(df.m) <- c("length","motifs","motifs per kb")
rownames(df.m) <- m4
df.m <- df.m[order(-df.m$motifs),]


#optional filter criteria (e.g. <400bp length)
df.m <- df.m[which(df.m$length<400),]

#Remove strings from mirbase annotation
df.m$miR <- rownames(df.m)
rownames(df.m) <- c(1:nrow (df.m))
df.m <- df.m[,c(4,1,2,3)]
df.m$miR <- sub("Homo sapiens", "", df.m$miR)
df.m$miR <- sub("stem-loop", "", df.m$miR)
head(df.m)
#df.m <- df.m[which(df.m$motifs>0),]
df.m <- df.m[grep("hsa-",rownames(df.m)),]
df.m$miR  <- substr(df.m$miR,1,15)
df.m$miR <- gsub(" .*","",df.m$miR)
df.m$miR <- gsub("hsa-","",df.m$miR)
df.m$miR <- gsub("mir-","miR-",df.m$miR)

#Remove strings from mart annotation
df.m$GeneID <- rownames(df.m)
rownames(df.m) <- c(1:nrow (df.m))
df.m <- df.m[,c(4,1,2,3)]
df.m$ENSG  <- substr(df.m$GeneID,1,15)
df.m$ENST  <- substr(df.m$GeneID,33,53)
df.m$ENST <- sub(".*\\|(.+)\\|.*", "\\1", df.m$ENST)
df.m$GeneID <- df.m$ENSG 
head(df.m)



hist(df.m$`motifs per kb`,main=filepath,xlab="motifs",breaks=10,col='blue', xlim=c(0,5), prob = T)
lines(density(df.m$`motifs per kb`), col = "red", lwd = 2)
summary(df.m$`motifs per kb`)

#summarise dataset either by looking into ENST or ENSG
df.m <- df.m %>% group_by(GeneID) %>% summarise(length = median(length), all_motifs = sum(motifs), trans =n_distinct(ENST),  motifs_gene = all_motifs/trans, MPKB = median(`motifs per kb`),sdMPKB = sd(`motifs per kb`))
#df.m <- df.m %>% group_by(ENST) %>%  summarise(ENSG = ENSG, length = median(length), all_motifs = sum(motifs), trans =n_distinct(ENST),  motifs_gene = all_motifs/trans, MPKB = median(`motifs per kb`),sdMPKB = sd(`motifs per kb`))

fivep <- df.m
threep <- df.m

sites.5p <- m0
sites.3p <- m0 

#fix and filter fastaFile names from miRbase annotation
sites.3p <- sites.3p[grep("hsa-",sites.3p@NAMES)]
sites.3p@NAMES <- substr(sites.3p@NAMES,1,17)
sites.3p@NAMES <- gsub(" .*","",sites.3p@NAMES)
sites.3p@NAMES<- gsub("hsa-","",sites.3p@NAMES)
sites.5p <- sites.5p[grep("hsa-",sites.5p@NAMES)]
sites.5p@NAMES <- substr(sites.5p@NAMES,1,17)
sites.5p@NAMES <- gsub(" .*","",sites.5p@NAMES)
sites.5p@NAMES<- gsub("hsa-","",sites.5p@NAMES)
#fix and filter fastaFile names from mart annotation
sites.3p@NAMES <- substr(sites.3p@NAMES,1,15)
sites.5p@NAMES <- substr(sites.5p@NAMES,1,15)

#compare sites between fasta files (e.g. upregulated vs downregulated)
matup <- df.m$`motifs per kb`
matdown <- df.m$`motifs per kb`
matunch <- df.m$`motifs per kb`

t.test(matdown,matup)

#find all genes that contain both motifs at the same time from mirbase annotation
fivep <-  fivep[which(fivep$motifs>0),]
threep <-  threep[which(threep$motifs>0),]
t1 <- threep[rownames(fivep),]
t1 <- t1[which(!is.na(t1$motifs)),]
#geneID <- rownames(t1)
t1$miR <- sub("miR","mir",t1$miR)
geneID <- t1$miR

#find all genes that contain both motifs at the same time from mart annotation
fivep <-  fivep[which(fivep$motifs_gene>0),]
threep <-  threep[which(threep$motifs_gene>0),]
t1 <- merge(fivep,threep,by = "GeneID",all=F)
geneID <- t1$GeneID


#find gene postions where both motifs are close together
results <- vector(mode = "list", length = length(geneID))
sites.3p <- unlist(sites.3p)
sites.5p <- unlist(sites.5p)


#For mart annotation only: remove all fields except ENSG
sites.3p@NAMES <- substr(sites.3p@NAMES,1,15)
sites.5p@NAMES <- substr(sites.5p@NAMES,1,15)

#For mart annotation only: remove all fields except ENST
#sites.3p@NAMES <- substr(sites.3p@NAMES,33,53)
#sites.3p@NAMES <- sub(".*\\|(.+)\\|.*", "\\1", sites.3p@NAMES)
#sites.5p@NAMES <- substr(sites.5p@NAMES,33,53)
#sites.5p@NAMES <- sub(".*\\|(.+)\\|.*", "\\1", sites.5p@NAMES)


for (i in 1:length(geneID)) {
  
  store <- sites.3p[which(sites.3p@NAMES==geneID[i])]
  store2 <- sites.5p[which(sites.5p@NAMES==geneID[i])]
    
    dist <- outer(store@start,store2@start, "-")
    hit <- dist[abs(dist)<20]
    ind <- which(abs(dist)<20, arr.ind = T)
    pos3p <- as.data.frame(store@start[ind[,1]])
    ifelse(nrow(pos3p)>0, results[[i]] <- list(as.data.frame(pos3p)),results[[i]] <-NA)
    names(results[[i]]) <- geneID[i]
  
}

results <-  results[!sapply(results,is.na)]
genes  <-   sapply(results,names)

#get genomic coordinates
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("chromosome_name","start_position","end_position","ensembl_gene_id","strand","hgnc_symbol")
filters <- c("ensembl_gene_id")
values <- genes
all.mart <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
all.mart <- all.mart[order(all.mart$chromosome_name,all.mart$start_position),]

#translate position on ensemblID into genomic coordrinates
oldbed <- data.frame()
for (k in 1:length(results)) {
  rw <- grep(names(results[[k]]),all.mart[,4])
  chr <-all.mart[rw,1]
  start <-all.mart[rw,2]
  end <-all.mart[rw,3]
  strnd <-all.mart[rw,5]
  
  genpos <- sapply(results[[k]],function(x) start+x*strnd)
  newbed <- lapply(genpos, function(x) cbind(chr,x,x+len3p,strnd,names(results[[k]])))
  newbed <- as.data.frame(newbed[[1]])
  colnames(newbed) <- c("Chr","Start","End","Strand","GeneID")
  oldbed <- rbind(oldbed,newbed)
 
}
oldbed[,4]<- sapply(oldbed[,4], function(x) ifelse(x==1,x <- "+",x <- "-"))
oldbed <- makeGRangesFromDataFrame(oldbed,keep.extra.columns = T)
export.bed(oldbed,"5p_3p_hits.bed")

#all 3UTRs that have one or more motifs
select <- fasta[m2,]
dim(select)
#select 3'UTRs with many motifs
twenty <- m2[]>19
ten <- m2[]>9 & m2[]<20
five <- m2[]>4 & m2[]<10
one <- m2[]>0 & m2[]<5
none <- m2[]==0

n20 <-  fasta[twenty,]
n10 <-  fasta[ten,]
n5 <-  fasta[five,]
n1 <-  fasta[one,]
n0 <-  fasta[none,]

write(n0$X.name,"no32_5p_sites.txt")
write(n1$X.name,"four32_5p_sites.txt")
write(n5$X.name,"ten32_5p_sites.txt")

o20 <-  fasta[twenty,]
o10 <-  fasta[ten,]
o5 <-  fasta[five,]
o1 <-  fasta[one,]
o0 <-  fasta[none,]

p20 <-  fasta[twenty,]
p10 <-  fasta[ten,]
p5 <-  fasta[five,]
p1 <-  fasta[one,]
p0 <-  fasta[none,]

write(o0$X.name,"no32_3p_sites.txt")
write(o1$X.name,"four32_3p_sites.txt")
write(o5$X.name,"ten32_3p_sites.txt")

require("gplots")
require("VennDiagram")
vennlist = list (up = up.8mer , down = down.8mer, unchanged = un.8mer)

vennA <- venn(vennlist, show.plot=T)

#import tsv data (eg ENCODE) that have a mix of numeric and character columns
DE <-  read.csv2("ENCFF166HLT.tsv",sep = '\t')
DE[] <- lapply(DE, function(x) as.numeric(as.character(x)))
DE$GeneID <- rownames(DE)
rownames(DE) <-  c(1:nrow(DE))
str(DE)
signif <- DE[which(DE$padj<0.05),]
upregulated <- signif[which(signif$log2FoldChange>0.4),]
upregulated$GeneID <- substr(upregulated$GeneID,1,15)
downregulated <- signif[which(signif$log2FoldChange< -0.4),]
downregulated$GeneID <- substr(downregulated$GeneID,1,15)

notsignif <- DE[order(-DE$padj),][1:max(nrow(downregulated),nrow(upregulated)),]
notsignif$GeneID <- substr(notsignif$GeneID,1,15)
hist(signif$baseMean,breaks = 100,xlim=c(0,1000))


seq_name <- notsignif$GeneID
seq_name <- downregulated$GeneID
seq_name <- upregulated$GeneID

#***JULISCAN***
#select mature miRNAs from all fbe-pre-miRNAs take and find enrichment in gene set (JuliScan)
#first create fbe pre-mir list using above scripts save as data.frame gnuug, gcaug etc
filepath = "mature.fa"
fastaFile <- readRNAStringSet(filepath)
seq_name = names(fastaFile)


seeds <-  subseq(fastaFile, start = 2, width = 8)
seeds <- seeds[grep("hsa-",seeds@ranges@NAMES),]

seeds@ranges@NAMES <- gsub(" .*","",seeds@ranges@NAMES)
seeds@ranges@NAMES <- gsub("hsa-","",seeds@ranges@NAMES)

seeds <- reverseComplement(seeds)

filepath = "ENCFF166HLT_unchanged.txt"
fastaFile <- readDNAStringSet(filepath)
fastaFile <- fastaFile[which(seqlengths(fastaFile)>10)]
fastaFile <- RNAStringSet(fastaFile)
seq_name = names(fastaFile)

df.res <- data.frame()

for(s in 1:length(seeds)){
  
  s2 <- vcountPattern(seeds[[s]],fastaFile,fixed="subject")
  s4 <- ifelse(s2>0,TRUE,FALSE)
  s1 <- fastaFile@ranges@width
  s3 <- 1000*(s2/s1)
  s5 <- fastaFile@ranges@NAMES
  df.s <- as.data.frame(cbind(s5,s1,as.integer(s2),s3))
  colnames(df.s) <- c("GeneID","length","motifs","motifs per kb")
  df.s <- df.s[which(df.s$motifs>0),]
  df.s$mir <- c(rep(seeds@ranges@NAMES[s],nrow(df.s)))
  df.s$target <- c(rep(as.character(seeds[s]),nrow(df.s)))
  df.res <- rbind(df.res,df.s)
  }

df.res$ENSG  <- substr(df.res$GeneID,1,15)
df.res$ENST  <- substr(df.res$GeneID,33,53)
df.res$ENST <- sub(".*\\|(.+)\\|.*", "\\1", df.res$ENST)
df.res <-  df.res[,-1]
df.res <-  df.res[,c(6,7,1:5)]
df.res$length <- as.integer(df.res$length)
df.res$motifs <- as.integer(df.res$motifs)
df.res$`motifs per kb` <- as.numeric(df.res$`motifs per kb`)

df.res <- df.res[order(-df.res$motifs),]
df.res <- df.res[!is.na(df.res$ENSG),]

rm(s,df.s)

#old metrics
miRsites <- df.res %>%
  group_by(mir) %>%
  summarise(count = n_distinct(GeneID))

#new metrics
miRsites <- df.res %>%
  group_by(mir,target) %>% 
  summarise(count = sum(motifs), trans =n_distinct(ENST), gene=n_distinct(ENSG),  mpG = count/gene, MPKB = median(`motifs per kb`),sdMPKB = sd(`motifs per kb`))

#collapse identical 8mer seeds into one miRNA
miRsites <- miRsites %>%
  group_by(target) %>% mutate(count_family = sum(count)) %>% distinct(target, .keep_all=T)

miRsites <- miRsites[order(-miRsites$count_family),]

#run above on all data sets (up down unchanged)
un.8mer <- miRsites
up.8mer <- miRsites
down.8mer <- miRsites

#Continue to WAYPOINT1



#optional filter criteria (e.g. <400bp length)
df.m <- df.m[which(df.m$length<400),]

df.m <- df.m[which(df.m$motifs>0),]
df.m <- df.m[grep("hsa-",rownames(df.m)),]
df.m$miR <- rownames(df.m)
rownames(df.m) <- c(1:nrow (df.m))
df.m <- df.m[,c(4,1,2,3)]
df.m$miR <- sub("Homo sapiens", "", df.m$miR)
df.m$miR <- sub("stem-loop", "", df.m$miR)

#Remove strings from mirbase annotation
df.m$miR  <- substr(df.m$miR,1,17)
df.m$miR <- gsub(" .*","",df.m$miR)
df.m$miR <- gsub("hsa-","",df.m$miR)

hist(df.m$`motifs per kb`,main=filepath,xlab="motif per kb",breaks=1000,col='blue', xlim=c(0,1), prob = T)
lines(density(df.m$`motifs per kb`), col = "red", lwd = 2)
summary(df.m$`motifs per kb`)

fivep <- df.m
threep <- df.m

sites.3p <- m0
sites.5p <- m0

#compare sites between fasta files (e.g. upregulated vs downregulated)
matup <- df.m$`motifs per kb`
matdown <- df.m$`motifs per kb`
matunch <- df.m$`motifs per kb`

t.test(matdown,matup)

#fasta biostring to dataframe
dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))
seeds <- dss2df(test)


#BACKUP
#find gene postions where both motifs are close together
results <- vector(mode = "list", length = length(geneID))
for (i in 1:length(geneID)) {
  
  store <- sites.3p[[which(sites.3p@NAMES==geneID[i])]]
  store2 <- sites.5p[[which(sites.5p@NAMES==geneID[i])]]
  
  dist <- outer(store@start, store2@start, "-")
  hit <- dist[abs(dist)<10]
  ind <- which(abs(dist)<10, arr.ind = T)
  pos3p <- as.data.frame(store@start[ind[,1]])
  ifelse(nrow(pos3p)>0, results[[i]] <- list(as.data.frame(pos3p)),results[[i]] <-NA)
  names(results[[i]]) <- geneID[i]
  
}

results <-  results[!sapply(results,is.na)]
genes  <-   sapply(results,names)
