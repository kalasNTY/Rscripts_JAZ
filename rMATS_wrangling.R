#Read substratified SE rMATS and separate significant events from background to use in STREME analysis
setwd("~/Downloads/NEW/")
df = read.delim("~/Downloads/NEW/PTB_vs_MATR_KD_comparisons/KD_co-Reg_edit.txt",stringsAsFactors=FALSE)
head(df)
#Define background splice events as FDR > 0.1 and |dPSI| < 0.05
bg = df[df$FDR>0.1&df$IncLevelDifference<abs(0.05),]

#Now get the upregulated and downregulated sets

up = df[df$FDR<0.1&df$IncLevelDifference>0.05,]
down =df[df$FDR<0.1&df$IncLevelDifference<(-0.05),]

#Now prepare BED files of each set
colnames(up)
up_bed = up[,c(4,6,7,1,2,5)]
up_bed$GeneID <- 0

down_bed = down[,c(4,6,7,1,2,5)]
down_bed$GeneID <- 0

bg_bed = bg[,c(4,6,7,1,2,5)]
bg_bed$GeneID <- 0

nrow(up_bed)
nrow(down_bed)
nrow(bg_bed)

#Downsample background to n = 1500 
bg1500_bed = bg_bed[sample(nrow(bg_bed),1500),]

#Create new beds for intronic regions flanking the exon of interest in a given window size up/downstream
#Take into consideration the strandedness upstream on + equals downstream on - in bed format
window = 500

#Function to find upstream start
ups = function(x,output){
  upstart = 0
  estart = as.numeric(x[2])
  eend = as.numeric(x[3])
  strand = x[6]
  if (strand == "+") {
    upstart = estart - 1 - window
  } else {
    upstart = eend + 1
  }
  return(upstart)
}

#Function to find upstream end
upe = function(x,out){
  upend = 0
  estart = as.numeric(x[2])
  eend = as.numeric(x[3])
  strand = x[6]
  if (strand == "+") {
    upend = estart - 1
  } else {
    upend = eend + 1 + window
  }
  return(upend)
}

#Function to find downstream start
ds = function(x,ou){
  dstart = 0
  estart = as.numeric(x[2])
  eend = as.numeric(x[3])
  strand = x[6]
  if (strand == "+") {
    dstart = eend + 1
  } else {
    dstart = estart - 1 - window
  }
  return(dstart)
}

#Function to find downstream end
de = function(x,outp){
  dend = 0
  estart = as.numeric(x[2])
  eend = as.numeric(x[3])
  strand = x[6]
  if (strand == "+") {
    dend = eend + 1 + window
  } else {
    dend = estart - 1
  }
  return(dend)
}

#Calculate upstream and downstream windows using above functions for all beds
bg1500_bed$upstart = apply(bg1500_bed, 1, ups)
bg1500_bed$upend = apply(bg1500_bed, 1, upe)
bg1500_bed$dstart = apply(bg1500_bed, 1, ds)
bg1500_bed$dend = apply(bg1500_bed, 1, de)

head(bg1500_bed)

up_bed$upstart = apply(up_bed, 1, ups)
up_bed$upend = apply(up_bed, 1, upe)
up_bed$dstart = apply(up_bed, 1, ds)
up_bed$dend = apply(up_bed, 1, de)

head(up_bed)

down_bed$upstart = apply(down_bed, 1, ups)
down_bed$upend = apply(down_bed, 1, upe)
down_bed$dstart = apply(down_bed, 1, ds)
down_bed$dend = apply(down_bed, 1, de)

head(down_bed)

#Reorder/delete columns and export beds
dPSIdown_upstream = down_bed[,c(1,7,8,4,5,6)]
dPSIdown_downstream = down_bed[,c(1,9,10,4,5,6)]
dPSIup_upstream = up_bed[,c(1,7,8,4,5,6)]
dPSIup_downstream = up_bed[,c(1,9,10,4,5,6)]
bg_upstream = bg1500_bed[,c(1,7,8,4,5,6)]
bg_downstream = bg1500_bed[,c(1,9,10,4,5,6)]
head(dPSIup_downstream)


write.table(dPSIdown_upstream,"dPSIdown_upstream.bed",row.names = FALSE,col.names = FALSE, sep = "\t",quote = FALSE)
write.table(dPSIdown_downstream,"dPSIdown_downstream.bed",row.names = FALSE,col.names = FALSE, sep = "\t",quote = FALSE)
write.table(dPSIup_upstream,"dPSIup_upstream.bed",row.names = FALSE,col.names = FALSE, sep = "\t",quote = FALSE)
write.table(dPSIup_downstream,"dPSIup_downstream.bed",row.names = FALSE,col.names = FALSE, sep = "\t",quote = FALSE)
write.table(bg_upstream,"bg_upstream.bed",row.names = FALSE,col.names = FALSE, sep = "\t",quote = FALSE)
write.table(bg_downstream,"bg_downstream.bed",row.names = FALSE,col.names = FALSE, sep = "\t",quote = FALSE)


#upload beds to CREATE/NEMO and use bedtools getfasta to create primary sequence list from genome of choice (streme_bed.sh)
#bedtools getfasta -fi ~/Documents/hcrsite-local-shared/uploaded/media/trxpt/GRCh38.primary_assembly.genome.fa -fo streme_fasta2 -bed test.bed -s


#next use meme suite streme to carry out de novo motif discovery. 

#streme --p ./streme_fasta/dPSIdown_downstream.bed.fasta --n ./streme_fasta/bg_downstream.bed.fasta --rna --nmotifs 100 --o streme_output/dPSIdown_downstream --verbosity 1

#Carry on using sequences and motifs and input them into centrimo to find positional enrichment of discovered motis
#centrimo --neg ./streme_fasta/bg_upstream.bed.fasta --o ./centrimo_out --local --verbosity 1 ./streme_fasta/dPSIup_upstream.bed.fasta ./streme_output/dPSIup_upstream/streme.txt




