asoiter <- function(file,asolen =25) {

library(Biostrings)
  
seq <-readChar(file, file.info(file)$size)
seq <- gsub("\n","",seq)
seq <- DNAString(seq)
seq <- reverseComplement(seq)
seq <- as.character(seq)
length <- nchar(seq)
nrmer <-length-asolen+1
k <- 1
n <- 1


results <- matrix(c(1:nrmer),nrmer,1)
dim(results) <-  c(nrmer,1)
  for(n in 1:nrmer){
    mer <- substr(seq,k,k+24)
    results[k,1] <- mer
    k <- k+1
  }
results <- data.frame(I(results))
gccont <- matrix(c(1:nrmer),nrmer,1)
gccont <- as.data.frame(gccont)

k <- 1
  for(n in 1:nrmer){
    aso <-results[k,1]
    asog <- gsub("G","",aso)
    asoc <- gsub("C","",aso)
    nrg <<- asolen-nchar(asog)
    nrc <<- asolen-nchar(asoc)
    gcs <- 100*((nrg+nrc)/asolen)
    gccont[k,1] <- gcs
    k <- k +1  
  }
  
  results <- cbind(results,gccont)

  gstretch <- matrix(c(1:nrmer),nrmer,1)
  gstretch <- data.frame(gstretch)
  
  k <- 1
  for(n in 1:nrmer){
    aso <-results[k,1]
    asogggg <- gsub("GGGG","",aso)
    gstr <- asolen-nchar(asogggg)
    if(gstr>0){
      gstr <- gstr/4
    }
    gstretch[k,1] <- gstr
    k <- k +1  
  }  

  results <- cbind(results,gstretch)
  colnames(results) <- c("Sequence","GCs","G_stretches")
  print (results)
  filtered <- results[which(results$GCs >= 40 &results$GCs <= 60 & results$G_stretches ==0),]
  
  write.csv(filtered,paste("ASOs_",gsub(".txt","",file)))

}
