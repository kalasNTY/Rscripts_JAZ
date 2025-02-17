bedtofasta <- function(bedtable, output = "my.fasta") {
  setwd("~/Desktop/Argh")
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(rtracklayer)
  
  humbug <- import(bedtable, format="bed")

  seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, humbug)
  #names(seq) = paste0("SEQUENCE_", seq_along(seq))
  names(seq) = humbug$name
  Biostrings::writeXStringSet(seq, output)
  
  
}
