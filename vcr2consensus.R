# script to take a VCF and return a fasta consensus
# ignores indels

args = commandArgs(trailingOnly=TRUE)
library(stringr)

args = c("/archive/21_2_ancient_dna/SRR5581849/refbased/SRR5581849_allpos_variants.vcf",
         "/archive/21_2_ancient_dna/SRR5581849/refbased/SRR5581849_allpos_variants.vcf.fasta")

vcf = readLines(args[1])
t = grep("#CHROM", vcf)

vcf = read.table(args[1],sep = "\t", skip = t-1)
vcf$freq = str_extract(vcf[,10], "[0-9]{1,3}%")
vcf$freq = as.numeric(str_extract(vcf[,11], "[0-9]{1,3}"))
vcf = vcf[,c(2,4,5,11)]
names(vcf) = c("pos", "ref", "alt", "freq")



#build fasta
vcf2chrvec = function(i){
  if(i %in% vcf$pos){
    row = vcf[vcf$pos == i,]
    if(row$alt == "."){
      append = row$ref
    }else if(row$alt %in% c("A", "T", "C", "G") & row$freq > 50){
      append = row$alt
    }else{
      append = "-"
    }
  }else{
    append = "-"
  }
  append = as.character(append)
  return(append)
}
fasta.seq = sapply(1:max(vcf$pos), vcf2chrvec)
fasta.seq = paste0(fasta.seq, collapse = "")


accession = str_split(basename(args[1]),"_",simplify = T)[,1]
fasta.header = paste0(">", accession)
fasta.out = paste0(fasta.header, "\n", fasta.seq)
fasta.file = args[2]
cat(fasta.out,file = fasta.file)
