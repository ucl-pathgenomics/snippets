# calculates nucleodite diversity measures for an alignment

### changeable
alignment = "results/figure1/all_add_best_nodup.fasta"
### end 


library(Bios2cor)
library(ape)


hzyg <- function(i, seq){
  # .6 secs per 1000
  # heterozygocity
  t = seq[,i]
  #n = length(t)
  t = as.data.frame(table(t))
  t = t[t$t != "n",] #remove n sequences
  t = t[t$t != "-",] #remove indel sequences, hahn book says.
  n = sum(t$Freq) # updated n only non-ambiguous bses
  t.nuc = sum((t$Freq / n)^2)
  t.hzyg = (n/(n-1))*(1-t.nuc)
  return(t.hzyg)
}





mds = read.dna(alignment,format = "fasta",as.matrix = T,as.character = T) #not add
# may want to run this on add


for(pos in 1:length(mds[1,])){ # for each pos
  
  #----------- shannon entropy
  align = as.list(mds[,pos])
  entropy = Bios2cor::entropy(align)
  
  
  
  #----------- heterozygocty
  t.hzyg = hzyg(pos, mds)
  
  
  
  if(pos == 1){cat(paste("pos,entropy,hzyg", alignment),sep = "\n",file = "results/figure1/nucdiv-add.csv")}
  cat(paste(pos, entropy,t.hzyg, sep = ","),sep = "\n",file = "results/figure1/nucdiv-add.csv",append = T)
  
}
