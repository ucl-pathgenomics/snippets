# script to optimise multiple sequence alignments, by their translated protein sequences. For a whole genome.

library(ape)
library(Biostrings)       # Provides DNAString, DNAStringSet, etc
library(GenomicRanges)    # Provides GRanges, etc
library(GenomicFeatures)  # provides Txdb
library(stringr)
library(DECIPHER) # amino acid alignment handling

### inputs
in.msa = "all_raw_best_msa_man.fasta"
in.ref.match = "erlin"    # reference should be aligned, this is a unique string to identify the sequence in the msa.
in.ref_seq = "ref/NC_006273.2.fasta"
in.ref_gff = "ref/NC_006273.2.gff3"
### inputs end

### load
msa = ape::read.dna(in.msa, format = "fasta",as.matrix = T)
ref.seq = as.character(ape::read.dna(in.ref_seq, "fasta"))
ref.seq.string = paste(ref.seq, collapse = "")
### end load

# get transcripts
txdb <- makeTxDbFromGFF(file=in.ref_gff, format="gff3")
#gn <- genes(txdb)
gn = transcriptsBy(txdb, by = "gene") # as we
gn = as.data.frame(gn)
gn = gn[order(gn$start),]
gn = gn[!grepl("rna", gn$tx_name),] # remove rna


# check transcripts do not overlap
for(i in 2:length(gn)){ 
  if(gn[i,4] < gn[i-1,5]){
    print(paste(i, "error"))
    warning("transripts are overlapping")
    break
  }
}

pattern = "[^ACTG-]"


#for each transcript
# get reference start and end kmer
# locate start and end of reference in msa
# extract msa chunk
# RC if needed
# align by transcript & reverse back
# bind by row to part before. using rbind of matrix
# next

ref.msa.index = grep(pattern = in.ref.match,labels(msa)) # index of ref in msa
ref.msa.string = as.character(msa[ref.msa.index,])
ref.msa.string = paste(ref.msa.string, collapse = "")
for( t1 in 1:length(gn)){ # for each transcript
  print(paste(gn$tx_name[t1], gn$tx_id[t1]))
  
  #------------------------------ get reference start and end kmer
  t.start = gn$start[t1]
  t.end = gn$end[t1] + (gn$width[t1] -1)
  t.ref = ref.seq[t.start:t.end]
  t.ref.collapse = paste(t.ref,collapse = "") # seq of transcript
  
  # buld regex allowing for insertions for start and end of transcript
  t.pattern = ""  
  for(t.char in t.ref){
    t.pattern = paste0(t.pattern, t.char, "[-]{0,1000}") # allow multiple insertions
  }
  
  # ----------------------------- locate start and end of reference in msa
  match.loc = str_locate(ref.msa.string, t.pattern) 
  match.start = match.loc[1]
  match.end = match.loc[2]
  
  
  # ----------------------------- extract msa regions
  temp.msa = msa[,match.start:match.end]
  
  
  # ----------------------------- RC if needed
  if(gn$strand[t1] == "-"){
    # reverse complement
    temp.msa = ape::complement(temp.msa)
  }else{
    # do nothing
  }
  
  #write.FASTA(temp.msa,file = "temp.msa.fasta")
  
  
  #---------------------------- insert between transcript blocks
  
  if(t1 == 1){ # if the first transcript block
    new.msa = msa[,1:(match.start - 1)] # all the alignment up until the first transcript
  }else{
    # take last transcript + 1 till new transcript - 1
    new.msa = as.character(cbind(new.msa, msa[, (match.end.prev + 1) : (match.start - 1)])) # append msa between last transcript and current transcript
  }
  match.end.prev = match.end # update for next loop
  
  
  
  
  
  # ----------------------------- align by transcript & reverse back
  # as this is all definitely now a transcript block we can remove all the gaps before re-alignment
  # http://www.bioconductor.org/packages/release/bioc/vignettes/DECIPHER/inst/doc/ArtOfAlignmentInR.pdf
  a = del.gaps(temp.msa)
  a = a %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet # format conversion
  
  if(min(a@ranges@width) < 3){ # if we have issues with sequences
    new.msa = as.character(cbind(new.msa, temp.msa)) # just append as is.
    next
  }
  

  a = DECIPHER::AlignTranslation(a) # translate to AA, align, then translate back 
  write.dna(a,format = "fasta", "temp_a.fasta")
  b = ape::read.dna("temp_a.fasta", format = "fasta",as.matrix = T) # bodge, cant get this to work well without
  
  
  # -----------------------------  bind by col to part before.
  

  
  #write.dna(new.msa,format = "fasta", "temp_new.msa.fasta")
  new.msa = as.character(cbind(new.msa, b)) # all new columns to existing new 
  #write.dna(new.msa,format = "fasta", "temp_new.new.msa.fasta")
  
  if(t1 == length(gn)){ # if this is the last transcript
    # append rest of msa
    new.msa = as.character(cbind(new.msa, msa[, (match.end + 1) : length(msa[1,])]))
  }
  
}

write.FASTA(new.msa, file = "allignment_optimised_msa.fasta")

print("complete")










