###
# Author: Charles OJ
# Date April 2020
# What: Calculates FST for a number of sequence files & prints FST & bootstrapped p-values
# calculats FST by heterozygosity, uses hahn (2019) eq 3.2, eq 5.9 for multi locus Gst Estimate of nei 1979
###
# input - 1 or many sequence alignment files in fasa format
# output - csv printed results


##options
opt.remove_n = T
opt.bootstraps = 1000
gene_chunk = 1:173 #as slow, make many streams





library(ape)
library(stringr)


# actual
hzyg = function(i, seq){
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

fin.out = data.frame(gene = "", u95 = 1, l95 = 1, u99 = 1, l99 = 1)[-1,]

pop = as.factor(c(rep("af",32), rep("eu",32))) # assign pop
genes = list.files("data/geo_cluster/af_gene")
genes = tools::file_path_sans_ext(genes)
fst.dat = data.frame(gene = "", gamma_st = 0.0001, signif = "")[-1,]
for(gene in genes){
#for(gene in genes[gene_chunk]){
  #gene = genes[1]
  #gene = "test"
  
  ### af+eu
  class = "eu+af"
  infile = paste("data/geo_cluster/",class,"_gene/",gene,".fasta", sep = "")
  seq = read.dna(infile, format = "fasta", as.matrix = T, as.character = T)
  
  n_all = length(seq[,1])
  n_bases = length(seq[1,])
  
  # # remove n letters - may affect results
  cols = colSums(seq == "n")
  cols2 = cols = colSums(seq == "-")
  cols = cols + cols2 # count per position, number of n or - characters
  cols = ifelse(cols > 4,F,T)  #//change
  t.cols = as.data.frame(table(cols))
  t.cols = t.cols[t.cols$cols == F,]
  if(nrow(t.cols) > 0){
    o.frac_n = t.cols$Freq[1] / n_bases
  }else{
    o.frac_n = 0
  }
  
  if(opt.remove_n == T){
    seq = seq[,cols]
  }
  
  if(o.frac_n > 0.8){  ##hardcoded alter
    
    out = paste(gene, "NA", "too many n bases", "l99", "l95", "u95", "u99", o.frac_n, sep = ", ")
    cat(out, sep = "\n")
    
   next 
  }
  
  dat_all = lapply(1:length(seq[1,]), hzyg, seq=seq)
  pi_t = sum(unlist(dat_all))/n_bases
  
  
  seq_eu = seq[pop == "eu",]
  seq_af = seq[pop == "af",]
  
  dat_eu = lapply(1:length(seq[1,]), hzyg, seq=seq_eu)
  pi_eu = unlist(dat_eu)
  pi_eu[is.na(pi_eu)] = 0
  pi_eu = sum(pi_eu)/n_bases
  dat_af = lapply(1:length(seq[1,]), hzyg, seq=seq_af)
  pi_af = unlist(dat_af)
  pi_af[is.na(pi_af)] = 0
  pi_af = sum(pi_af)/n_bases
  
  pi_s = (pi_eu + pi_af)/2
  
  gamma_st = (pi_t - pi_s)/ pi_t # the Genome sequence estimate of Gst
  #print(gamma_st)
  
  
  ##### bootstrapping ######
  boot_dat = data.frame(gene = "gene", rep = 0, boot_fst = 0 )[-1,]
  for(rep in 1:opt.bootstraps){
    # randomly assign sequences
    vec_all = 1:n_all
    ran_eu = sample(1:n_all, n_all/2, replace=FALSE) #should be the right algorithm
    ran_af = vec_all[! vec_all %in% ran_eu]
    
    boot_eu_seq = seq[ran_eu,]
    boot_af_seq = seq[ran_af,]
    
    #debug
    #boot_eu_seq = seq[pop == "eu",]
    #boot_af_seq = seq[pop == "af",]
    
    
    boot_dat_eu = lapply(1:length(boot_eu_seq[1,]), hzyg, seq=boot_eu_seq)
    boot_pi_eu = unlist(boot_dat_eu)
    boot_pi_eu[is.na(boot_pi_eu)] = 0
    boot_pi_eu = sum(boot_pi_eu)/n_bases
    boot_dat_af = lapply(1:length(boot_af_seq[1,]), hzyg, seq=boot_af_seq)
    boot_pi_af = unlist(boot_dat_af)
    boot_pi_af[is.na(boot_pi_af)] = 0
    boot_pi_af = sum(boot_pi_af)/n_bases
    
    boot_pi_s = (boot_pi_eu + boot_pi_af)/2
    
    boot_gamma_st = (pi_t - boot_pi_s)/ pi_t # the Genome sequence estimate of Gst
    
    
    t.boot_dat = data.frame(gene = gene, rep = rep, boot_fst = boot_gamma_st )
    
    
    ##### end bootstrap #####
    
    
    boot_dat = rbind(boot_dat, t.boot_dat) # lots of stats per pos, per rep
    #print(rep)
    
  }
  #have checked this thoroughly so happy not to paste
  #write.csv(boot_dat, paste("data/geo_cluster/fst/", gene, "_boot_",opt.bootstraps,"_rep.csv", sep = "")) #do want it here as a sanity check for running
  #write.csv(boot_dat, paste("data/geo_cluster/fst_non/", gene, "_boot_",opt.bootstraps,"_rep.csv", sep = "")) #do want it here as a sanity check for running
  
  t1 = boot_dat[order(boot_dat$boot_fst),] # order by fst lowest at top
  l99 = t1$boot_fst[.01 * opt.bootstraps]
  u99 = t1$boot_fst[0.99 * opt.bootstraps]
  l95 = t1$boot_fst[0.05 * opt.bootstraps]
  u95 = t1$boot_fst[0.95 * opt.bootstraps]
  
  
  
  if(gamma_st < l99){t.signif = "l99"
  }else if(gamma_st < l95){t.signif = "l95"
  }else if(gamma_st > u99){t.signif = "u99"
  }else if(gamma_st > u95){t.signif = "u95"
  }else{t.signif= "not"}
  
  
  
  
  
  t.fst = data.frame(gene = gene, gamma_st = gamma_st, signif = t.signif, l99 = l99, l95 = l95, u99 = u99, u95 = u95)
  #fst.dat = rbind(fst.dat, t.fst)
  #print(gene)
  out = paste(gene, gamma_st, t.signif, l99, l95, u95, u99, o.frac_n, sep = ", ")
  cat(out, sep = "\n") # cat is clean, print add's extra text
  t.signif = "" # just a sanity check to refresh it
  
}


