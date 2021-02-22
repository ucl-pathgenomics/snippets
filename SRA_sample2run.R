# example of rading SRA sample accession data, look up the SRA Run accession, perhaps for a gneomics pipeline, and return a relational table.

get_runs = function(acc = "SRS049712"){
  # takes an SRA sample, returns a character vector of SRA runs
  command = paste0("wget -qO- ", "'","http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=", acc,"'", " | grep ", acc, " | cut -f1 -d")
  command = paste0(command, '","')
  a = system(command,intern = T)
  return(a)
}

library(stringr)
library(dplyr)
library(tidyr)

# read in github table
system("wget https://github.com/SPAAM-community/AncientMetagenomeDir/blob/master/ancientmetagenome-hostassociated/ancientmetagenome-hostassociated.tsv")
df = read.table("ancientmetagenome-hostassociated.tsv",header = T,sep = "\t")
df = df[df$archive == "SRA",]

# comma separated rows -> rows
df = df %>% 
  mutate(archive_accession = strsplit(as.character(archive_accession), ",")) %>% 
  unnest(archive_accession)


df2 = data.frame()
for(i in 1: nrow(df)){
  sample = df$archive_accession[i]
  
  runs = get_runs(sample)
  for(run in runs){
    df2 = rbind(df2, data.frame(sample, runs))
    
  }
}


df3 = merge(df,df2, by.x = "archive_accession", by.y = "sample")

write.csv(df3, row.names = F, "sra_data.csv")
write.csv(df3$runs, row.names = F, "sra_list2.tab", col.names = NULL, quote = F)

