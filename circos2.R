# circular plot

# emulating Juanita code

#---------------data
# https://cran.r-project.org/web/packages/RCircos/vignettes/Using_RCircos.pdf


genestartstop = data.frame(Chromosome = "CMV_DNA", 
                           chromStart = gn@ranges@start,
                           chromEnd = gn@ranges@start + gn@ranges@width, 
                           GeneName = gn$gene_id)
genestartstop = genestartstop[!grepl("RNA|TRS|IRS",x = genestartstop$GeneName),] # remove RNA, repeats
write.table(genestartstop,file = "results/figure1/genestartstop.txt",  sep = "\t", row.names = F,quote = F)


dat = read.csv("results/figure1/nucdiv-add.csv")
 

hzyg = data.frame("CMV_DNA",
                  dat$pos,
                  dat$pos + 1,
                  1,
                  dat$hzyg)[-(235644:235646),] # remove last pos
hzyg<-hzyg[which(!is.na(hzyg$dat.hzyg)),] # remove na 
write.table(hzyg,file = "results/figure1/hzyg.txt",  sep = "\t", row.names = F,quote = F,col.names = F)

entropy = data.frame("CMV_DNA",
                  dat$pos,
                  dat$pos + 1,
                  1,
                  dat$entropy)[-(235644:235646),] # remove last pos
entropy<-entropy[which(!is.na(entropy$dat.entropy)),] # remove na 
write.table(entropy,file = "results/figure1/entropy.txt",  sep = "\t", row.names = F,quote = F,col.names = F)


kmer = read.csv("results/figure1/kmer_fvalues.csv",header = F)
kmer = data.frame("CMV_DNA",
                  kmer$V1,
                  kmer$V1+1,
                  1,
                  kmer$V2)
kmer<-kmer[which(!is.na(kmer$kmer.V2)),] # remove na 
write.table(kmer,file = "results/figure1/kmer.txt",  sep = "\t", row.names = F,quote = F,col.names = F)

#-------------- plot

library(RCircos)
data=read.table(file="results/figure1/data.txt",sep="\t",stringsAsFactors = FALSE,header = T)
data1=read.table(file="results/figure1/Circos_plot/BS_position.txt",sep="\t",stringsAsFactors = FALSE,header = F)
data2=read.table(file="results/figure1/genestartstop.txt",sep="\t",stringsAsFactors = FALSE,header=T)
data3=read.table(file="results/figure1/Circos_plot/BM_position.txt",sep="\t",stringsAsFactors = FALSE,header = F)
data4=read.table(file="results/figure1/hzyg.txt",sep="\t",stringsAsFactors = FALSE,header = F)
data5=read.table(file="results/figure1/entropy.txt",sep="\t",stringsAsFactors = FALSE,header = F)
data6=read.table(file="results/figure1/kmer.txt",sep="\t",stringsAsFactors = FALSE,header = F)

chr.exclude<-NULL
cyto.info<-data
tracks.inside<-6
tracks.outside<-1
RCircos.Set.Core.Components(cyto.info,chr.exclude,tracks.inside,tracks.outside)

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$base.per.unit <- 1
rcircos.params$track.background = "white"
rcircos.params$text.size = 4
rcircos.params$hist.width = 2 # histogram width
RCircos.Reset.Plot.Parameters(rcircos.params)

out.file<-"results/figure1/RCircos-plot.png"
#pdf(file=out.file,height=8,width=8,compress=TRUE)
png(file=out.file,height=10000,width=10000)
RCircos.Set.Plot.Area()

RCircos.Chromosome.Ideogram.Plot()

## genes will need to make custom function https://github.com/cran/RCircos/blob/master/R/RCircosPlotDataTracks.R
tile.data <- data2[,1:3]
tile.colors <- "blue";
tile.data["PlotColor"] <- tile.colors;
track.num <- 1;
RCircos.Tile.Plot(tile.data, track.num, "in");


name.col<-4
side<-"in"
track.num<-3
RCircos.Gene.Connector.Plot(data2,track.num,side)
track.num<-2
RCircos.Gene.Name.Plot(data2,name.col,track.num,side)



### hzyg
data4.col<-5 
track.num<-4
side<-"in"
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$hist.color <- "blue"
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.Histogram.Plot(data4,data4.col,track.num,side)

### entropy
data5.col<-5 
track.num<-5
side<-"in"
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$hist.color <- "green"
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.Histogram.Plot(data5,data5.col,track.num,side)

### kmer frac
data6.col<-5 
track.num<-6
side<-"in"
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$hist.color <- "red"
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.Histogram.Plot(data6,data6.col,track.num,side)

dev.off()


