#---------------------- convert sra to fastq


#for SEQ in $(ls -d VN/* | grep [0-9] |cut -d "/" -f2)
for SEQ in `cat sra_list`; do
	#get sra record from NCBI, write to folder of sra code
	toolkit/bin/prefetch ${SEQ} 
	
	toolkit/bin/fastq-dump --outdir ${SEQ}/ --split-files ${SEQ}/${SEQ}.sra


 	#paste - - - - < *_1.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > 1.fa

done



#--------------------- merge fastq
bbmerge.sh in1=$qc_file1 in2=$qc_file2 out=combined.fastq
