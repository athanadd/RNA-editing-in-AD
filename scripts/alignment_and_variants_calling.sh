#The sequences of the same phenotype were merged in one file using the cat command:
cat SRR609422-26.fastq.gz > AD_phenotype.fastq.gz
cat SRR609427-30.fastq.gz > Control_phenotype.fastq.gz

#The alignment was then performed using tophat2 (https://ccb.jhu.edu/software/tophat/index.shtml):

tophat genes.gtf -g 1 --no-novel-indels --no-coverage-search --mate-inner-dist=-46 --mate-std-dev=40 --transcriptome-index=hg38/transcriptome_data --read-mismatches 3 --read-edit-dist 3 --read-realign-edit-dist 0 -a 12 -m 1 -o output_directory hg38/Bowtie2Index/genome input_file_1 input_file_2

#The bam files were sorted and indexed using samtools (http://samtools.sourceforge.net/):

samtools sort -o sorted_file.bam accepted_hits.bam
samtools index sorted_file.bam

#A Blat Correction analysis was performed using the original script provided from the REDItools website (http://reditools.sourceforge.net):

python REDItoolBlatCorrection.py -i input_bam_file -o output_directory -V -f reference_fasta_file -F reference_2bit_file

#The main analysis was performed by REDItool Denovo, a script of the REDItools suite:

python REDItoolDenovo.py -i input_bam_file -o output_directory snps_file -p -f reference_fasta_file -e -l -u -B blat_correction_directory -W splice_sites

#The snps and splice_sites files were prepared according to the official REDItools documentation.

#The output table with the significant sites was converted in order to be processed by ANNOVAR using the following R script:

R --vanilla --quiet <<RSCRIPT
.libPaths("path_to_library")
library(dplyr)


  reditools = read.table(skip = 1, "table_sig_from_REDItools")
  Alt_base = sapply(as.character(reditools[,11]), function(x){
    unlist(strsplit(x, split = ""))[[2]]
    })
  reditools%>%
    (function(x){
      cbind(x[,1], x[,2], x[,2], x[,3], Alt_base, x[,4:ncol(x)])
      })%>%
    (function(y){
      write.table(y, file = "converted_file" , 
                  row.names = FALSE, col.names = FALSE,quote = FALSE)
    })
RSCRIPT

#The converted file is used as an input to ANNOVAR (http://annovar.openbioinformatics.org/en/latest/):

ANNOVAR/table_annovar.pl converted_file annovar_files/humandb/hg38 -buildver hg38 -remove -protocol refGene,snp144 -operation g,f -nastring . -csvout -out annotated_file
