#The .fastq files with the following SRA accession numbers were downloaded:
#SRR609422
#SRR609423
#SRR609424
#SRR609425
#SRR609426
#SRR609427
#SRR609428
#SRR609429
#SRR609430

#The files were analysed using FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) without any extra options.

#The following trimmings were performed using trim_galore! (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/):

#First trimming: Illumina standard

#Second trimming:	
trim_galore --fastqc -o output_directory -a2 TGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -a TGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT --paired file_1.fastq.gz file_2.fastq.gz

#Third trimming:
trim_galore --fastqc -o output_directory -a2 TGGTATCAACGCAGAGTACATGGGATGGCACATGCAGCGCAAGTAGGTCT -a TGGTATCAACGCAGAGTACATGGGATGGCACATGCAGCGCAAGTAGGTCT --paired file_1.fastq.gz file_2.fastq.gz

#Fourth trimming(for SRR609422,SRR609423,SRR609425 only):
trim_galore --fastqc -o output_directory -a2 TGGTATCAACGCAGAGTCCGGTAATCGCATAAAACTTAAAACTTTACAGT -a TGGTATCAACGCAGAGTCCGGTAATCGCATAAAACTTAAAACTTTACAGT --paired file_1.fastq.gz file_2.fastq.gz

#Fifth trimming(for all):	
trim_galore --fastqc -o output_directory -a2 TGGTATCAACGCAGAGTACATGGGCAACTTCTTAGAGGGACAAGTGGCGT -a TGGTATCAACGCAGAGTACATGGGCAACTTCTTAGAGGGACAAGTGGCGT --trim1 --clip_R1 8 --clip_R2 8 --three_prime_clip_R1 3 --three_prime_clip_R2 3 --paired file_1.fastq.gz file_2.fastq.gz
