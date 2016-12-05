#For the differential expression analysis the workflow described at the Bioconductor website was followed. The comments can also be found at the same workflow (some of them are copied below to facilitate reading). It can be found at this link http://www.bioconductor.org/help/workflows/rnaseqGene/


#Create the table that contains info about the samples and load it.
sample.table <- read.csv("sample_table.csv", sep = ";", row.names = 1)

#Using the Run column in the sample table, we construct the full paths to the files we want to perform the counting operation on

filenames <- file.path ("G:", "satoh_sorted_bam_files",paste0( "sorted_", sample.table$Run, ".bam"))

#We indicate in Bioconductor that these files are BAM files using the BamFileList function from the Rsamtools package:
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)

#Change the names of the bamfiles because they should not be all the same:
names(bamfiles) <- sample.table$Run

#Next, we need to read in the gene model that will be used for counting reads/fragments. We will read the gene model from UCSC GTF file, using makeTxDbFromGFF from the GenomicFeatures package:
library("GenomicFeatures")

#We indicate that none of our sequences (chromosomes) are circular using a 0-length character vector:
genesfile <- file.path ( "C:", "R_wd", "required_files", "known_genes", "hg38_genes", "Homo_sapiens.GRCh38.84.gtf")
txdb <- makeTxDbFromGFF (genesfile, format="gtf", circ_seqs=character())

#The following line produces a GRangesList of all the exons grouped by gene:
ebg <- exonsBy(txdb, by="gene")

#Change the chromosomes names to match the gtf file:
first, get the names of the chromosomes:

library("GenomeInfoDb")
chrUCSC <- extractSeqlevels(species="Homo_sapiens", style="UCSC")
chrENSEMBL <- extractSeqlevels(species="Homo_sapiens", style="Ensembl")

#Create a vector with the "from-to" names:

newnames <- chrUCSC #To
names(newnames) <- chrENSEMBL #From

#Rename:
ebg2 <- renameSeqlevels(ebg,newnames)

#The function summarizeOverlaps from the GenomicAlignments package will do the counting. This produces a SummarizedExperiment object that contains a variety of information about the experiment:

library("GenomicAlignments")

se <- summarizeOverlaps(features=ebg2, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        fragments=FALSE,
                        ignore.strand=TRUE )


#The colData slot, so far empty, should contain all the metadata.
#Because we used a column of sampleTable to produce the bamfiles vector, we know the columns of se are in the same order as the rows of sampleTable. We can assign the sampleTable as the colData of the summarized experiment, by converting it into a DataFrame and using the assignment function:
colData(se) <- DataFrame(sample.table)

#Note: it is prefered in R that the first level of a factor be the reference level (e.g. control, or untreated samples), so we can relevel the dex factor:

se$Phenotype <- relevel(se$Phenotype, "normal")

#Once we have our fully annotated SummarizedExperiment object, we can construct a DESeqDataSet object from it that will then form the starting point of the analysis. We add an appropriate design for the analysis:

library("DESeq2")
dds <- DESeqDataSet(se, design = ~ Phenotype)

#We remove rows of the DESeqDataSet that have no counts, or only a single count across all samples:

dds <- dds[ rowSums(counts(dds)) > 1, ]

#We perform the rlog transformation:

rld <- rlog(dds)

#We can run the differential expression pipeline on the raw counts with a single call to the function DESeq:

dds2 <- DESeq(dds)

#Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula:

res <- results(dds2)

#We keep results with a p-value lower than 0.05:

res.sig <- res[!is.na(res$padj),]
res.sig <- res.sig[res.sig$padj < 0.05, ]

#Use Ensembl ids to add gene symbols
library("AnnotationDbi")  #Base package
library("org.Hs.eg.db")  #Human

res.sig$SYMBOL <- mapIds(org.Hs.eg.db,
                          keys=row.names(res.sig),
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")


res$SYMBOL <- mapIds(org.Hs.eg.db,
                                keys=row.names(res),
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first")

#Export the results to be used at the pathway analysis:

dir.create("outputs")
library("xlsx")
write.xlsx(res.sig, "outputs/sig.results.GE.xlsx")
write.table(res, "outputs/results.group1")
