#To perform the pathway analysis the pathview package was used (http://pathview.r-forge.r-project.org/).


#Copy and load the file with the gene expression results from the expression analysis
sig.results.GE <- read.xlsx("input_files/sig.results.GE.xlsx", sheetIndex = 1)

#Copy and load the tables of the genes (only UTR3) from the diff.editing analysis
load("input_files/tables_for_pathways")

#Merge the tables into one table
#change the names of the frequency columns for the unique
names(unique.adar.ad.UTR3)[11] <- "ad.UTR3"
names(unique.adar.normal.UTR3)[11] <- "normal.UTR3"
names(unique.apobec.ad.UTR3)[11] <- "ad.UTR3"
names(unique.apobec.normal.UTR3)[11] <- "normal.UTR3"
names(common.adar.UTR3)[c(11,12)] <- c("normal.UTR3", "ad.UTR3")
names(common.apobec.UTR3)[c(11,12)] <- c("normal.UTR3", "ad.UTR3")
#add new column for frequencies
unique.adar.ad.UTR3$normal.UTR3 <- NA
unique.apobec.ad.UTR3$normal.UTR3 <- NA
unique.adar.normal.UTR3$ad.UTR3 <- NA
unique.apobec.normal.UTR3$ad.UTR3 <- NA
#reorder the columns
unique.adar.ad.UTR3 <- unique.adar.ad.UTR3[c(1:10,14,11,12,13)]
unique.apobec.ad.UTR3 <- unique.apobec.ad.UTR3[c(1:10,14,11,12,13)]
unique.adar.normal.UTR3 <- unique.adar.normal.UTR3[c(1:11,14,12,13)]
unique.apobec.normal.UTR3 <- unique.apobec.normal.UTR3[c(1:11,14,12,13)]
#remove the 20% change column from the common
common.adar.UTR3 <- common.adar.UTR3[,-13]
common.apobec.UTR3 <- common.apobec.UTR3[,-13]
#Add a new column for the enzyme
unique.adar.ad.group2.UTR3$enzyme <- "adar"
unique.apobec.ad.group2.UTR3$enzyme <- "apobec"
unique.adar.normal.group2.UTR3$enzyme <- "adar"
unique.apobec.normal.group2.UTR3$enzyme <- "apobec"
common.adar.group2.UTR3$enzyme <- "adar"
common.apobec.group2.UTR3$enzyme <- "apobec"
#Bind the table
all <- rbind(common.apobec.group2.UTR3, common.adar.group2.UTR3,
             unique.apobec.normal.group2.UTR3, unique.adar.normal.group2.UTR3,
             unique.apobec.ad.group2.UTR3, unique.adar.ad.group2.UTR3)
#Add a new column in the tables for ad, one for common and another
#for normal that will have a value of:
#0 for the commons
#-1 for the AD
#1 for the Normal
t.common <- nrow(all) - nrow(rbind(unique.apobec.normal.group2.UTR3, unique.adar.normal.group2.UTR3,
                                   unique.apobec.ad.group2.UTR3, unique.adar.ad.group2.UTR3))

t.normal <- nrow(all) - t.common - nrow(rbind(unique.apobec.ad.group2.UTR3, unique.adar.ad.group2.UTR3))

t.ad <- nrow(all) - t.common - t.normal
#commons
all$colcode1 <- c(rep(0, times= t.common), rep(NA, times=t.normal), rep(NA, times = t.ad))
#unique.normal
all$colcode2 <- c(rep(NA, times= t.common), rep(1, times=t.normal), rep(NA, times = t.ad))
#unique.ad
all$colcode3 <- c(rep(NA, times= t.common), rep(NA, times=t.normal), rep(-1, times = t.ad))


#prepare the matrix to be used by pathview
names.all <- as.character(all$Gene.refGene)
mat.all <- as.matrix(all[,c(16,17,18)])
row.names(mat.all) <- names.all  

################Visualization

####LOOP for all the pathways
ids <- c("05014", "00190", "03013",
         "03050", "04064", "04120", "04141",
         "04142", "04210", "04721", "04724",
         "05010")
suffixes <- c("ALS", "Oxidative phosphorylation",
              "RNA transport", "Proteasome", "NFKB",
              "ubiquitin mediated proteolysis",
              "Protein processing in ER", "Lysosome",
              "Apoptosis", "Synaptic vesicle cycle",
              "Glutamatergic synapse", "Alzheimer's Disease")

library("pathview")

for (i in seq_along(ids)){
        pv.out <- pathview(gene.data = mat.all[,1:3], pathway.id = ids[i], species = "hsa",
                           out.suffix = suffixes[i], kegg.native = T, gene.idtype = "SYMBOL",
                           res=600, kegg.dir = "./kegg",
                           low = list(gene = "red", cpd = "blue"),
                           mid = list(gene = "#00FFFF", cpd = "gray"),
                           high = list(gene = "green", cpd = "yellow"))
}
