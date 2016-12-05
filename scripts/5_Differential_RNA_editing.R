#Place the final files from the editing meta analysis in the folder called "input_files"

#Load the files
adar.ad <- read.table("input_files/ADAR.AD", header = T)
adar.normal <- read.table("input_files/ADAR.normal", header = T)
apobec.ad <- read.table("input_files/APOBEC.AD", header = T)
apobec.normal <- read.table("input_files/APOBEC.normal", header = T)

#Function to remove the first column and the snp144 column and filter by editing frequency
#threshold is the minimum frequency to keep
#freq.col is the column containing the frequencies
#rm.col is a numerical vector containing the columns to remove
filter.frequency <- function(table, min.freq = 0.00, max.freq = 0.9, freq.col = 16, rm.col = c(1,13)){
        table <- table[which(table[,freq.col] >= min.freq),]
        table <- table[which(table[,freq.col] <= max.freq),]
        table <- table[,-rm.col]
        row.names(table) <- seq_len(nrow(table))
        table
}

adar.ad <- filter.frequency(adar.ad)
adar.normal <- filter.frequency(adar.normal)
apobec.ad <- filter.frequency(apobec.ad)
apobec.normal <- filter.frequency(apobec.normal)

#Function to generate a table that shows the distribution of the rna editing
#positions
editing.dist <- function(table){
        table[,"Func.refGene"] <- as.character(table[,"Func.refGene"])
        types <- unique(table[, "Func.refGene"])
        new.table <- as.data.frame(matrix(nrow = 1, ncol = length(types)))
        names(new.table) <- types
        for (i in seq_along(types)){
                rows <- nrow(table[which(table[,"Func.refGene"] == types[i]), ])
                new.table[1, types[i]] <- rows
        }
        new.table
        
}

#Manually make the barplot of the distribution using the function above

#Function to:
#1. remove unwanted columns
#2. rename "Start" and "End" to "Position"
reorganize <- function(table){
        new.table <- table[,-c(3,12,13,15)]
        names(new.table)[2] <- "Position"
        new.table
}

adar.ad <- reorganize(adar.ad)
adar.normal <- reorganize(adar.normal)
apobec.ad <- reorganize(apobec.ad)
apobec.normal <- reorganize(apobec.normal)

#Function to select common and unique sites
#the prepositions "common" "normal" and "ad" will automaticaly be added
# The function adds a final column with a logical value that tests whether
#the two frequences differ for a certain percentage.
select.sites <- function (normal.table, ad.table, normal.name, ad.name, common.name, percentage){
        #make an index with the common positions
        index1 <- logical(length = nrow(normal.table))
        for (i in seq_len(nrow(normal.table))){
                index1[i] <- normal.table[i, "Position"] %in% ad.table[,"Position"]
        }
        
        #make the uniques to normal.table
        uniq.norm <- normal.table[!index1,]
        row.names(uniq.norm) <- seq_len(nrow(uniq.norm))
        assign(paste0("unique.", normal.name), uniq.norm, .GlobalEnv)
        
        #Do the same for ad.table
        index2 <- logical(length = nrow(ad.table))
        for (i in seq_len(nrow(ad.table))){
                index2[i] <- ad.table[i, "Position"] %in% normal.table[,"Position"]
        }
        
        #make the uniques to normal.table
        uniq.ad <- ad.table[!index2,]
        row.names(uniq.ad) <- seq_len(nrow(uniq.ad))
        assign(paste0("unique.", ad.name), uniq.ad, .GlobalEnv)
        
        
        #Make the table with the common positions
        common <- normal.table[index1,]
        names(common)[ncol(common)] <- normal.name
        common[,ad.name] <- NA
        
        for (i in seq_len(nrow(common))){
                common[i,ad.name] <- ad.table[which(ad.table[,"Position"] == 
                                                            common[i, "Position"]), ncol(ad.table)]
        }
        
        #Insert a final column with a logic value that tests whether
        #the two frequencies differ for a certain percentage
        perc <- paste0(as.character(percentage*100), "%")
        common[,paste0(perc, " difference")] <- NA
        for (i in seq_len(nrow(common))){
                common[i,ncol(common)] <- abs(common[i,ncol(common) - 1] -
                        common[i,ncol(common) - 2]) >= percentage
        }
        
        row.names(common) <- seq_len(nrow(common))
        assign(paste0("common.", common.name), common, .GlobalEnv)
}


#Select UTR3
select.sites(adar.normal[which(adar.normal[,"Func.refGene"] == "UTR3"),],
             adar.ad[which(adar.ad[,"Func.refGene"] == "UTR3"),],
             "adar.normal.UTR3",
             "adar.ad.UTR3",
             "adar.UTR3",
             0.2)

#select exonic
select.sites(adar.normal[which(adar.normal[,"Func.refGene"] == "exonic"),],
             adar.ad[which(adar.ad[,"Func.refGene"] == "exonic"),],
             "adar.normal.exonic",
             "adar.ad.exonic",
             "adar.exonic",
             0.2)

#APOBEC Select UTR3
select.sites(apobec.normal[which(apobec.normal[,"Func.refGene"] == "UTR3"),],
             apobec.ad[which(apobec.ad[,"Func.refGene"] == "UTR3"),],
             "apobec.normal.UTR3",
             "apobec.ad.UTR3",
             "apobec.UTR3",
             0.2)

#APOBEC select exonic
select.sites(apobec.normal[which(apobec.normal[,"Func.refGene"] == "exonic"),],
             apobec.ad[which(apobec.ad[,"Func.refGene"] == "exonic"),],
             "apobec.normal.exonic",
             "apobec.ad.exonic",
             "apobec.exonic",
             0.2)
 
#Copy the results from the gene expression analysis and load them
results.GE <- read.table("input_files/gene expression/results.GE", header = T)

#You need to split double genes so that the add.change function can work
#Function to split double genes in the ref.Gene field in seperate rows
split.genes <- function(table) {
        #make an index with the rows that have double genes
        index <- grep(",", table$Gene.refGene, fixed = T)
        #if there are no double genes, then return the table
        if (length(index) == 0) {table}
        #else, build two tables and then combine them
        else {
                table2 <- table[-c(index),]
                table3 <- table[index,]
                #build a table with twice as many rows to put the double genes
                double <- table3[rep(seq_len(nrow(table3)), each=2),]
                double$Gene.refGene <- as.character(double$Gene.refGene)
                for (i in seq_len(nrow(table3))) {
                        data.row1 <- i*2 -1
                        data.row2 <- i*2
                        gene1 <- strsplit(as.character(table3$Gene.refGene[i]),
                                          ",", fixed = T)[[1]][1]
                        gene2 <- strsplit(as.character(table3$Gene.refGene[i]),
                                          ",", fixed = T)[[1]][2]
                        
                        double[data.row1,"Gene.refGene"] <- gene1
                        double[data.row2,"Gene.refGene"] <- gene2
                }
                #combine the two tables
                final <- rbind(table2,double)
                row.names(final) <- seq_len(nrow(final))
                final
        }
}

common.adar.exonic <- split.genes(common.adar.exonic)
common.adar.UTR3 <- split.genes(common.adar.UTR3)
unique.adar.ad.exonic <- split.genes(unique.adar.ad.exonic)
unique.adar.ad.UTR3 <- split.genes(unique.adar.ad.UTR3)
unique.adar.normal.exonic <- split.genes(unique.adar.normal.exonic)
unique.adar.normal.UTR3 <- split.genes(unique.adar.normal.UTR3)

common.apobec.exonic <- split.genes(common.apobec.exonic)
common.apobec.UTR3 <- split.genes(common.apobec.UTR3)
unique.apobec.ad.exonic <- split.genes(unique.apobec.ad.exonic)
unique.apobec.ad.UTR3 <- split.genes(unique.apobec.ad.UTR3)
unique.apobec.normal.exonic <- split.genes(unique.apobec.normal.exonic)
unique.apobec.normal.UTR3 <- split.genes(unique.apobec.normal.UTR3)


#Function to add a column to the tables that will contain the log2fold change.
#symbol.col is the column of the results table that contains the gene symbols.
#gene.col is the column of the table that contains the gene symbols
add.change <- function (table, results, gene.col = "Gene.refGene") {
        library("AnnotationDbi")  #Base package
        library("org.Hs.eg.db")  #Human
        table[,"log2FoldChange"] <- NA
        table[,"significant change"] <- NA
        table[,gene.col] <- as.character(table[,gene.col])
        for (i in seq_len(nrow(table))) {                
                gene <- function(x){
                        tryCatch(mapIds(org.Hs.eg.db,
                                                keys=x,
                                                column="ENSEMBL",
                                                keytype="SYMBOL",
                                                multiVals="first"), error = function(e) {NA})
                }
                
                if (gene(table[i, gene.col]) %in% row.names(results)) {
                        table[i, "log2FoldChange"] <-
                                results[which(row.names(results) == gene(table[i, gene.col])), "log2FoldChange"]
                        
                        if (!is.na(results[which(row.names(results) == gene(table[i, gene.col])), "padj"])) {
                                if (results[which(row.names(results) == gene(table[i, gene.col])), "padj"] <= 0.05) {
                                        table[i, "significant change"] <- TRUE
                                }
                                else {
                                        table[i, "significant change"] <- FALSE
                                }
                        }
                        else {
                                table[i, "significant change"] <- FALSE
                        }
                }
                else {
                        print("NAs returned")
                }
        }
        table
}


common.adar.exonic <- add.change(common.adar.exonic, results.GE)
common.adar.UTR3 <- add.change(common.adar.UTR3, results.GE)
unique.adar.ad.exonic <- add.change(unique.adar.ad.exonic, results.GE)
unique.adar.ad.UTR3 <- add.change(unique.adar.ad.UTR3, results.GE)
unique.adar.normal.exonic <- add.change(unique.adar.normal.exonic, results.GE)
unique.adar.normal.UTR3 <- add.change(unique.adar.normal.UTR3, results.GE)

common.apobec.exonic <- add.change(common.apobec.exonic, results.GE)
common.apobec.UTR3 <- add.change(common.apobec.UTR3, results.GE)
unique.apobec.ad.exonic <- add.change(unique.apobec.ad.exonic, results.GE)
unique.apobec.ad.UTR3 <- add.change(unique.apobec.ad.UTR3, results.GE)
unique.apobec.normal.exonic <- add.change(unique.apobec.normal.exonic, results.GE)
unique.apobec.normal.UTR3 <- add.change(unique.apobec.normal.UTR3, results.GE)

#Save the tables containing only UTR3 positions to be used in the pathway analysis
save(common.apobec.UTR3, common.adar.UTR3, unique.apobec.normal.UTR3,
     unique.adar.ad.UTR3, unique.apobec.ad.UTR3, unique.adar.normal.UTR3,
     file = "tables_for_pathways")
