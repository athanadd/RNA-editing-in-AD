#The strand correction takes place using a custom R function that detects the position of each gene from a reference file and adds this information to the position of the editing event. The reference file is obtained from UCSC and is Homo_sapiens.GRCh38.83.csv for human datasets.

#To use the script, copy and paste the whole folder to a new directory.
#Place the output, annotated by ANNOVAR files in a subdirectory called "annotated",
#place the REDItools converted files in a subdirectory called "converted".


#Function to add the strand information. The output files will be placed
#in a subdirectory called "with_strand". File names will be preserved.
add.strand <- function(species) {
        library(dplyr)
        
        #Get the full paths of the input files
        annotated = list.files("./annotated", full.names = T)
        converted = list.files("./converted", full.names = T)
        
        #Set the variables of the annotation files according to the species
        if (species == "mouse") {
                RG = read.csv("./required_files/mouse/Mus_musculus_ensembl_dataset.csv")
                RG = RG[!duplicated(RG$gene_name), ]
        }
        else if (species == "human") {
                RG = read.csv("./required_files/human/Homo_sapiens.GRCh38.83.csv")
                RG = RG[!duplicated(RG$gene_name),]
        }
        
        #Create the ouputs folder
        dir.create("./with_strand")
        
        #Start the loop for each file
        for (y in seq_along(annotated)) {
                TA = read.csv(annotated[y])
                ind1 = nrow(TA)
                strand_vector = c()
                
                for (i in seq(ind1)) {
                        # if the position is located in intergenic, upstream or downstream region,
                        # then it is assigned directly with ".", because strand cannot be defined
                        
                        if (TA[i, "Func.refGene"] %in% c("intergenic",
                                                         "upstream",
                                                         "downstream")) {
                                print("no_strand")
                                strand_vector = c(strand_vector, ".")
                        }
                        else {
                                pos = unlist(strsplit(as.character(TA[i, "Gene.refGene"]), split = ","))
                                
                                # one gene name for this position
                                if (length(pos) == 1) {
                                        if (length(which(RG$gene_name == pos)) == 1) {
                                                ind2 = which(RG$gene_name == pos)
                                                strand_vector = c(
                                                        strand_vector,
                                                        as.character(RG$strand[ind2])
                                                )
                                        }
                                        else {
                                                # "missing" if the gene name is not matched
                                                strand_vector = c(strand_vector,
                                                                  "missing")
                                        }
                                }
                                else {
                                        # two gene names in this position
                                        if (length(which(RG$gene_name == pos[1])) ==
                                            1) {
                                                ind2 = which(RG$gene_name == pos[1])
                                                strand_vector = c(
                                                        strand_vector,
                                                        as.character(RG$strand[ind2])
                                                )
                                        }
                                        else if (length(which(RG$gene_name == pos[2])) ==
                                                 1) {
                                                ind2 = which(RG$gene_name == pos[2])
                                                strand_vector = c(
                                                        strand_vector,
                                                        as.character(RG$strand[ind2])
                                                )
                                        }
                                        else {
                                                strand_vector = c(strand_vector,
                                                                  "missing")
                                        }
                                }
                        }
                }
                
                # the strand information has been gathered in the strand_vector
                new_TA_1 = TA[, which(colnames(TA) == "Chr"):which(colnames(TA) ==
                                                                           "Alt")]
                new_TA_1 = mutate(new_TA_1, strand = strand_vector)
                new_TA_2 = cbind(new_TA_1, TA[, which(colnames(TA) == "Func.refGene"):ncol(TA)])
                output = list.files("./converted")[y]
                output = gsub("converted_", "", output)
                output = paste0("./with_strand/", output, ".csv")
                write.csv(new_TA_2, file = output)
        }
}

#Choose between "human" or "mouse" for species
#and run this function!
add.strand("human")

#Function to filter the RNA editing positions and keep edititng events
#that appear on the correct strand. The function also splits the tables
#according to the enzyme (ADAR, APOBEC).
#The output files will be placed in a subdirectory called "final_files".
### WARNING - The function does NOT remove SNPs, they should already be
#absent by providing them as input to REDItools.
filter.files <- function(){
        
        library(xlsx)
        library(dplyr)
        
        #Get the full paths of the input files
        converted = list.files("./converted", full.names = T)
        with.strand = list.files("./with_strand/", full.names = T)
        
        #Create the ouputs folders
        dir.create("./final_files/ADAR_all_positions/tables", recursive = T)
        dir.create("./final_files/APOBEC_all_positions/tables", recursive = T)
        dir.create("./final_files/ADAR_all_positions/xlsx", recursive = T)
        dir.create("./final_files/APOBEC_all_positions/xlsx", recursive = T)
        
        #Start the loop for each file
        for (y in seq_along(converted)) {
                
                output = list.files("./converted")[y]
                output = gsub("converted_", "", output)
                table_annovar = read.csv(with.strand[y])
                index=read.table(converted[y])
                data1 = table_annovar%>%
                        mutate(Read_depth = index[,7])%>%
                        mutate(Quality_score = index[,8])%>%
                        mutate(Variant_allele_frequency = index[,14])%>%
                        mutate(P_value = index[,15])
                
                #For APOBEC
                
                data2 = data1 %>%
                        (function(x){
                                index1 = c(which(data1$Ref == "C" & data1$Alt == "T" & data1$strand %in% c(".","+","missing")), 
                                           which(data1$Ref == "G" & data1$Alt == "A" & data1$strand %in% c(".","-","missing")))
                                x[index1,]
                        }) %>%
                        (function(x){x})
                
                data2 <- data2[which(data2$snp144 == "."), ]
                
                write.table(row.names = FALSE,data2,file=paste0("./final_files/APOBEC_all_positions/tables/APOBEC.",output))
                write.xlsx(row.names = FALSE,data2,file=paste0("./final_files/APOBEC_all_positions/xlsx/APOBEC.",output,".xlsx"))
                
                #For ADAR
                
                data2 = data1%>%
                        (function(x){
                                index1 = c(which(data1$Ref == "A" & data1$Alt == "G" & data1$strand %in% c(".","+","missing")), 
                                           which(data1$Ref == "T" & data1$Alt == "C" & data1$strand %in% c(".","-","missing")))
                                x[index1,]
                        }) %>%
                        (function(x){x})
                data2 <- data2[which(data2$snp144 == "."), ]
                
                write.table(row.names = FALSE,data2,file=paste0("./final_files/ADAR_all_positions/tables/ADAR.",output))
                write.xlsx(row.names = FALSE,data2,file=paste0("./final_files/ADAR_all_positions/xlsx/ADAR.",output,".xlsx"))
        }
}

#Just run the function
filter.files()
