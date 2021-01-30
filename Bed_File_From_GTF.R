library(tidyverse)

#Example of how to make a bed file from a gtf for downstream use e.g. metagene and readcoverage heatmaps in deepTools

#Fetch your primary assembly GTF file
gtf <- read.table("/Your/GTF/file", sep = "\t", stringsAsFactors = F)
gtf <- filter(gtf, grepl("gene", gtf$V3))

#Fetch your PROMPT or genomic features list.  This example assumes your list is a data.frame with a single column of names
prompts <- read.table("/Your/genomic/feature/list",
                      sep = "\t", stringsAsFactors = F)

#The gene names Marchf6, Marchf9, and Septin9 do not match other databases and cause problems with set operations in other analyses.  
#In general, this is true for the septins and march gene families and it is better to change the names accordingly. 
#My list of PROMPTs includes these three problem names

prompts[1504, ] = "March6"
prompts[1576, ] = "March9"
prompts[1954, ] = "Sept9"

#Make a vector of the PROMPT names for downstream filtering
p.vector <- as.vector(prompts$V1)
p.vector <- paste0(p.vector, "\\b")

#Filter the gtf based on the vector of PROMPT names.  The vector used for this type of filtering must have <2500 entries.
p.gtf <- filter(gtf, grepl(paste(p.vector, collapse = "|"), gtf$V9))

#Make a new column and extract the gene name from column 9
p.gtf$symbol <- str_extract(p.gtf$V9, "\\bgene_name.+?\\;")
p.gtf$symbol <- str_remove(p.gtf$symbol, "\\bgene_name\\s")
p.gtf$symbol <- str_remove(p.gtf$symbol, "\\;")

#Grep based filtering base on gene names will also return pseudogenes.  This step is to find and filter out these added genes
z <- as.vector(setdiff(p.gtf$symbol, prompts$V1))
z <- paste0(z, "\\b")
p.gtf2 <- filter(p.gtf, !grepl(paste(z, collapse = "|"), p.gtf$symbol))

#This is to check that all of the pseudogenes have been removed.  It should return an empty vector.
y <- setdiff(p.gtf2$symbol, prompts$V1)

#Now subset the GTF into bed format
p.gtf2 <- p.gtf2[, c(1, 4, 5, 10)]
p.gtf2$V6 <- 0
p.gtf2$V7 <- "."

write.table(p.gtf2, "/Your/new/bed/file", sep = "\t", col.names = F, row.names = F, quote = F)
