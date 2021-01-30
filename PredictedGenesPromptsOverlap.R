library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)

#############################################################
# Overlapping PROMPTs and Predicted Genes, RIKEN cDNAs, etc.

#open prompts gtf created earlier
p.gtf <- read.delim(file = "/vM24.all.prompts.gtf", header = T, stringsAsFactors = F)


#open file of shared prompts from deseq2, edgeR, and limma-voom in gtf format
NIH.prompts <- read.delim(file = "/SharedPrompts.gtf", header = T, stringsAsFactors = F)


#open primary assembly gtf
vM24 <- read.delim("/gencode.vM24.primary_assembly.annotation.gtf", header = F, stringsAsFactors = F, skip = 5)

#open shared predicted genes deseq2 file
pred.genes <- read.delim("/SharedPredictedGenes.txt", header = T)

genes <- filter(vM24, grepl("gene", vM24$V3))
gene_vector <- as.vector(pred.genes$Row.names)
gene_vector <- paste0(gene_vector, "\\b")

pred.gene.coord <- filter(genes, grepl(paste(gene_vector, collapse = "|"), genes$V9))

gene.ranges <- makeGRangesFromDataFrame(pred.gene.coord, keep.extra.columns = T, start.field = "V4", end.field = "V5",
                                        seqnames.field = "V1", strand.field = "V7")                     
prompt.ranges <- makeGRangesFromDataFrame(NIH.prompts, keep.extra.columns = T, start.field = "V4", end.field = "V5",
                                          seqnames.field = "V1", strand.field = "V7")

overlaps <- findOverlaps(gene.ranges, prompt.ranges, type = "any", select = "all", ignore.strand = T)
genes.overlap <- overlaps@from
prompts.overlap <- overlaps@to

genes.subset <- pred.gene.coord[genes.overlap, ]
prompts.subset <- NIH.prompts[prompts.overlap, ]

#########Check for overlap accuracy and number of sense/antisense relationships############

genes.subset$symbol <- str_extract(genes.subset$V9, "\\bgene_name.+?;")
genes.subset$symbol <- str_remove(genes.subset$symbol, "gene_name\\s")
genes.subset$symbol <- str_remove(genes.subset$symbol, ";")
genes.subset <- genes.subset[, c(1, 4, 5, 7, 10)]
names(genes.subset) <- c("g.chr", "g.start", "g.end", "g.strand", "g.symbol")

#gene symbols from either set work; can use either column 12 or 13
prompts.subset <- prompts.subset[, c(1, 4, 5, 7, 12)]
names(prompts.subset) <- c("p.chr", "p.start", "p.end", "p.strand", "p.symbol")

genes.subset <- genes.subset[order(genes.subset$g.start), ]
prompts.subset <- prompts.subset[order(prompts.subset$p.start), ]
combo <- cbind(genes.subset, prompts.subset)

for (i in seq_along(combo$g.chr)) {
  if(combo$g.chr[i] == combo$p.chr[i]) 
    {combo$chr.match[i] = "Match"}
  else(combo$match[i] = "No Match")
  
  if(combo$g.strand[i] == combo$p.strand[i]) 
    {combo$strand.match[i] = "TRUE"}
  else(combo$strand.match[i] = "False")
}

table(unlist(combo[, c(11, 12)]))

combo$g.length <- combo$g.end - combo$g.start


#####################################
#Overlapping PROMPTs and enhancers

eRNA.gff <- read.table(file = "/enhancers.GRCm38.regulatory.build.gff", sep = "\t", header = F)
eRNAs <- read.table(file = "shared_enhancers.txt", sep = "\t", header = F)
eRNAs <- as.vector(eRNAs$V1)
eRNAs <- paste0(eRNAs, "\\b")
enhancers <- filter(eRNA.gff, grepl(paste(eRNAs, collapse = "|"), eRNA.gff$V9))

eRNA.ranges <- makeGRangesFromDataFrame(enhancers, keep.extra.columns = T, start.field = "V4", end.field = "V5",
                                        seqnames.field = "V1")

EPoverlaps <- findOverlaps(prompt.ranges, eRNA.ranges, type = "any", select = "all", ignore.strand = T)
prompt.eRNA.from <- EPoverlaps@from
prompt.eRNA.to <- EPoverlaps@to
p.subset <- NIH.prompts[prompt.eRNA.from, ]
eRNA.subset <- enhancers[prompt.eRNA.to, ]

#####################################
#Overlapping PROMPTs and promoters

promoters.gff <- read.table(file = "/promoters.GRCm38.regulatory.build.gff", sep = "\t", header = F)
NIH.promoters <- read.table(file = "shared_promoters.txt", sep = "\t", header = F)
NIH.promoters <- as.vector(NIH.promoters$V1)
NIH.promoters <- paste0(NIH.promoters, "\\b")
promoters <- filter(promoters.gff, grepl(paste(NIH.promoters, collapse = "|"), promoters.gff$V9))

promoter.ranges <- makeGRangesFromDataFrame(promoters, keep.extra.columns = T, start.field = "V4", end.field = "V5",
                                        seqnames.field = "V1")
PromptPromoter <- findOverlaps(prompt.ranges, promoter.ranges, type = "any", select = "all", ignore.strand = T)
prompt.promoter.from <- PromptPromoter@from
prompt.promoter.to <- PromptPromoter@to
pp.subset <- NIH.prompts[prompt.promoter.from, ]

#####################################
#Overlapping PROMPTs and promoter flanking regions
pfr.gff <- read.table(file = "/pfr.GRCm38.regulatory.build.gff", sep = "\t", header = F)
NIH.pfr <- read.table(file = "shared_pfr.txt", sep = "\t", header = F)
NIH.pfr <- as.vector(NIH.pfr$V1)
NIH.pfr <- paste0(NIH.pfr, "\\b")
pfr <- filter(pfr.gff, grepl(paste(NIH.pfr, collapse = "|"), pfr.gff$V9))

pfr.ranges <- makeGRangesFromDataFrame(pfr, keep.extra.columns = T, start.field = "V4", end.field = "V5",
                                            seqnames.field = "V1")
PromptPFR <- findOverlaps(prompt.ranges, pfr.ranges, type = "any", select = "all", ignore.strand = T)
prompt.pfr.from <- PromptPFR@from
prompt.pfr.to <- PromptPFR@to
pfr.subset <- NIH.prompts[prompt.pfr.from, ]

#####################################
#Overlapping PROMPTs and transcription factor binding sites
tfbs.gff <- read.table(file = "/tfbs.GRCm38.regulatory.build.gff", sep = "\t", header = F)
NIH.tfbs <- read.table(file = "shared_tfbs.txt", sep = "\t", header = F)
NIH.tfbs <- as.vector(NIH.tfbs$V1)
NIH.tfbs <- paste0(NIH.tfbs, "\\b")
tfbs <- filter(tfbs.gff, grepl(paste(NIH.tfbs, collapse = "|"), tfbs.gff$V9))

tfbs.ranges <- makeGRangesFromDataFrame(tfbs, keep.extra.columns = T, start.field = "V4", end.field = "V5",
                                       seqnames.field = "V1")
Prompttfbs <- findOverlaps(prompt.ranges, tfbs.ranges, type = "any", select = "all", ignore.strand = T)
prompt.tfbs.from <- Prompttfbs@from
prompt.tfbs.to <- Prompttfbs@to
tfbs.subset <- NIH.prompts[prompt.tfbs.from, ]

#####################################
#Overlapping PROMPTs and open chromatin regions
ocr.gff <- read.table(file = "/ocr.GRCm38.regulatory.build.gff", sep = "\t", header = F)
NIH.ocr <- read.table(file = "shared_ocr.txt", sep = "\t", header = F)
NIH.ocr <- as.vector(NIH.ocr$V1)
NIH.ocr <- paste0(NIH.ocr, "\\b")
ocr <- filter(ocr.gff, grepl(paste(NIH.ocr, collapse = "|"), ocr.gff$V9))

ocr.ranges <- makeGRangesFromDataFrame(ocr, keep.extra.columns = T, start.field = "V4", end.field = "V5",
                                        seqnames.field = "V1")
Promptocr <- findOverlaps(prompt.ranges, ocr.ranges, type = "any", select = "all", ignore.strand = T)
prompt.ocr.from <- Promptocr@from
prompt.ocr.to <- Promptocr@to
ocr.subset <- NIH.prompts[prompt.ocr.from, ]

#####################################
#Overlapping PROMPTs and CTCF binding sites
CTCF.gff <- read.table(file = "/CTCF.GRCm38.regulatory.build.gff", sep = "\t", header = F)
NIH.CTCF <- read.table(file = "shared_ctcf.txt", sep = "\t", header = F)
NIH.CTCF <- as.vector(NIH.CTCF$V1)
NIH.CTCF <- paste0(NIH.CTCF, "\\b")
CTCF <- filter(CTCF.gff, grepl(paste(NIH.CTCF, collapse = "|"), CTCF.gff$V9))

CTCF.ranges <- makeGRangesFromDataFrame(CTCF, keep.extra.columns = T, start.field = "V4", end.field = "V5",
                                       seqnames.field = "V1")
PromptCTCF <- findOverlaps(prompt.ranges, CTCF.ranges, type = "any", select = "all", ignore.strand = T)
prompt.CTCF.from <- PromptCTCF@from
prompt.CTCF.to <- PromptCTCF@to
ctcf.subset <- NIH.prompts[prompt.CTCF.from, ]

