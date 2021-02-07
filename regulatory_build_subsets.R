library(tidyverse)

#Open the regulatory build gff file
gff <- read.delim("/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff", header = F)

#This gff file lacks the customary "chr" prefix for the chromosomes so we must add one.
gff <- transform(gff, V1 = sprintf('chr%s', V1))

#Filter the gff file based on the feature
eRNA.gff <- filter(gff, grepl("Enhancer", gff$V9))
eRNA.gff$V3 <- "enhancer"
write.table(eRNA.gff, file = "/enhancers.GRCm38.regulatory.build.gff", row.names = F, col.names = F, quote = F, sep = "\t")

promoter.gff <- filter(gff, grepl("Promoter", gff$V9) & !grepl("Flanking", gff$V9))
promoter.gff$V3 <- "promoter"
write.table(promoter.gff, file = "/promoters.GRCm38.regulatory.build.gff", row.names = F, col.names = F, quote = F, sep = "\t")

pfr.gff <- filter(gff, grepl("Promoter Flanking Region", gff$V9))
pfr.gff$V3 <- "pfr"
write.table(pfr.gff, file = "/pfr.GRCm38.regulatory.build.gff", row.names = F, col.names = F, quote = F, sep = "\t")

ctcf.gff <- filter(gff, grepl("CTCF", gff$V9))
ctcf.gff$V3 <- "ctcf"
write.table(ctcf.gff, file = "/ctcf.GRCm38.regulatory.build.gff", row.names = F, col.names = F, quote = F, sep = "\t")

tfbs.gff <- filter(gff, grepl("Transcription", gff$V9))
tfbs.gff$V3 <- "tfbs"
write.table(tfbs.gff, file = "/tfbs.GRCm38.regulatory.build.gff", row.names = F, col.names = F, quote = F, sep = "\t")

ocr.gff <- filter(gff, grepl("Open chromatin", gff$V9))
ocr.gff$V3 <- "ocr"
write.table(ocr.gff, file = "/ocr.GRCm38.regulatory.build.gff", row.names = F, col.names = F, quote = F, sep = "\t")


