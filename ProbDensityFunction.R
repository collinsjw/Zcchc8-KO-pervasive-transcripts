library(tidyverse)
library(data.table)
library(scales)

#Prep gtf files for creating rug plots by adding a column with the start site of all genes and prompts.
mus.gtf <- fread("/Your/GTF/of/expressed/genes")
sims.prompts <-fread("/Your/GTF/of/PROMPTs")

gtf <- mus.gtf
#gtf <- gtf[grepl("chr", gtf$V1) & !grepl("chrM", gtf$V1) & !grepl("chrY", gtf$V1) & grepl("gene", gtf$V3), ]
for (i in seq_along(gtf$V1)) {
  if(gtf$V7[i] == "+")
    {gtf$start[i] = gtf$V4[i]}
  else if(gtf$V7[i] == "-")
  {gtf$start[i] = gtf$V5[i]}
}

sims.prompts <- sims.prompts[grepl("chr", sims.prompts$V1), ]
for (i in seq_along(sims.prompts$V1)) {
  if(sims.prompts$V7[i] == "+")
  {sims.prompts$start[i] = sims.prompts$V4[i]}
  else if(sims.prompts$V7[i] == "-")
  {sims.prompts$start[i] = sims.prompts$V5[i]}
}

#########################################################
#Genes and PROMPTs per Mb per chromosome

j <- table(unlist(gtf[, 1]))
j <- as.data.frame(j, stringsAsFactors = F)
j <- j[gtools::mixedorder(j$Var1), ]

k <- table(unlist(sims.prompts[, 1]))
k <- as.data.frame(k, stringsAsFactors = F)
k <- k[gtools::mixedorder(k$Var1), ]

l <- merge(j, k, by = "Var1")
l <- rename(l, chr = "Var1", genes = "Freq.x", prompts = "Freq.y")
l <- l[gtools::mixedorder(l$chr), ]
#Mouse chromosome lengths as of January 2021.  Follow the RefSeq hyperlinks https://www.ncbi.nlm.nih.gov/genome?term=mus%20musculus
#for the exact chromosome size.
l$chrLength <- c(195154279, 181755017, 159745316, 156860686, 151758149, 149588044, 144995196, 130127694, 124359700, 130530862,
          121973369, 120092757, 120883175, 125139656, 104073951, 98008968, 95294699, 90720763, 61420004, 169476592)

l <- l %>% pivot_longer(cols = c(genes, prompts), names_to = "type", values_to = "number")

l$PerMb <- (l$number/l$chrLength)*10^6 

##For scientific notation in axis ticks
fancy_scientific <- function(x) {
  # turn in to character string in scientific notation
  x <- format(x, scientific = TRUE)
  # keep zero looking nice
  x <- gsub("0e\\+00","0", x)
  # quote the part before the exponent to keep all the digits
  x <- gsub("^(.*)e", "'\\1'e", x)
  # turn the 'e+/-' into math exponential format
  x <- gsub("e-", "%*%10^-", x)
  x <- gsub("e\\+", "%*%10^", x)
  # return this as an expression
  parse(text=x)
}

##When plotting the barchart, this will put the chromosomes in correct order 
chr.order <- j$Var1
l$chr <- factor(l$chr, levels = chr.order)

##Create vector of chromosome names for facet labelling
chr.vect <- c("Chromosome 1", "Chromosome 2", "Chromosome 3", "Chromosome 4", "Chromosome 5", "Chromosome 6", "Chromosome 7",
              "Chromosome 8", "Chromosome 9", "Chromosome 10", "Chromosome 11", "Chromosome 12", "Chromosome 13", "Chromosome 14",
              "Chromosome 15", "Chromosome 16", "Chromosome 17", "Chromosome 18", "Chromosome 19", "Chromosome X")

m <- ggplot(data = l) +
  geom_bar(aes(x = chr, y = PerMb, fill = type), stat = "identity", position = position_dodge(-0.4))

m +
  scale_fill_manual(values = c("grey60", "firebrick"), labels = c("Genes", "PROMPTs")) + 
  labs(fill = "") +
  xlab("") + ylab("Frequency (per Mb)") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = c(0.2, 0.9),
    legend.key.size = unit(0.25, "in"),
    legend.direction = "vertical",
    legend.title.align = 0,
    legend.background = element_rect(fill="transparent"),
    legend.key = element_blank(),
    axis.text.x=element_text(angle=45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "black", linetype = "dashed", size = 0.1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    # Change axis line
    axis.line = element_line(color = "black")
  )

#Kernel Density Estimation for Probability Density Function
sims.prompts$V1 <- factor(sims.prompts$V1, levels = chr.order, labels = chr.vect)
gtf$V1 <- factor(gtf$V1, levels = chr.order, labels = chr.vect)

#Sometimes this throws a "polygon edge not found" error.  There are many reasons for this error so 
#I suggest searching stackoverflow.com for solutions.

tiff("/Users/Joshua/JHU_dissertation/Figures/OnlyExpressedGenesNrd0.tiff", units = "in", width = 8.5, height = 11, res = 600)
f <- ggplot() + 
  labs(x = "Chromosome Position (bp)", y = "Probability Density Function") +
  geom_density(data = sims.prompts, 
               aes(x = start), color = "red", bw = "nrd0", n = 200, kernel = "g") + 
  geom_density(data = gtf, 
                aes(x = start), color = "black", bw = "nrd0", n = 200, kernel = "g") +
  geom_rug(data = sims.prompts, aes(x = start), color = "red", sides = "b", 
           alpha = 0.5, length = unit(0.1, "npc")) +
  geom_rug(data = gtf, aes(x = start), color = "black", sides = "t", 
           alpha = 0.05, length = unit(0.1, "npc")) +
  facet_wrap(~V1, scales = "free_x", ncol = 4) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)), labels = fancy_scientific) +
  scale_x_continuous(labels = fancy_scientific) +
  theme(
    text=element_text(family = "Arial"),
    axis.text.x = element_text(angle = 20, vjust = 1.2, hjust = 0.8, size = 6),
    axis.text.y = element_text(angle = 45, vjust = 1.8, hjust = 0.4, size = 6),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x=unit(0.1, "cm"),
    panel.border = element_rect(color = "black", fill=NA, size=0.5)
  )
f
dev.off()

#Checking for strand bias in PROMPTs
f2 <- ggplot() + 
  geom_density(data=sims.prompts[sims.prompts$V7 == "+", ], 
               aes(x = start), color = "red", bw = "SJ") + 
  geom_density(data=sims.prompts[sims.prompts$V7 == "-", ], 
               aes(x = start), color = "black", bw = "SJ") +
  geom_rug(data = sims.prompts[sims.prompts$V7 == "+", ], 
           aes(x = start), color = "red", sides = "b", 
           alpha = 0.5, length = unit(0.1, "npc")) +
  geom_rug(data = sims.prompts[sims.prompts$V7 == "-", ], 
           aes(x = start), color = "black", sides = "b", 
           alpha = 0.5, length = unit(0.1, "npc")) +
  facet_wrap(~V1, scales = "free_x", ncol = 4)
f2



