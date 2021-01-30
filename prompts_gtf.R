library(tidyverse)

#Import primary assembly gtf file
VM24.gtf <- read.table("/gencode.vM24.primary_assembly.annotation.gtf", sep = "\t", header = F, stringsAsFactors = F)

#Make working copy
prompts <- VM24.gtf

#Filter for gene start site
prompts <- filter(VM24.gtf, grepl("gene", VM24.gtf$V3))

#Define PROMPT end coordinates
prompts$V3 = "prompt"
for (i in seq_along(prompts$V1)) {
  if(prompts$V7[i] == "+") prompts$V5[i] = prompts$V4[i]
  if(prompts$V7[i] == "-") prompts$V4[i] = prompts$V5[i]
}

#Define PROMPT start coordinates
for (i in seq_along(prompts$V1)) {
  if (prompts$V7[i] == "+") prompts$V4[i] = prompts$V4[i] - 3000
  if (prompts$V7[i] == '-') prompts$V5[i] = prompts$V5[i] + 3000
}

#Write GTF file
write.table(prompts, file = "/vM24.all.prompts.gtf", quote = F, row.names = FALSE, col.names = FALSE, sep = "\t")


