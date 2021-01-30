library(tidyverse)
library(venn)
library(ggpolypath)

setwd("/Your/working/directory/")

ngd <- read.delim("GSE165689/genes/ZcKOvsWT_genes_deseq2.txt")
nge <- read.delim("GSE165689/genes/ZcKOvsWT_genes_edgeR.txt")
ngl <- read.delim("GSE165689/genes/ZcKOvsWT_genes_limmaVoom.txt")

hgd <- read.delim("GSE126108/genes/ZcKOvsWT_genes_deseq2.txt")
hge <- read.delim("GSE126108/genes/ZcKOvsWT_genes_edgeR.txt")
hgl <- read.delim("GSE126108/genes/ZcKOvsWT_genes_limmaVoom.txt")

sgd <- read.delim("GSE127790/genes/ZcKOvsWT_genes_deseq2.txt")
sge <- read.delim("GSE127790/genes/ZcKOvsWT_genes_edgeR.txt")
sgl <- read.delim("GSE127790/genes/ZcKOvsWT_genes_limmaVoom.txt")

npd <- read.delim("GSE165689/prompts/ZcKOvsWT_prompts_deseq2.txt")
npe <- read.delim("GSE165689/prompts/ZcKOvsWT_prompts_edgeR.txt")
npl <- read.delim("GSE165689/prompts/ZcKOvsWT_prompts_limmaVoom.txt")

hpd <- read.delim("GSE126108/prompts/ZcKOvsWT_prompts_deseq2.txt")
hpe <- read.delim("GSE126108/prompts/ZcKOvsWT_prompts_edgeR.txt")
hpl <- read.delim("GSE126108/prompts/ZcKOvsWT_prompts_limmaVoom.txt")

spd <- read.delim("GSE127790/prompts/ZcKOvsWT_prompts_deseq2.txt")
spe <- read.delim("GSE127790/prompts/ZcKOvsWT_prompts_edgeR.txt")
spl <- read.delim("GSE127790/prompts/ZcKOvsWT_prompts_limmaVoom.txt")

#################################################
# Filter sets based on log2FC > 0.58 or < -0.58 #
#################################################
#GSE165689 Zcchc8 KO SIMS Cells

npd_pos <- filter(npd, npd$log2FoldChange > 0.58 & npd$padj <= 0.01)
npd_neg <- filter(npd, npd$log2FoldChange < -0.58 & npd$padj <= 0.01)

npe_pos <- filter(npe, npe$logFC > 0.58 & npe$FDR <= 0.01)
npe_neg <- filter(npe, npe$logFC < -0.58 & npe$FDR <= 0.01)

npl_pos <- filter(npl, npl$logFC > 0.58 & npl$adj.P.Val <= 0.01)
npl_neg <- filter(npl, npl$logFC < -0.58 & npl$adj.P.Val <= 0.01)

ngd_pos <- filter(ngd, ngd$log2FoldChange > 0.58 & ngd$padj <= 0.01)
ngd_neg <- filter(ngd, ngd$log2FoldChange < -0.58 & ngd$padj <= 0.01)

nge_pos <- filter(nge, nge$logFC > 0.58 & nge$FDR <= 0.01)
nge_neg <- filter(nge, nge$logFC < -0.58 & nge$FDR <= 0.01)

ngl_pos <- filter(ngl, ngl$logFC > 0.58 & ngl$adj.P.Val <= 0.01)
ngl_neg <- filter(ngl, ngl$logFC < -0.58 & ngl$adj.P.Val <= 0.01)

#GSE126108 KO mouse
hpd_pos <- filter(hpd, hpd$log2FoldChange > 0.58 & hpd$padj <= 0.01)
hpd_neg <- filter(hpd, hpd$log2FoldChange < -0.58 & hpd$padj <= 0.01)

hpe_pos <- filter(hpe, hpe$logFC > 0.58 & hpe$FDR <= 0.01)
hpe_neg <- filter(hpe, hpe$logFC < -0.58 & hpe$FDR <= 0.01)

hpl_pos <- filter(hpl, hpl$logFC > 0.58 & hpl$adj.P.Val <= 0.01)
hpl_neg <- filter(hpl, hpl$logFC < -0.58 & hpl$adj.P.Val <= 0.01)

hgd_pos <- filter(hgd, hgd$log2FoldChange > 0.58 & hgd$padj <= 0.01)
hgd_neg <- filter(hgd, hgd$log2FoldChange < -0.58 & hgd$padj <= 0.01)

hge_pos <- filter(hge, hge$logFC > 0.58 & hge$FDR <= 0.01)
hge_neg <- filter(hge, hge$logFC < -0.58 & hge$FDR <= 0.01)

hgl_pos <- filter(hgl, hgl$logFC > 0.58 & hgl$adj.P.Val <= 0.01)
hgl_neg <- filter(hgl, hgl$logFC < -0.58 & hgl$adj.P.Val <= 0.01)

#GSE127790 KO ES cells
spd_pos <- filter(spd, spd$log2FoldChange > 0.58 & spd$padj <= 0.01)
spd_neg <- filter(spd, spd$log2FoldChange < -0.58 & spd$padj <= 0.01)

spe_pos <- filter(spe, spe$logFC > 0.58 & spe$FDR <= 0.01)
spe_neg <- filter(spe, spe$logFC < -0.58 & spe$FDR <= 0.01)

spl_pos <- filter(spl, spl$logFC > 0.58 & spl$adj.P.Val <= 0.01)
spl_neg <- filter(spl, spl$logFC < -0.58 & spl$adj.P.Val <= 0.01)

sgd_pos <- filter(sgd, sgd$log2FoldChange > 0.58 & sgd$padj <= 0.01)
sgd_neg <- filter(sgd, sgd$log2FoldChange < -0.58 & sgd$padj <= 0.01)

sge_pos <- filter(sge, sge$logFC > 0.58 & sge$FDR <= 0.01)
sge_neg <- filter(sge, sge$logFC < -0.58 & sge$FDR <= 0.01)

sgl_pos <- filter(sgl, sgl$logFC > 0.58 & sgl$adj.P.Val <= 0.01)
sgl_neg <- filter(sgl, sgl$logFC < -0.58 & sgl$adj.P.Val <= 0.01)

################################
# GSE165689 SIMS Cells Venn Diagrams #
################################
A <- rbind(npd_pos, npd_neg)
B <- rbind(npe_pos, npe_neg)
C <- rbind(npl_pos, npl_neg)
A1 <- A$external_gene_name
B1 <- B$external_gene_name
C1 <- C$external_gene_name

a <- npd_pos$external_gene_name
b <- npe_pos$external_gene_name
c <- npl_pos$external_gene_name

d <- npd_neg$external_gene_name
e <- npe_neg$external_gene_name
f <- npl_neg$external_gene_name

G <- rbind(ngd_pos, ngd_neg)
H <- rbind(nge_pos, nge_neg)
I <- rbind(ngl_pos, ngl_neg)
G1 <- G$external_gene_name
H1 <- H$external_gene_name
I1 <- I$external_gene_name

g <- ngd_pos$external_gene_name
h <- nge_pos$external_gene_name
i <- ngl_pos$external_gene_name

j <- ngd_neg$external_gene_name
k <- nge_neg$external_gene_name
l <- ngl_neg$external_gene_name

#GSE165689 PROMPTs
nih.prompts <- list("LimmaVoom" = C1, "DESeq2" = A1, "edgeR" = B1)
nih.prompts.venn <- venn(nih.prompts, zcolor = "style", ilcs = 1, sncs = 1.2)
#GSE165689 Positive Prompts
nih.pos.prompts <- list("DESeq2" = a, "edgeR" = b, "Limma-Voom" = c)
npp.venn <- venn(nih.pos.prompts, zcolor = "style", ilcs = 1, sncs = 1.2)

npp.int <- attributes(npp.venn)
npp.int <- npp.int[["intersections"]][[7]]
z<- as.data.frame(npp.int)

#GSE165689 Negative Prompts
nih.neg.prompts <- list("DESeq2" = d, "edgeR" = e, "limmaVoom" = f)
nnp.venn <- venn(nih.neg.prompts, zcolor = "style")

nnp.int <- attributes(nnp.venn)
nnp.int <- nnp.int[["intersections"]][[5]]
y <- as.data.frame(nnp.int)


#GSE165689 Genes
nih.genes <- list("LimmaVoom" = I1, "DESeq2" = G1, "edgeR" = H1)
nih.genes.venn <- venn(nih.genes, zcolor = "style", ilcs = 1, sncs = 1.2)

#GSE165689 Positive Genes
nih.pos.genes <- list("DESeq2" = g, "edgeR" = h, "limmaVoom" = i)
npg.venn <- venn(nih.pos.genes, zcolor = "style")

npg.int <- attributes(npg.venn)
npg.int <- npg.int[["intersections"]][[6]]
x <- as.data.frame(npg.int)

#GSE165689 Negative Genes
nih.neg.genes <- list("DESeq2" = j, "edgeR" = k, "limmaVoom" = l)
nng.venn <- venn(nih.neg.genes, zcolor = "style")

nng.int <- attributes(nng.venn)
nng.int <- nng.int[["intersections"]][[5]]
y <- as.data.frame(nng.int)

#################################
# GSE126108 KO MouseVenn Diagrams #
#################################

hop.a <- hpd_pos$external_gene_name
hop.b <- hpe_pos$external_gene_name
hop.c <- hpl_pos$external_gene_name

hop.d <- hpd_neg$external_gene_name
hop.e <- hpe_neg$external_gene_name
hop.f <- hpl_neg$external_gene_name

hop.g <- hgd_pos$external_gene_name
hop.h <- hge_pos$external_gene_name
hop.i <- hgl_pos$external_gene_name

hop.j <- hgd_neg$external_gene_name
hop.k <- hge_neg$external_gene_name
hop.l <- hgl_neg$external_gene_name

#GSE126108 Positive Prompts
hop.pos.prompts <- list("DESeq2" = hop.a, "edgeR" = hop.b, "limmaVoom" = hop.c)
hpp.venn <- venn(hop.pos.prompts, zcolor = "style")

hpp.int <- attributes(hpp.venn)
hpp.int <- hpp.int[["intersections"]][[7]]

#GSE126108 Negative Prompts
hop.neg.prompts <- list("DESeq2" = hop.d, "edgeR" = hop.e, "limmaVoom" = hop.f)
hnp.venn <- venn(hop.neg.prompts, zcolor = "style")

hnp.int <- attributes(hnp.venn)
hnp.int <- hnp.int[["intersections"]][[7]]

#GSE126108 Positive Genes
hop.pos.genes <- list("DESeq2" = hop.g, "edgeR" = hop.h, "limmaVoom" = hop.i)
hpg.venn <- venn(hop.pos.genes, zcolor = "style")

hpg.int <- attributes(hpg.venn)
hpg.int <- hpg.int[["intersections"]][[6]]


#GSE126108 Negative Genes
hop.neg.genes <- list("DESeq2" = hop.j, "edgeR" = hop.k, "limmaVoom" = hop.l)
hng.venn <- venn(hop.neg.genes, zcolor = "style")

hng.int <- attributes(hng.venn)
hng.int <- hng.int[["intersections"]][[7]]

#####################################
#GSE127790 KO ES Cells Venn Diagrams #
#####################################

sh.a <- spd_pos$external_gene_name
sh.b <- spe_pos$external_gene_name
sh.c <- spl_pos$external_gene_name

sh.d <- spd_neg$external_gene_name
sh.e <- spe_neg$external_gene_name
sh.f <- spl_neg$external_gene_name

sh.g <- sgd_pos$external_gene_name
sh.h <- sge_pos$external_gene_name
sh.i <- sgl_pos$external_gene_name

sh.j <- sgd_neg$external_gene_name
sh.k <- sge_neg$external_gene_name
sh.l <- sgl_neg$external_gene_name

#GSE127790 Positive Prompts
sh.pos.prompts <- list("DESeq2" = sh.a, "edgeR" = sh.b, "limmaVoom" = sh.c)
spp.venn <- venn(sh.pos.prompts, zcolor = "style")

spp.int <- attributes(spp.venn)
spp.int <- spp.int[["intersections"]][[6]]

#GSE127790 Negative Prompts
sh.neg.prompts <- list("DESeq2" = sh.d, "edgeR" = sh.e, "limmaVoom" = sh.f)
snp.venn <- venn(sh.neg.prompts, zcolor = "style")

snp.int <- attributes(snp.venn)
snp.int <- snp.int[["intersections"]][[5]]

#GSE127790 Positive Genes
sh.pos.genes <- list("DESeq2" = sh.g, "edgeR" = sh.h, "limmaVoom" = sh.i)
spg.venn <- venn(sh.pos.genes, zcolor = "style")

spg.int <- attributes(spg.venn)
spg.int <- spg.int[["intersections"]][[7]]


#GSE127790 Negative Genes
sh.neg.genes <- list("DESeq2" = sh.j, "edgeR" = sh.k, "limmaVoom" = sh.l)
sng.venn <- venn(sh.neg.genes, zcolor = "style")

sng.int <- attributes(sng.venn)
sng.int <- sng.int[["intersections"]][[5]]



########################################################
# Correlation between Feature Expression and Cell Type #
########################################################

#Venn diagram for GSE165689 and GSE126108 Prompts
nh.prompts <- list("GSE165689 PROMPTs" = npp.int, "GSE126108 PROMPTs" = hpp.int)
nh.prompts.venn <- venn(nh.prompts, zcolor = "style")
nh.prompts.int <- attributes(nh.prompts.venn)
nh.prompts.int <- nh.prompts.int[["intersections"]][[3]]

#Venn diagram for Negative GSE165689 and GSE126108 Prompts
nh.neg.prompts <- list("GSE165689 PROMPTs" = nnp.int, "GSE126108 PROMPTs" = hnp.int)
nh.neg.prompts.venn <- venn(nh.neg.prompts, zcolor = "style")


#Venn diagram for GSE165689 and GSE127790 Prompts
ns.prompts <- list("GSE165689 PROMPTs" = npp.int, "GSE127790 PROMPTs" = spp.int)
ns.prompts.venn <- venn(ns.prompts, zcolor = "style")
ns.prompts.int <- attributes(ns.prompts.venn)
ns.prompts.int <- ns.prompts.int[["intersections"]][[7]]

#Venn diagram for Negative GSE165689 and GSE127790 Prompts
ns.neg.prompts <- list("GSE165689 PROMPTs" = nnp.int, "GSE127790 PROMPTs" = snp.int)
ns.neg.prompts.venn <- venn(ns.neg.prompts, zcolor = "style")

#Venn diagram for GSE165689, GSE126108, and GSE127790 Prompts
nhs.prompts <- list("GSE165689 PROMPTs" = npp.int, "GSE127790 PROMPTs" = spp.int, "GSE126108 PROMPTs" = hpp.int)
nhs.prompts.venn <- venn(nhs.prompts, zcolor = "style", ilcs = 1)
nhs.prompts.int <- attributes(nhs.prompts.venn)
nhs.prompts.int <- nhs.prompts.int[["intersections"]][[7]]
x <- as.data.frame(nhs.prompts.int)

#Venn diagram for Negative GSE165689, GSE126108, and GSE127790 Prompts
nhs.neg.prompts <- list("GSE165689 PROMPTs" = nnp.int, "GSE127790 PROMPTs" = snp.int, "GSE126108 PROMPTs" = hnp.int)
nhs.neg.prompts.venn <- venn(nhs.neg.prompts, zcolor = "style", ilcs = 1)

#Venn diagram for combined GSE165689, GSE126108, GSE127790 postive genes
pos.gn <- list("NPG" = npg.int, "HPG" = hpg.int, "SPG" = spg.int)
pos.int.venn <- venn(pos.gn, zcolor = "style")
pos.gn.int <- attributes(pos.int.venn)
pos.gn.int <- pos.gn.int[["intersections"]][[7]]
###################Filter out predicted genes#############################
df.a <- as.data.frame(pos.gn.int)
df.a <- filter(df.a, !grepl("Rik", df.a$pos.gn.int))
df.a <- filter(df.a, !grepl("Gm", df.a$pos.gn.int))

#Venn diagram for GSE165689 and GSE126108 postive genes
nh.pos.gn <- list("NPG" = npg.int, "HPG" = hpg.int)
nh.pos.int.venn <- venn(nh.pos.gn, zcolor = "style")
nh.pos.gn.int <- attributes(nh.pos.int.venn)
nh.pos.gn.int <- nh.pos.gn.int[["intersections"]][[3]]
df1 <- as.data.frame(nh.pos.gn.int)
df2 <- filter(df1, !grepl("Rik", df1$nh.pos.gn.int))
df2 <- filter(df2, !grepl("Gm", df2$nh.pos.gn.int))

#Venn diagram for combined GSE165689, GSE126108, GSE127790 negative genes
neg.gn <- list("NPG" = nng.int, "HPG" = hng.int, "SPG" = sng.int)
neg.int.venn <- venn(neg.gn, zcolor = "style")

#Venn diagram for GSE165689 and GSE126108 negative genes
nh.neg.gn <- list("NPG" = nng.int, "HPG" = hng.int)
nh.neg.int.venn <- venn(nh.neg.gn, zcolor = "style")
nh.neg.gn.int <- attributes(nh.neg.int.venn)
nh.neg.gn.int <- nh.neg.gn.int[["intersections"]][[3]]

#Venn diagram for GSE165689 and GSE127790 positive genes
ns.pos.gn <- list("NPG" = npg.int, "SPG" = spg.int)
ns.pos.int.venn <- venn(ns.pos.gn, zcolor = "style")
ns.pos.gn.int <- attributes(ns.pos.int.venn)
ns.pos.gn.int <- ns.pos.gn.int[["intersections"]][[3]]

#Venn diagram for GSE165689 and GSE127790 negative genes
ns.neg.gn <- list("NNG" = nng.int, "SNG" = sng.int)
ns.neg.int.venn <- venn(ns.neg.gn, zcolor = "style")
ns.neg.gn.int <- attributes(ns.neg.int.venn)
ns.neg.gn.int <- ns.neg.gn.int[["intersections"]][[3]]


##################################################################
###################### Btbd7 Expression ##########################
##################################################################

btbd7.ngd <- read.delim("GSE165689/genes/BtKOvsWT_genes_deseq2.txt")
btbd7.nge <- read.delim("GSE165689/genes/BtKOvsWT_genes_edgeR.txt")
btbd7.ngl <- read.delim("GSE165689/genes/BtKOvWT_genes_limmaVoom.txt")

btbd7.npd <- read.delim("GSE165689/prompts/BtKOvsWT_prompts_deseq2.txt")
btbd7.npe <- read.delim("GSE165689/prompts/BtKOvsWT_prompts_edgeR.txt")
btbd7.npl <- read.delim("GSE165689/prompts/BtKOvWT_prompts_limmaVoom.txt")

btbd7.npd_pos <- filter(btbd7.npd, btbd7.npd$log2FoldChange > 0.58 & btbd7.npd$padj <= 0.01)
btbd7.npd_neg <- filter(btbd7.npd, btbd7.npd$log2FoldChange < -0.58 & btbd7.npd$padj <= 0.01)

btbd7.npe_pos <- filter(btbd7.npe, btbd7.npe$logFC > 0.58 & btbd7.npe$FDR <= 0.01)
btbd7.npe_neg <- filter(btbd7.npe, btbd7.npe$logFC < -0.58 & btbd7.npe$FDR <= 0.01)

btbd7.npl_pos <- filter(btbd7.npl, btbd7.npl$logFC > 0.58 & btbd7.npl$adj.P.Val <= 0.01)
btbd7.npl_neg <- filter(btbd7.npl, btbd7.npl$logFC < -0.58 & btbd7.npl$adj.P.Val <= 0.01)

btbd7.ngd_pos <- filter(btbd7.ngd, btbd7.ngd$log2FoldChange > 0.58 & btbd7.ngd$padj <= 0.01)
btbd7.ngd_neg <- filter(btbd7.ngd, btbd7.ngd$log2FoldChange < -0.58 & btbd7.ngd$padj <= 0.01)

btbd7.nge_pos <- filter(btbd7.nge, btbd7.nge$logFC > 0.58 & btbd7.nge$FDR <= 0.01)
btbd7.nge_neg <- filter(btbd7.nge, btbd7.nge$logFC < -0.58 & btbd7.nge$FDR <= 0.01)

btbd7.ngl_pos <- filter(btbd7.ngl, btbd7.ngl$logFC > 0.58 & btbd7.ngl$adj.P.Val <= 0.01)
btbd7.ngl_neg <- filter(btbd7.ngl, btbd7.ngl$logFC < -0.58 & btbd7.ngl$adj.P.Val <= 0.01)

btbd7.a <- btbd7.npd_pos$external_gene_name
btbd7.b <- btbd7.npe_pos$external_gene_name
btbd7.c <- btbd7.npl_pos$external_gene_name

btbd7.d <- btbd7.npd_neg$external_gene_name
btbd7.e <- btbd7.npe_neg$external_gene_name
btbd7.f <- btbd7.npl_neg$external_gene_name

btbd7.g <- btbd7.ngd_pos$external_gene_name
btbd7.h <- btbd7.nge_pos$external_gene_name
btbd7.i <- btbd7.ngl_pos$external_gene_name

btbd7.j <- btbd7.ngd_neg$external_gene_name
btbd7.k <- btbd7.nge_neg$external_gene_name
btbd7.l <- btbd7.ngl_neg$external_gene_name

#GSE165689 Positive Prompts
btbd7.pos.prompts <- list("DESeq2" = btbd7.a, "edgeR" = btbd7.b, "limmaVoom" = btbd7.c)
btbd7.npp.venn <- venn(btbd7.pos.prompts, zcolor = "style")

btbd7.npp.int <- attributes(btbd7.npp.venn)
btbd7.npp.int <- btbd7.npp.int[["intersections"]][[6]]

#GSE165689 Negative Prompts
btbd7.neg.prompts <- list("DESeq2" = btbd7.d, "edgeR" = btbd7.e, "limmaVoom" = btbd7.f)
btbd7.nnp.venn <- venn(btbd7.neg.prompts, zcolor = "style")

btbd7.nnp.int <- attributes(btbd7.nnp.venn)
btbd7.nnp.int <- btbd7.nnp.int[["intersections"]][[5]]


#GSE165689 Positive Genes
btbd7.pos.genes <- list("DESeq2" = btbd7.g, "edgeR" = btbd7.h, "limmaVoom" = btbd7.i)
btbd7.npg.venn <- venn(btbd7.pos.genes, zcolor = "style")

btbd7.npg.int <- attributes(btbd7.npg.venn)
btbd7.npg.int <- btbd7.npg.int[["intersections"]][[6]]
x <- as.data.frame(btbd7.npg.int)

#GSE165689 Negative Genes
btbd7.neg.genes <- list("DESeq2" = btbd7.j, "edgeR" = btbd7.k, "limmaVoom" = btbd7.l)
btbd7.nng.venn <- venn(btbd7.neg.genes, zcolor = "style")

btbd7.nng.int <- attributes(btbd7.nng.venn)
btbd7.nng.int <- btbd7.nng.int[["intersections"]][[7]]
y <- as.data.frame(btbd7.nng.int)


