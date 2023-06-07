log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(fgsea)
library(msigdbr)
library(tidyverse)

msigdbr_df <- msigdbr(species = snakemake@params[["species"]], category = "H")
pathways.hallmark <- split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

res <- read_tsv(snakemake@input[["table"]])

res2 <- res %>%
  dplyr::select(gene, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene) %>%
  summarize(stat = mean(stat))

ranks <- deframe(res2)
fgseaRes <- fgsea(pathways = pathways.hallmark, stats = ranks)

#svg('gsea_plot.svg')
#ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
#  geom_col(aes(fill=padj<0.05)) +
#  coord_flip() +
#  labs(x="Pathway", y="Normalized Enrichment Score",
#       title="Hallmark pathways NES from GSEA") +
#  theme_minimal()
#dev.off()

pdf(snakemake@output[["plot"]])
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.hallmark[topPathways], ranks, fgseaRes, gseaParam=0.5)
dev.off()

save.image()