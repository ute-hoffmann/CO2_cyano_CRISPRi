library(tidyverse)
library(ggpubr)
library(patchwork)
library(openxlsx)

library(clusterProfiler)
library(enrichplot)


dat <- read_tsv('data/2-crispri.tsv')
paths <- read_tsv('data/2-pathways.tsv')


################################################################################

p1 <-
  dat |>
  select(gene, contains('CO2')) |>
  pivot_longer(- gene) |>
  mutate_at('name', str_remove, 'CO2_') |>
  mutate_at('name', str_replace, 'percent', '%') |>
  ggdensity('value', fill = 'name') +
  xlab('weighted mean fitness') +
  ggsci::scale_fill_jama(name = 'CO2') +
  theme_pubr(18)

p2 <-
  dat |>
  select(gene, contains('CO2')) |>
  pivot_longer(- gene) |>
  mutate_at('name', str_remove, 'CO2_') |>
  mutate_at('name', str_replace, 'percent', '%') |>
  ggecdf('value', color = 'name') +
  xlab('weighted mean fitness') +
  ggsci::scale_color_jama(name = 'CO2') +
  theme_pubr(18)

p3 <-
  dat |>
  ggplot(aes(CO2_30percent, CO2_4percent, color = type)) +
  geom_point() +
  ggsci::scale_color_d3() +
  xlab('30% CO2') +
  ylab('4% CO2') +
  theme_pubr(18)
  
  
((p1 / p2) | p3) +
  plot_annotation(tag_levels = 'A')

ggsave('results/3-overview.svg', width = 16, height = 9)

################################################################################

sets <-
  paths |>
  select(pathway, gene)

res <- list()
set.seed(1234)
res$CO2_30percent <- dat |> select(gene, i = CO2_30percent) |> arrange(desc(i)) |> with(set_names(i, gene)) |> GSEA(TERM2GENE = sets, pvalueCutoff = .01, seed = T)
res$CO2_4percent <- dat |> select(gene, i = CO2_4percent) |> arrange(desc(i)) |> with(set_names(i, gene)) |> GSEA(TERM2GENE = sets, pvalueCutoff = .01, seed = T)

res.tbl <-
  res |>
  map(as_tibble) |>
  map(select, - c(ID, leading_edge, core_enrichment)) %>%
  map2(names(.), ~ mutate(.x, condition = .y)) |>
  bind_rows() |>
  select(condition, everything())

write_tsv(res.tbl, 'results/3-gsea.tsv')

res |>
  map(function(i) {
    # i <- res$CO2_30percent
    p <- gseaplot2(i, 1:nrow(i), 
                   base_size = 14,
                   color = RColorBrewer::brewer.pal(nrow(i), 'Paired'))
    
    # Change axis text for clarity
    p[[3]] <-  p[[3]] + ylab('weighted mean fitness')
    wrap_elements(full = print(p))
  }) -> ps

ps %>%
  map2(names(.), ~ .x + ggtitle(.y |>
                                  str_remove('CO2_') |>
                                  str_replace('percent', '%'))) |>
  reduce(.f = `|`) &
  theme_pubr(18)


ggsave('results/3-gsea.svg', width = 20, height = 12)
       

#return ranked list of genes in significant pathways for 4% CO2
get.ranked.list <- function(i) {
  pw.genes <- paths %>% filter(pathway %in% i) %>% select(gene)
  all.genes <- dat %>% select(gene, CO2_4percent)
  ranked.list <- inner_join(pw.genes, all.genes, by = "gene") %>% arrange(desc(CO2_4percent))
  ranked.list <- inner_join(ranked.list, dat, by = c("gene","CO2_4percent"))
}
l.ranked.lists <- lapply(res$CO2_4percent@result$ID, get.ranked.list)
write.xlsx(l.ranked.lists, "results/3-ranked.lists-CO2_4percent.xlsx", sheetName = substr(res$CO2_4percent@result$ID, 1, 28))

#return ranked list of genes in significant pathways for 30% CO2
get.ranked.list <- function(i) {
  pw.genes <- paths %>% filter(pathway %in% i) %>% select(gene)
  all.genes <- dat %>% select(gene, CO2_30percent)
  ranked.list <- inner_join(pw.genes, all.genes, by = "gene") %>% arrange(desc(CO2_30percent))
  ranked.list <- inner_join(ranked.list, dat, by = c("gene","CO2_30percent"))
}
l.ranked.lists <- lapply(res$CO2_30percent@result$ID, get.ranked.list)
write.xlsx(as.list(l.ranked.lists), "results/3-ranked.lists-CO2_30percent.xlsx", sheetName = substr(res$CO2_30percent@result$ID, 1, 28))
