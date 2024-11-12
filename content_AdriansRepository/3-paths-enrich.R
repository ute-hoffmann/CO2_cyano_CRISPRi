library(tidyverse)
library(ggpubr)
library(patchwork)

library(clusterProfiler)
library(enrichplot)


dat <- read_tsv('2-crispri.tsv')
paths <- read_tsv('2-pathways.tsv')


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

ggsave('3-overview.jpeg', width = 16, height = 9)

################################################################################

sets <-
  paths |>
  select(pathway, gene)

res <-
  c('CO2_30percent', 'CO2_4percent') %>%
  set_names() |>
  map(function(i) {
    dat |>
      select(gene, i = CO2_30percent) |>
      arrange(desc(i)) |>
      with(set_names(i, gene)) |>
      GSEA(TERM2GENE = sets, pvalueCutoff = .01)
  })

res.tbl <-
  res |>
  map(as_tibble) |>
  map(select, - c(ID, leading_edge, core_enrichment)) %>%
  map2(names(.), ~ mutate(.x, condition = .y)) |>
  bind_rows() |>
  select(condition, everything())

write_tsv(res.tbl, '3-gsea.tsv')

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


ggsave('3-gsea.jpeg', width = 20, height = 12)
       