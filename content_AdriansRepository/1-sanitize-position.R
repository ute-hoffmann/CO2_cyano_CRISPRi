# Lookup all exact matches for guide sequences in the genome
# The positions should match the claimed info (ideally uniquely)

library(tidyverse)
library(ggpubr)

library(patchwork)

library(furrr)

plan(multisession)

###############################################################################

guides <-
  'raw-data/guides.fasta' |>
  Biostrings::readDNAStringSet()

genome <-
  'GCF_000009725.1_ASM972v1_genomic.fna.gz' |>
  Biostrings::readDNAStringSet()

###############################################################################
# Find locations of all aligned sequences following the
# Efficient genome searching tutorial, section
# "Finding an arbitrary nucleotide pattern in an entire genome"
# https://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/GenomeSearching.pdf
# Idea: reverse complement pattern before search

# search but including PAM
g2 <- paste0(guides, 'NGG') |>
  Biostrings::DNAStringSet()
names(g2) <- names(guides)

patterns <-
  bind_rows(
    tibble(
      guide = names(g2),
      seq = g2 |> as.character(),
      strand = '+'
    ),
    tibble(
      guide = names(g2),
      seq = g2 |>
        Biostrings::reverseComplement() |>
        as.character(),
      strand = '-'
    )
  )

search.helper <- function(x) {
  # fixed FALSE for N in NGG pattern
  Biostrings::vmatchPattern(x, genome, fixed = FALSE) |>
    as.data.frame() |>
    transmute(
      seq = x,
      start,
      end,
      chr = names(genome)[group]
    )
}

    

guide.search <-
  patterns |>
  pull(seq) |>
  unique() |>
  future_map(search.helper) |>
  bind_rows()



###############################################################################

guide.pos <-
  patterns |>
  left_join(guide.search, 'seq',
            relationship = "many-to-many") |>
  drop_na()

write_tsv(guide.pos, '1-seq-search.tsv.gz')

###############################################################################
# guide.pos |>
#   count(guide) |>
#   filter(n > 1)
# 0

###############################################################################

p1 <-
  guides |>
  IRanges::width() |>
  tibble(i = _) |>
  count(i) |>
  ggbarplot('i', 'n') +
  # scale_y_log10() +
  xlab('Guide sequence length') +
  ylab('No. guides') +
  theme_pubr(18) +
  ggtitle(
    paste(guides |> length() |> prettyNum(big.mark = ','), 'guides'),
    'input/Synechocystis_v2_trimmed.fasta'
  )

p2 <-
  guide.pos |>
  count(chr) |>
  mutate(len = set_names(
    IRanges::width(genome),
    names(genome)
  )[chr]) |>
  ggscatter(
    'len', 'n',
    add = 'reg.line',
    cor.coef = TRUE
  ) +
  scale_x_log10() +
  scale_y_log10() +
  xlab('Replicon length, nt') +
  ylab('Guide sequence matches') +
  theme_pubr(18)

p1 + p2 +
  plot_annotation(tag_levels = 'A')
ggsave('1-guide-seq-lookup.jpeg', width = 12, height = 6)

###############################################################################
# > w <- IRanges::width(guides)
# > table(w)
# w
# 17   18   19   20   21   22   23 
# 4178 3132 3398 3863 3126 3771    2 
# > guides[w > 22]
# DNAStringSet object of length 2:
#   width seq                                                                                                          names               
# [1]    23 ATGTTCTTCTCAGGGTGTAAACC                                                                                      ctrl3|0
# [2]    23 TATGCTTGGTTTGGACTGACTCA                                                                                      ctrl6|0

###############################################################################
# Sanitize distance between nearest gene and the one targeted

gff <-
  'GCF_000009725.1_ASM972v1_genomic.gff.gz' |>
  rtracklayer::import.gff() |>
  filter(type != 'region')

gr <-
  guide.pos |>
  transmute(
    seqnames = str_remove(chr, ' .*$'),
    start, end, strand,
    guide,
    row = 1:n()
  ) |>
  plyranges::as_granges()

gg.pairs <-
  plyranges::join_nearest(
    gr,
    gff |> filter(type == 'gene'),
    distance = TRUE
  )

###############################################################################

gg <-
  gr |>
  as_tibble() |>
  mutate(guide.gene = str_remove(guide, '\\|.*$')) |>
  left_join(
    'raw-data/names.tsv' |> read_tsv(),
    c('guide.gene' = 'gene')
  ) |>
  select(guide, row, locus) |>
  unique() |>
  inner_join(
    gff |>
      as_tibble() |>
      select(ID, old_locus_tag), 
    c('locus' = 'old_locus_tag')
  ) |>
  drop_na() |>
  mutate(gff.row = set_names(1:length(gff), gff$ID)[ID])

dat.target <-
  gg |>
  mutate(
    dist =  GenomicRanges::distance(gr[row], gff[gff.row], ignore.strand = TRUE),
    same.strand = as.character(gr@strand[row]) == as.character(gff@strand[gff.row])
  )

pa <-
  dat.target |>
  left_join(
    guide.pos |>
      count(guide) |>
      filter(n > 1) |>
      transmute(guide, guide.multi = TRUE),
    'guide'
  ) |>
  mutate_at('guide.multi', replace_na, FALSE) |>
  mutate(
    same.strand = ifelse(
      same.strand,
      'guide sequence on same strand as target',
      'guide sequence anti-sense to target'
    ),
    guide.multi = ifelse(guide.multi, 
                         'guide sequence with multiple occurences',
                         'guide sequence with single occurence')
  ) |>
  ggplot(aes(dist + 1, group = paste(same.strand, guide.multi),
             color = same.strand)) +
             # linetype = guide.multi)) +
  stat_ecdf() +
  scale_x_log10() +
  xlab('distance to target, nt\n(incl. pseudo count)') +
  theme_pubr(18)
  
###############################################################################
# Comparison nearest gene

dat.nearest <-
  gg.pairs |>
  as_tibble() |>
  transmute(
    guide,
    guide.gene = str_remove(guide, '\\|.*$'),
    guide.row = row,
    ID,
    Name,
    locus_tag,
    old_locus_tag,
    distance
  )  |>
  left_join(
    'raw-data/names.tsv' |> read_tsv(),
    c('guide.gene' = 'gene')
  ) |>
  left_join(
    guide.pos |>
      count(guide) |>
      filter(n > 1) |>
      transmute(guide, guide.multi = TRUE),
    'guide'
  ) |>
  mutate_at('guide.multi', replace_na, FALSE) |>
  mutate(
    match = case_when(
      old_locus_tag == locus ~ 'nearest gene as targeted',
      is.na(locus) ~ 'target not in gff (non-coding?)',
      TRUE ~ 'nearest gene different than target'
    ),
    guide.multi = ifelse(guide.multi, 
                         'guide sequence with multiple occurences',
                         'guide sequence with single occurence')
  )

pb <-
  dat.nearest |>
  ggplot(aes(distance, group = paste(match, guide.multi),
             color = match)) +
             # linetype = guide.multi)) +
  stat_ecdf() +
  xlab('distance to nearest genomic neighbor') +
  theme_pubr(18)

dat.nearest |>
  filter(
    match != 'nearest gene as targeted',
    ! str_starts(guide.gene, 'nc_')
  ) |>
  View()

################################################################################
pa + pb +
  plot_annotation(tag_level = 'A') &
  theme(legend.direction = 'vertical')

ggsave('1-guide-dist.jpeg', width = 18, height = 8)

################################################################################
dat.target |>
  filter(dist > 100) |>
  with(bind_cols(
    gr[row] |>
      select(- row) |>
      as_tibble() |>
      transmute(
        guide,
        guide.chr = seqnames,
        guide.start = start,
        guide.end = end
      ),
    gff[gff.row] |>
      as_tibble() |>
      transmute(
        gene = ID,
        old_locus_tag,
        Name,
        gene.chr = seqnames,
        gene.start = start,
        gene.end = end,
        gene.strand = strand
      )
  ) |>
    mutate(diatance = dist)) |>
  write_tsv('1-pairs4ute.tsv')
