library(tidyverse)

################################################################################
# tidy up CyanoCyc data

genes <-
  # combine "frame" id and human readible names
  # ! make sure the smart table export from the website has the exact same
  #   row order
  'cyanocyc/Gene-IDS.txt' |> 
  read_tsv()  |>
  select(gene = `All-Genes`) |>
  bind_cols(
    'cyanocyc/Gene-Name-from-All-genes-of-Synechocystis-sp.-PCC-6803-substr.-Kazusa.txt' |>
      read_tsv() |>
      select(
        seqnames = `Replicon of Gene`,
        start = `Left-End-Position`,
        end = `Right-End-Position`,
        strand = Direction,
        name = `All-Genes`,
        locus = `Accession-1`,
        locus_alt = `Accession-2`,
        synonyms = Synonyms,
        product = Product,
      )
  ) |>
  arrange(seqnames, start, end)

pathways <-
  'cyanocyc/Pathways-from-All-pathways-of-Synechocystis-sp.-PCC-6803-substr.-Kazusa.txt' |>
  read_tsv() |>
  select(pathid = Pathways, pathway = `Common-Name`, gene = `Genes of pathway`) |>
  separate_rows(gene, sep = ' // ') |>
  unique()


tus <-
  'cyanocyc/Transcription-Units-from-Transcription-Units-from-All-transcription-units-of-Synechocystis-sp.-PCC.txt' |>
  read_tsv() |>
  select(
    tu = `Transcription-Units`,
    gene =`Genes of transcription unit`  
  ) |>
  separate_rows(gene, sep = ' // ') |>
  unique() |>
  mutate_at('gene', str_replace, '^ncl', 'NCL') |>
  mutate_at('gene', str_replace, '-as', '-AS')

# assure that gene id is unique
assertthat::are_equal(
  genes |>
    count(gene) |>
    filter(n > 1) |>
    nrow(),
  0
)

# assure that locus_alt matches all genes
assertthat::are_equal(
  pathways |>
    anti_join(genes, 'gene') |>
    nrow(),
  0
)
# assure that locus_alt matches all genes
assertthat::are_equal(
  tus |>
    anti_join(genes, 'gene') |>
    nrow(),
  0
)

write_tsv(genes, '2-genes.tsv')
write_tsv(pathways, '2-pathways.tsv')
write_tsv(tus, '2-tus.tsv')

################################################################################

guides <-
  '1-seq-search.tsv.gz' |>
  read_tsv()

fit <-
  'raw-data/fitness.tsv' |>
  read_tsv()

################################################################################
# export guide fitnesses per condition

dat <-
  guides |>
  mutate_at('chr', str_remove, ' .*$') |>
  left_join(
    fit |>
      select(guide = sgRNA, condition, fitness) |>
      unique(),
    'guide'
  ) |>
  transmute(
    chr,
    start = start - 1,
    end,
    name = paste(guide, seq, sep = '-'),
    fitness,
    strand,
    condition
  ) |>
  arrange(chr, start, end, strand)

assertthat::are_equal(
  dat |>
    filter(start > end) |>
    nrow(),
  0
)

# bed of general guide positions
dat |>
  select(- condition) |>
  mutate(fitness = 0) |>
  unique() |>
  write_tsv('2-guides.bed')

dat |>
  pull(condition) |>
  unique() |>
  map(function(i) {
  })

################################################################################
# Build some sort of wig file

wig.dat <-
  dat |>
  select(chr, start, end, fitness, condition) |>
  group_by_all() |>
  reframe(i = start:end) |>
  ungroup() |>
  arrange(condition, chr, i) |>
  group_by(chr, condition, i) |>
  summarize(f = fitness[which.max(abs(fitness))]) |>
  ungroup()


wig.dat |>
  group_by(condition) |>
  group_split() |>
  map(function(genome) {
    co2 <- first(genome$condition)
    genome |>
      group_by(chr) |>
      group_split() |>
      map(function(chr) {
        c(
          sprintf('variableStep\tchrom=%s',
                  first(chr$chr)),
          with(chr, paste(i, f, sep = '\t'))
        )
      }) |>
      unlist() |>
      write_lines(sprintf('2-cripri-%s.wig', co2))
  })

################################################################################

genome.code <- c(                    # included in RefSeq ref
  'chromosome' =	'NC_000911.1',     # yes
  'plasmid pSYSM' =	'NC_005229.1',   # yes
  'plasmid pSYSX' = 'NC_005232.1',   # yes
  'plasmid pSYSA' =	'NC_005230.1',   # yes
  'plasmid pSYSG' =	'NC_005231.1',   # yes
  # 'plasmid pCC5.2' = 'NC_020290.1',  # (only 5 genes) \
  # 'plasmid pCA2.4' = 'NC_020289.1',  # (only 3 genes)  + - not targeted by CRISPRi
  # 'plasmid pCB2.4' = 'NZ_L25424.1'   # (only 2 genes) /
  # -> Exlcude
)

genes |>
  transmute(
    chr = genome.code[seqnames],
    start = start,
    end,
    name = sprintf(
      '%s %s (%s)',
      locus, name, product
    ),
    score = 0,
    strand
  ) |>
  # exclude small plasmids
  drop_na(chr) |>
  mutate(
    tmp = start,
    start = ifelse(start > end, end, start) - 1,
    end = ifelse(tmp > end, tmp, end)
  ) |>
  select(-tmp) |>
  drop_na() |>
  arrange(chr, start, end, strand) |>
  write_tsv('2-genes.bed', col_names = FALSE)

tus |>
  left_join(genes) |>
  select(tu, gene, seqnames, strand, start, end) |>
  group_by(tu, seqnames) |>
  # ignore strand, since not fully consistent
  summarize(
    start = min(c(start, end)),
    end = max(c(start, end))
  ) |>
  ungroup() |>
  mutate(chr = genome.code[seqnames]) |>
  drop_na(chr) |>
  transmute(
    chr,
    start = start - 1,
    end,
    name = tu,
    score = 0
  ) |>
  write_tsv('2-tu.bed')

################################################################################

# try to match genes by any synonym
# name.list <-
#   genes |>
#   select(gene, name, locus, locus_alt, synonyms) |>
#   separate_rows(synonyms, sep = ' // ') |>
#   pivot_longer(- gene) |>
#   select(gene, locus = value) |>
#   drop_na()

rnas <-
  'cyanocyc/All-instances-of-RNAs-in-Synechocystis-sp.-PCC-6803-substr.-Kazusa.txt' |>
  read_tsv() |>
  select(gene = Gene) |>
  mutate(type = 'non-coding')

fit.dat <-
  fit |>
  select(sgRNA_target, condition, wmean_fitness) |>
  unique() |>
  spread(condition, wmean_fitness) |>
  mutate_at('sgRNA_target', str_remove, 'nc_') |>
  left_join(
    'raw-data/names.tsv' |>
      read_tsv(),
    c('sgRNA_target' = 'gene')
  ) |>
  # fill in, non-coding are not part of list ?
  mutate(locus = ifelse(is.na(locus), sgRNA_target, locus)) |>
  # check which matches fail
  # anti_join(genes, 'locus') |>
  # pull(locus) |> unique()
# [1] "ctrl1"    "ctrl10"   "ctrl2"    "ctrl3"    "ctrl4"    "ctrl5"    "ctrl6"    "ctrl7"    "ctrl8"    "ctrl9"    "sll1572"  "slr0603"  "ncl1680" 
# [14] "ncr0700"  "sll0375"  "slr0163"  "slr0366"  "slr0743a" "ssr2317"  "smr0010"  "smr0003"  "smr0006"  "sll0752"  "sll0763"  "sll0811"  "sll1157" 
# [27] "sll1250"  "sll1575"  "sll1792"  "sll5072"  "sll6052"  "slr0800"  "slr1246"  "slr1283"  "slr1586"  "slr1903"  "slr2096"  "slr7076"  "ssl1507" 
# [40] "ssl1922"  "ssl5045"  "ssl5068"  "ssl5108"  "ssl7019"  "ssl7020"  "ssl7021"  "ssl7022"  "ssl8010"  "ssl8039"  "ssl8041"  "ssr2708"  "ssr2898" 
# [53] "ssr6032"  "ssr7018"  "ssr7079"  "ssr8047"  "sml0007" 
# contorl genes and some weired extra ones -> ignore
  inner_join(genes, 'locus') |>
  mutate(chr = genome.code[seqnames]) |>
  drop_na(chr) |>
  left_join(rnas, 'gene') |>
  select(chr, start, end, strand, type, gene, name, product, locus, CO2_30percent, CO2_4percent) |>
  arrange(chr, start, end) |>
  mutate_at('type', replace_na, 'protein-coding')

write_tsv(fit.dat, '2-crispri.tsv')


################################################################################

fit.dat |>
  ggplot(aes(CO2_30percent, CO2_4percent, color = type)) +
  geom_point()

