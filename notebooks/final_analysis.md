Analysis of prophages in 10 Klebsiella pneumoniae genomes
================
Eden Black

10 *Klebsiella pneumoniae* genomes isolated in Nanchang, China were
chosen for analysis. 6 of the genomes are from the dominant ST11 strain
(which accounts for 60% of infections in China (SOURCE)) and one each of
the other four strain types with NCBI entries for isolates from
Nanchang. With one genome per lineage for these STs, statistical testing
is not plausible - this is instead a descriptive analysis.

``` r
samples <- read_tsv(here("results", "tables", "samples.tsv"))
```

    ## Rows: 10 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): sample
    ## dbl (5): n_prophages, n_prophage_cds, n_defense_systems, n_defense_in_propha...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
prophages <- read_tsv(here("results", "tables", "prophages.tsv"))
```

    ## Rows: 60 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): sample, seq_name, source_seq, integrases
    ## dbl (5): start, end, length, n_genes, v_vs_c_score
    ## lgl (1): in_seq_edge
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# First 1000 rows of CARD are N/A, so read_tsv assigns boolean type by default
# Fix with guess_max = Inf
prophage_genes <- read_tsv(here("results", "tables", "prophage_genes.tsv"),
                           guess_max = Inf)
```

    ## Rows: 2970 Columns: 44
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (29): sample, gene, strand, prophage, mmseqs_phrog, mmseqs_alnScore, mms...
    ## dbl  (8): start, stop, score, partial, CARD_alnScore, CARD_seqIdentity, CARD...
    ## lgl  (7): vfdb_hit, vfdb_alnScore, vfdb_seqIdentity, vfdb_eVal, vfdb_short_n...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
defence_systems <- read_tsv(here("results", "tables", "defense_systems.tsv"))
```

    ## Rows: 305 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (10): sample, sys_id, type, subtype, activity, sys_beg, sys_end, protein...
    ## dbl  (1): genes_count
    ## lgl  (1): in_prophage
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
defence_genes <- read_tsv(here("results", "tables", "defense_genes.tsv"))
```

    ## Rows: 497 Columns: 31
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (14): sample, replicon, hit_id, gene_name, model_fqn, sys_id, hit_gene_r...
    ## dbl (15): hit_pos, sys_loci, locus_num, sys_wholeness, sys_score, sys_occ, h...
    ## lgl  (2): counterpart, in_prophage
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# Treat strain as a factor and order by lineage
samples <- samples |>
  mutate(st = factor(paste0("ST", st))) |>
  arrange(st, sample) |>
  mutate(sample = fct_inorder(sample))

#Table for looking up strain of each sample
st_lookup <- select(samples, sample, st)
```

``` r
prophage_size <- prophages |>
  group_by(sample) |>
  summarise(prophage_bp = sum(length), .groups = "drop")

burden <- samples |>
  left_join(prophage_size, by = "sample") |>
  select(sample, st, n_prophages, n_prophage_cds, prophage_bp)

burden
```

<div class="kable-table">

| sample          | st    | n_prophages | n_prophage_cds | prophage_bp |
|:----------------|:------|------------:|---------------:|------------:|
| GCF_029167725.1 | ST11  |           7 |            318 |      278478 |
| GCF_038977565.1 | ST11  |           7 |            394 |      306941 |
| GCF_039515315.1 | ST11  |           9 |            428 |      443011 |
| GCF_041022025.1 | ST11  |           6 |            305 |      245672 |
| GCF_043950185.1 | ST11  |           7 |            363 |      283814 |
| GCF_045344295.1 | ST11  |           8 |            368 |      359701 |
| GCF_039515985.1 | ST15  |           2 |            105 |       76767 |
| GCF_004118955.1 | ST23  |           3 |            163 |      138080 |
| GCF_040029055.1 | ST485 |           5 |            272 |      209454 |
| GCF_002944845.1 | ST86  |           6 |            254 |      228070 |

</div>

``` r
burden_long <- burden |>
  select(sample, st, `Number of prophages` = n_prophages,
         `Prophage sequence length (kb)` = prophage_bp) |>
         pivot_longer(c(`Number of prophages`, `Prophage sequence length (kb)`))


p_burden <- ggplot(burden_long, aes(x = fct(sample), y = value, fill = st)) +
  geom_col() +
  facet_wrap(~name, scales = "free") +
  labs(title = "ST11 genomes have a higher prophage burden",
    subtitle = "Ten complate K. pneumoniae genomes; prophages called with geNomad",
    x = NULL, y = NULL, fill = "Sequence type"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.4, vjust = 0.6))

p_burden
```

<img src="final_analysis_files/figure-gfm/burden-chart-1.png" alt="" style="display: block; margin: auto;" />

The plots show an expected strong correlation between the number of
prophages in a bacterial genome and the sequence length of the
prophages. More interestingly, all isolates of the ST11 sequence type
contained as many or more integrated prophages as the ST86 isolate,
which contained the most prophages of any of the other sequence types.
ST11 is both highly virulent and frequently displays AMR characteristics
(SOURCE), both traits often encoded by prophage morons. To investigate,
the prophage gene information from pharokka was examined.

``` r
categories <- prophage_genes |>
  mutate(category = replace_na(category, "unknown function"),
         category = as.factor(category)) |>
  count(category, sort = TRUE) |>
  mutate(percent = round(100 * n / sum(n), 1))

categories
```

<div class="kable-table">

| category                                          |    n | percent |
|:--------------------------------------------------|-----:|--------:|
| unknown function                                  | 1304 |    43.9 |
| tail                                              |  381 |    12.8 |
| head and packaging                                |  304 |    10.2 |
| DNA, RNA and nucleotide metabolism                |  254 |     8.6 |
| transcription regulation                          |  203 |     6.8 |
| lysis                                             |  163 |     5.5 |
| connector                                         |  121 |     4.1 |
| integration and excision                          |  100 |     3.4 |
| other                                             |   92 |     3.1 |
| moron, auxiliary metabolic gene and host takeover |   48 |     1.6 |

</div>

``` r
p_categories <- ggplot(categories, aes(y = fct_rev(category), x = n,
                                       fill = category == "unknown function")) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(title = "43.9% of prophage genes have no known function",
       subtitle = sprintf("%d CDS across %d prophages; annotated with Pharokka (PHROG categories)",
                          nrow(prophage_genes), nrow(prophages)),
       y = NULL, x = "Number of CDS") +
  theme(axis.text.y = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 6))

p_categories
```

<img src="final_analysis_files/figure-gfm/categories-chart-1.png" alt="" style="display: block; margin: auto;" />
