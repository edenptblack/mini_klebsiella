Selection of 10 Klebsiella pneumoniae genomes for analysis
================
Eden Black

``` r
library(tidyverse)
library(here)

here::i_am("notebooks/candidate_selection.Rmd")
```

candidates.tsv was assembled from the NCBI database searching for
complete Klebsiella pneumoniae genomes from RefSeq using the following
code in the NCBI CLI:

datasets summary genome taxon “Klebsiella pneumoniae”\
–assembly-level complete –assembly-source refseq\
–as-json-lines\
\| dataformat tsv genome –fields\
accession,source_database,organism-name,organism-infraspecific-strain,\
assminfo-level,assmstats-total-sequence-len,assmstats-number-of-contigs,\
assmstats-total-number-of-chromosomes,assmstats-gc-percent,\
checkm-completeness,checkm-contamination,assminfo-sequencing-tech,\
assminfo-biosample-accession,assminfo-biosample-isolation-source,\
assminfo-biosample-host,assminfo-biosample-host-disease,\
assminfo-biosample-geo-loc-name,assminfo-biosample-collection-date,\
assminfo-bioproject,assminfo-bioproject-lineage-title,assminfo-release-date\
\> data/candidates.tsv

As a case study for this minimal analysis, the aim is to analyse the
genomes of different samples from a single geographical location to
identify similarities and differences, particularly intra- and
inter-strain.

``` r
candidates = read_tsv(here('data', 'candidates.tsv')) %>%
    mutate(`Assembly BioSample Accession` = as.factor(`Assembly BioSample Accession`))
```

    ## Rows: 3752 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (14): Assembly Accession, Source Database, Organism Name, Organism Infr...
    ## dbl   (6): Assembly Stats Total Sequence Length, Assembly Stats Number of Co...
    ## date  (1): Assembly Release Date
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
colnames(candidates)
```

    ##  [1] "Assembly Accession"                        
    ##  [2] "Source Database"                           
    ##  [3] "Organism Name"                             
    ##  [4] "Organism Infraspecific Names Strain"       
    ##  [5] "Assembly Level"                            
    ##  [6] "Assembly Stats Total Sequence Length"      
    ##  [7] "Assembly Stats Number of Contigs"          
    ##  [8] "Assembly Stats Total Number of Chromosomes"
    ##  [9] "Assembly Stats GC Percent"                 
    ## [10] "CheckM completeness"                       
    ## [11] "CheckM contamination"                      
    ## [12] "Assembly Sequencing Tech"                  
    ## [13] "Assembly BioSample Accession"              
    ## [14] "Assembly BioSample Isolation source"       
    ## [15] "Assembly BioSample Host"                   
    ## [16] "Assembly BioSample Host disease"           
    ## [17] "Assembly BioSample Geographic location"    
    ## [18] "Assembly BioSample Collection date"        
    ## [19] "Assembly BioProject Accession"             
    ## [20] "Assembly BioProject Lineage Title"         
    ## [21] "Assembly Release Date"

``` r
candidate_locations = candidates %>%
    group_by(`Assembly BioSample Geographic location`) %>%
    summarise(count = n()) %>%
    filter(count >= 10) %>%
    arrange(desc(count))

head(candidate_locations, 10)
```

<div class="kable-table">

| Assembly BioSample Geographic location | count |
|:---------------------------------------|------:|
| Taiwan                                 |   334 |
| China                                  |   181 |
| missing                                |   146 |
| Norway                                 |   138 |
| USA                                    |   108 |
| Australia                              |   101 |
| Switzerland                            |    88 |
| South Korea: Seoul                     |    85 |
| China:Nanchang                         |    83 |
| Germany                                |    57 |

</div>

While many samples have their locations logged only by country, two
cities appear on the top ten list of locations: Seoul (85) and Nanchang
(83).

``` r
seoul_candidates = candidates %>%
    filter(`Assembly BioSample Geographic location` == 'South Korea: Seoul' &
    `Assembly Stats Number of Contigs` <= 10 &
    `CheckM completeness` >= 98 &
    `CheckM contamination` < 1)

nrow(seoul_candidates)
```

    ## [1] 10

``` r
seoul_candidates %>% count(`Assembly BioProject Accession`, sort = TRUE)
```

<div class="kable-table">

| Assembly BioProject Accession |   n |
|:------------------------------|----:|
| PRJNA1230115                  |   4 |
| PRJNA396106                   |   3 |
| PRJNA428137                   |   1 |
| PRJNA594105                   |   1 |
| PRJNA594107                   |   1 |

</div>

``` r
sprintf('For samples from Seoul, the mean completeness is %.1f and mean contamination is %.1f',
    mean(seoul_candidates$`CheckM completeness`),
    mean(seoul_candidates$`CheckM contamination`))
```

    ## [1] "For samples from Seoul, the mean completeness is 98.6 and mean contamination is 0.6"

``` r
nanchang_candidates = candidates %>%
    filter(`Assembly BioSample Geographic location` == 'China:Nanchang' &
    `Assembly Stats Number of Contigs` <= 10 &
    `CheckM completeness` >= 98 &
    `CheckM contamination` < 1)

nrow(nanchang_candidates)
```

    ## [1] 71

``` r
nanchang_candidates %>% count(`Assembly BioProject Accession`, sort = TRUE)
```

<div class="kable-table">

| Assembly BioProject Accession |   n |
|:------------------------------|----:|
| PRJNA1474683                  |  27 |
| PRJNA672246                   |  11 |
| PRJNA1469284                  |   3 |
| PRJNA1344441                  |   2 |
| PRJNA1104315                  |   1 |
| PRJNA1104368                  |   1 |
| PRJNA1104369                  |   1 |
| PRJNA1104370                  |   1 |
| PRJNA1104371                  |   1 |
| PRJNA1105644                  |   1 |
| PRJNA1105645                  |   1 |
| PRJNA1105646                  |   1 |
| PRJNA1105649                  |   1 |
| PRJNA1116330                  |   1 |
| PRJNA1116332                  |   1 |
| PRJNA1139232                  |   1 |
| PRJNA1139245                  |   1 |
| PRJNA1139252                  |   1 |
| PRJNA1139253                  |   1 |
| PRJNA1139309                  |   1 |
| PRJNA1139310                  |   1 |
| PRJNA1139315                  |   1 |
| PRJNA1139322                  |   1 |
| PRJNA1153625                  |   1 |
| PRJNA1153628                  |   1 |
| PRJNA1175394                  |   1 |
| PRJNA1187795                  |   1 |
| PRJNA1187798                  |   1 |
| PRJNA1187800                  |   1 |
| PRJNA432625                   |   1 |
| PRJNA516531                   |   1 |
| PRJNA940206                   |   1 |

</div>

``` r
sprintf('For samples from Nanchang, the mean completeness is %.1f and mean contamination is %.1f',
    mean(nanchang_candidates$`CheckM completeness`),
    mean(nanchang_candidates$`CheckM contamination`))
```

    ## [1] "For samples from Nanchang, the mean completeness is 98.8 and mean contamination is 0.5"

71 candidates from Nanchang passed QC, compared to only 10 from Seoul.
As such, the samples from Nanchang were taken forward.

The final step for candidate selection was sample independence. Genomes
submitted in the same project are more likely to be the same or very
closely-related isolate sampled repeatedly. As 28 samples were submitted
as individuals, they were selected to be taken forward for multilocus
sequence typing.

``` r
nanchang_candidates = nanchang_candidates %>%
    filter(!(`Assembly BioProject Accession` %in% c('PRJNA1474683', 'PRJNA672246',
    'PRJNA1469284', 'PRJNA1344441')))

nrow(nanchang_candidates)
```

    ## [1] 28

Finally, the accession numbers of the candidate genomes were exported to
a .txt file to download from NCBI.

``` r
write_lines(nanchang_candidates$`Assembly Accession`, here('data', 'accession.txt'))
```
