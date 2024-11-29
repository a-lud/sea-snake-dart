## Load metadata and VCF files

Here we’re simply loading the population-map files and storing them into
the object `meta`. This object is a list.

``` r
# Load data
meta <- fs::dir_ls(
    here("data", "popmaps"),
    glob = "*.txt"
) |>
    (\(x) set_names(x, str_remove(basename(x), "-.*")))() |>
    map(
        read_delim, 
        delim = " ", col_names = c("sample", "population"), col_types = cols()
    )
```

Next we load the VCF files into `R` and `snpRdata` objects. Similarly,
we load each species file separately and store it in a list. Some
filtering is applied which somewhat mirrors what we’ve already already
applied. I’m including it here as another approach to filtering a VCF
file.

``` r
# List VCF files and name vector elements
vcfs <- fs::dir_ls(
    here("results", "ipyrad"),
    glob = "*stringent.vcf.gz",
    recurse = TRUE
) |>
    (\(x) set_names(x, str_remove(basename(x), "-stringent.*")))()

# Load each VCF file separately (with corresponding metadata)
vcfs <- imap(vcfs, \(v, i) {
    tmp_vcf <- read_vcf(
        file = v, 
        sample.meta = meta[[i]]
    )
    
    filter_snps(
        tmp_vcf,
        remove_garbage = 0.2,
        min_ind = 0.8,
        min_loci = 0.8,
        re_run = TRUE
    )
})
```

## Calculate statistics

Next we’ll calculate a range of diversity statistics. The program `snpR`
has a really convenient function `calc_basic_snp_stats()` which
practically estimates all the statistics we want.

Below we estimate maf, pi, ho, he, pairwise Fst, Fis, HWE in addition to
`Fst`. We estimate all statistics by population.

``` r
# Run calc_basic_snp_stats on each VCF object
basic_stats <- map(
    vcfs, 
    calc_basic_snp_stats, 
    # Below are aguments to 'calc_basic_snp_stats'
    facets = "population", fst.method = "WC"
)

# Get the statistics as list of tables
basic_stats_list_tables <- map(
    basic_stats, 
    get.snpR.stats,
    facets = "population", 
    stats = c("pi", "hwe", "ho", "he", "fis", "fst") 
)
```

## External calculation of *π*: Why invariant data matters

Nucleotide diversity calculated across multiple sites is given by the
equation below.

$$
\pi = \frac{1}{L} \sum\_{i=1}^{L} \left( \frac{2n_i(n_i-1)}{N(N-1)} \right)
$$

An important component of this equation is $\frac{1}{L}$. This is the
reciprocal of the **total** sequence length/number of sites being
considered. We know that most sites in the genome don’t vary - e.g. in
humans we have approximately **one SNP per kilobase (1kb)**. Below is a
schematic showing SNPs along a sequence of length 10.

    5' -  1 2 3 4 5 6 7 8 9 10 - 3'  Position
    *         *
    5' -  A T A C C G T A A C  - 3'  Seq1
    5' -  A T A G C G T A A C  - 3'  Seq2
    5' -  A T A G C G T A C C  - 3'  Seq3
    5' -  A T A C C G T A A C  - 3'  Seq4

Using the example above, if we considered **only** SNP positions, the
length would equal 2 (*L* = 2). This is not the actual sequence length,
but if we only provided variable sites to the equation, the program
calculating *π* wont know that and will simply take the number of
variants as the length. As the number of variant sites is **much
shorter** than the actual length, the value of *π* is scaled
incorrectly, resulting in a higher estimate due to a small value of *L*.

$$
\pi = \frac{1}{2} \times \frac{7}{6} \approx 0.5833
$$

-   Length (*L*) = 2 (2 variant sites)
-   Total number of differences (∑d<sub>*i**j*</sub>) = 7
-   Total number of comparisons $(\frac{4}{2})$ = 6 (1..2, 1..3, 1..4,
    2..3 etc…)

Instead, we should calculate length (*L*) using both variant **and**
invariant sites. The variant sites will still be used to calculate
divergence, but the combined length is used to properly scale the level
of differentiation. Using the correct length (*L* = 10) would scale the
value of *π* down, resulting in a more accurate reflection of divergence
within the population.

$$
\pi = \frac{1}{10} \times \frac{7}{6} \approx 0.1167
$$

The program `pixy` calculates *π* the second way when given a VCF with
both variant and invariant data. The authors of `Ipyrad` wrote a script
that can return an **all-sites** VCF file from the `.loci` files. I’ve
passed that to `pixy` to get a more accurate estimation of *π* in each
population.

``` r
pixy_pi <- fs::dir_ls(
    here("results"),
    glob = "*pi.txt",
    recurse = TRUE
) |>
    read_tsv(col_types = cols(), id = "species") |>
    mutate(species = str_remove(basename(species), "_.*")) |>
    summarise(
        pixy_pi = sum(count_diffs)/sum(count_comparisons),
        .by = c(species, pop)
    ) |>
    rename(subfacet = pop)
```

The *π* estimates are included in the next table.

## Within population statistics

Below we simply extract the weighted averages for each population. The
table contains four relevant columns - weighted means for nucleotide
diversity, observed heterozygosity, expected heterozygosity and `Fis`.
I’ve appended the `pixy` *π* estimates so you can compare the impact of
calculating nucleotide diversity with and without invariant sites.

``` r
l_ala <- length(unique(meta$ALA$population))
l_hma <- length(unique(meta$HMA$population))
l_hst <- length(unique(meta$HST$population))

within <- basic_stats_list_tables |>
    imap(
        \(x, i) {
            x |> 
                pluck("weighted.means") |> 
                as_tibble() |> 
                filter(is.na(weighted_mean_fst)) |>
                select(-c(snp.facet, snp.subfacet, weighted_mean_fst, mean_fst)) |>
                left_join(
                    meta[[i]] |> summarise(n_samples = n(), .by = population),
                    by = join_by(subfacet == population)
                )
            
        }) |>
    list_rbind(names_to = "species") |>
    left_join(pixy_pi)

# Write to file
within |>
    (\(x) split(x, x$species))() |>
    iwalk(\(x, i) {
        x |>
            select(
                species, 
                population = subfacet, 
                n_samples, 
                pixy_pi, 
                everything(),
                -facet
            ) |>
            write_csv(
                file = here("results", "diversity-statistics", glue::glue("{i}-within-stats.csv"))
            )
    })

# Example of output
within |>
    select(subfacet, n_samples, pixy_pi, everything(), -facet) |>
    gt(groupname_col = "species") |>
    fmt_number(
        columns = c(starts_with("weighted_"), "pixy_pi"), 
        n_sigfig = 3
    ) |>
    cols_label(
        subfacet = ""
    )
```

<div id="dvfzapqyhu" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#dvfzapqyhu table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#dvfzapqyhu thead, #dvfzapqyhu tbody, #dvfzapqyhu tfoot, #dvfzapqyhu tr, #dvfzapqyhu td, #dvfzapqyhu th {
  border-style: none;
}

#dvfzapqyhu p {
  margin: 0;
  padding: 0;
}

#dvfzapqyhu .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#dvfzapqyhu .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#dvfzapqyhu .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#dvfzapqyhu .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#dvfzapqyhu .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#dvfzapqyhu .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#dvfzapqyhu .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#dvfzapqyhu .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#dvfzapqyhu .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#dvfzapqyhu .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#dvfzapqyhu .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#dvfzapqyhu .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#dvfzapqyhu .gt_spanner_row {
  border-bottom-style: hidden;
}

#dvfzapqyhu .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#dvfzapqyhu .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#dvfzapqyhu .gt_from_md > :first-child {
  margin-top: 0;
}

#dvfzapqyhu .gt_from_md > :last-child {
  margin-bottom: 0;
}

#dvfzapqyhu .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#dvfzapqyhu .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#dvfzapqyhu .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#dvfzapqyhu .gt_row_group_first td {
  border-top-width: 2px;
}

#dvfzapqyhu .gt_row_group_first th {
  border-top-width: 2px;
}

#dvfzapqyhu .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#dvfzapqyhu .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#dvfzapqyhu .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#dvfzapqyhu .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#dvfzapqyhu .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#dvfzapqyhu .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#dvfzapqyhu .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#dvfzapqyhu .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#dvfzapqyhu .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#dvfzapqyhu .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#dvfzapqyhu .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#dvfzapqyhu .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#dvfzapqyhu .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#dvfzapqyhu .gt_left {
  text-align: left;
}

#dvfzapqyhu .gt_center {
  text-align: center;
}

#dvfzapqyhu .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#dvfzapqyhu .gt_font_normal {
  font-weight: normal;
}

#dvfzapqyhu .gt_font_bold {
  font-weight: bold;
}

#dvfzapqyhu .gt_font_italic {
  font-style: italic;
}

#dvfzapqyhu .gt_super {
  font-size: 65%;
}

#dvfzapqyhu .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#dvfzapqyhu .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#dvfzapqyhu .gt_indent_1 {
  text-indent: 5px;
}

#dvfzapqyhu .gt_indent_2 {
  text-indent: 10px;
}

#dvfzapqyhu .gt_indent_3 {
  text-indent: 15px;
}

#dvfzapqyhu .gt_indent_4 {
  text-indent: 20px;
}

#dvfzapqyhu .gt_indent_5 {
  text-indent: 25px;
}

#dvfzapqyhu .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#dvfzapqyhu div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id=""></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n_samples">n_samples</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="pixy_pi">pixy_pi</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="weighted_mean_pi">weighted_mean_pi</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="weighted_mean_ho">weighted_mean_ho</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="weighted_mean_he">weighted_mean_he</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="weighted_mean_fis">weighted_mean_fis</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="7" class="gt_group_heading" scope="colgroup" id="ALA">ALA</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="ALA  subfacet" class="gt_row gt_left">Ashmore</td>
<td headers="ALA  n_samples" class="gt_row gt_right">87</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00215</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0470</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0433</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0467</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0.0781</td></tr>
    <tr><td headers="ALA  subfacet" class="gt_row gt_left">Broome</td>
<td headers="ALA  n_samples" class="gt_row gt_right">28</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00193</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0492</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0448</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0482</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0.0913</td></tr>
    <tr><td headers="ALA  subfacet" class="gt_row gt_left">Exmouth_Gulf</td>
<td headers="ALA  n_samples" class="gt_row gt_right">14</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00192</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0489</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0453</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0471</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0.0757</td></tr>
    <tr><td headers="ALA  subfacet" class="gt_row gt_left">Gulf_of_Carpentaria</td>
<td headers="ALA  n_samples" class="gt_row gt_right">16</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00192</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0437</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0375</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0396</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0.159</td></tr>
    <tr><td headers="ALA  subfacet" class="gt_row gt_left">Heywood_Shoal</td>
<td headers="ALA  n_samples" class="gt_row gt_right">3</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00174</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0453</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0435</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0333</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0.0596</td></tr>
    <tr><td headers="ALA  subfacet" class="gt_row gt_left">New_Caledonia</td>
<td headers="ALA  n_samples" class="gt_row gt_right">12</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00122</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0324</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0317</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0308</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0.0222</td></tr>
    <tr><td headers="ALA  subfacet" class="gt_row gt_left">Pilbara</td>
<td headers="ALA  n_samples" class="gt_row gt_right">30</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00187</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0487</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0453</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0478</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0.0722</td></tr>
    <tr><td headers="ALA  subfacet" class="gt_row gt_left">Scott_Reef</td>
<td headers="ALA  n_samples" class="gt_row gt_right">1</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00143</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0415</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0415</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0208</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0</td></tr>
    <tr><td headers="ALA  subfacet" class="gt_row gt_left">Shark_Bay</td>
<td headers="ALA  n_samples" class="gt_row gt_right">10</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00181</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0470</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0435</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0445</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0.0772</td></tr>
    <tr><td headers="ALA  subfacet" class="gt_row gt_left">Townsville</td>
<td headers="ALA  n_samples" class="gt_row gt_right">4</td>
<td headers="ALA  pixy_pi" class="gt_row gt_right">0.00164</td>
<td headers="ALA  weighted_mean_pi" class="gt_row gt_right">0.0405</td>
<td headers="ALA  weighted_mean_ho" class="gt_row gt_right">0.0377</td>
<td headers="ALA  weighted_mean_he" class="gt_row gt_right">0.0350</td>
<td headers="ALA  weighted_mean_fis" class="gt_row gt_right">0.0827</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="7" class="gt_group_heading" scope="colgroup" id="HMA">HMA</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="HMA  subfacet" class="gt_row gt_left">Broome</td>
<td headers="HMA  n_samples" class="gt_row gt_right">16</td>
<td headers="HMA  pixy_pi" class="gt_row gt_right">0.00155</td>
<td headers="HMA  weighted_mean_pi" class="gt_row gt_right">0.0810</td>
<td headers="HMA  weighted_mean_ho" class="gt_row gt_right">0.0765</td>
<td headers="HMA  weighted_mean_he" class="gt_row gt_right">0.0783</td>
<td headers="HMA  weighted_mean_fis" class="gt_row gt_right">0.0578</td></tr>
    <tr><td headers="HMA  subfacet" class="gt_row gt_left">Exmouth_Gulf</td>
<td headers="HMA  n_samples" class="gt_row gt_right">8</td>
<td headers="HMA  pixy_pi" class="gt_row gt_right">0.00148</td>
<td headers="HMA  weighted_mean_pi" class="gt_row gt_right">0.0819</td>
<td headers="HMA  weighted_mean_ho" class="gt_row gt_right">0.0771</td>
<td headers="HMA  weighted_mean_he" class="gt_row gt_right">0.0767</td>
<td headers="HMA  weighted_mean_fis" class="gt_row gt_right">0.0635</td></tr>
    <tr><td headers="HMA  subfacet" class="gt_row gt_left">Fraser_Island</td>
<td headers="HMA  n_samples" class="gt_row gt_right">4</td>
<td headers="HMA  pixy_pi" class="gt_row gt_right">0.00125</td>
<td headers="HMA  weighted_mean_pi" class="gt_row gt_right">0.0778</td>
<td headers="HMA  weighted_mean_ho" class="gt_row gt_right">0.0753</td>
<td headers="HMA  weighted_mean_he" class="gt_row gt_right">0.0678</td>
<td headers="HMA  weighted_mean_fis" class="gt_row gt_right">0.0374</td></tr>
    <tr><td headers="HMA  subfacet" class="gt_row gt_left">Gulf_of_Carpentaria</td>
<td headers="HMA  n_samples" class="gt_row gt_right">29</td>
<td headers="HMA  pixy_pi" class="gt_row gt_right">0.00166</td>
<td headers="HMA  weighted_mean_pi" class="gt_row gt_right">0.0819</td>
<td headers="HMA  weighted_mean_ho" class="gt_row gt_right">0.0744</td>
<td headers="HMA  weighted_mean_he" class="gt_row gt_right">0.0801</td>
<td headers="HMA  weighted_mean_fis" class="gt_row gt_right">0.0938</td></tr>
    <tr><td headers="HMA  subfacet" class="gt_row gt_left">New_Caledonia</td>
<td headers="HMA  n_samples" class="gt_row gt_right">10</td>
<td headers="HMA  pixy_pi" class="gt_row gt_right">0.000916</td>
<td headers="HMA  weighted_mean_pi" class="gt_row gt_right">0.0555</td>
<td headers="HMA  weighted_mean_ho" class="gt_row gt_right">0.0497</td>
<td headers="HMA  weighted_mean_he" class="gt_row gt_right">0.0526</td>
<td headers="HMA  weighted_mean_fis" class="gt_row gt_right">0.111</td></tr>
    <tr><td headers="HMA  subfacet" class="gt_row gt_left">Pilbara</td>
<td headers="HMA  n_samples" class="gt_row gt_right">12</td>
<td headers="HMA  pixy_pi" class="gt_row gt_right">0.00161</td>
<td headers="HMA  weighted_mean_pi" class="gt_row gt_right">0.0823</td>
<td headers="HMA  weighted_mean_ho" class="gt_row gt_right">0.0719</td>
<td headers="HMA  weighted_mean_he" class="gt_row gt_right">0.0774</td>
<td headers="HMA  weighted_mean_fis" class="gt_row gt_right">0.133</td></tr>
    <tr><td headers="HMA  subfacet" class="gt_row gt_left">Shark_Bay</td>
<td headers="HMA  n_samples" class="gt_row gt_right">12</td>
<td headers="HMA  pixy_pi" class="gt_row gt_right">0.00155</td>
<td headers="HMA  weighted_mean_pi" class="gt_row gt_right">0.0831</td>
<td headers="HMA  weighted_mean_ho" class="gt_row gt_right">0.0669</td>
<td headers="HMA  weighted_mean_he" class="gt_row gt_right">0.0773</td>
<td headers="HMA  weighted_mean_fis" class="gt_row gt_right">0.208</td></tr>
    <tr><td headers="HMA  subfacet" class="gt_row gt_left">Weipa</td>
<td headers="HMA  n_samples" class="gt_row gt_right">1</td>
<td headers="HMA  pixy_pi" class="gt_row gt_right">0.00125</td>
<td headers="HMA  weighted_mean_pi" class="gt_row gt_right">0.0726</td>
<td headers="HMA  weighted_mean_ho" class="gt_row gt_right">0.0726</td>
<td headers="HMA  weighted_mean_he" class="gt_row gt_right">0.0363</td>
<td headers="HMA  weighted_mean_fis" class="gt_row gt_right">0</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="7" class="gt_group_heading" scope="colgroup" id="HST">HST</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="HST  subfacet" class="gt_row gt_left">Ashmore</td>
<td headers="HST  n_samples" class="gt_row gt_right">3</td>
<td headers="HST  pixy_pi" class="gt_row gt_right">0.00168</td>
<td headers="HST  weighted_mean_pi" class="gt_row gt_right">0.0472</td>
<td headers="HST  weighted_mean_ho" class="gt_row gt_right">0.0436</td>
<td headers="HST  weighted_mean_he" class="gt_row gt_right">0.0346</td>
<td headers="HST  weighted_mean_fis" class="gt_row gt_right">0.109</td></tr>
    <tr><td headers="HST  subfacet" class="gt_row gt_left">Broome</td>
<td headers="HST  n_samples" class="gt_row gt_right">64</td>
<td headers="HST  pixy_pi" class="gt_row gt_right">0.00152</td>
<td headers="HST  weighted_mean_pi" class="gt_row gt_right">0.0729</td>
<td headers="HST  weighted_mean_ho" class="gt_row gt_right">0.0472</td>
<td headers="HST  weighted_mean_he" class="gt_row gt_right">0.0722</td>
<td headers="HST  weighted_mean_fis" class="gt_row gt_right">0.355</td></tr>
    <tr><td headers="HST  subfacet" class="gt_row gt_left">Exmouth_Gulf</td>
<td headers="HST  n_samples" class="gt_row gt_right">17</td>
<td headers="HST  pixy_pi" class="gt_row gt_right">0.00199</td>
<td headers="HST  weighted_mean_pi" class="gt_row gt_right">0.0494</td>
<td headers="HST  weighted_mean_ho" class="gt_row gt_right">0.0471</td>
<td headers="HST  weighted_mean_he" class="gt_row gt_right">0.0479</td>
<td headers="HST  weighted_mean_fis" class="gt_row gt_right">0.0470</td></tr>
    <tr><td headers="HST  subfacet" class="gt_row gt_left">Gulf_of_Carpentaria</td>
<td headers="HST  n_samples" class="gt_row gt_right">8</td>
<td headers="HST  pixy_pi" class="gt_row gt_right">0.00113</td>
<td headers="HST  weighted_mean_pi" class="gt_row gt_right">0.0471</td>
<td headers="HST  weighted_mean_ho" class="gt_row gt_right">0.0403</td>
<td headers="HST  weighted_mean_he" class="gt_row gt_right">0.0384</td>
<td headers="HST  weighted_mean_fis" class="gt_row gt_right">0.179</td></tr>
    <tr><td headers="HST  subfacet" class="gt_row gt_left">Heywood_Shoal</td>
<td headers="HST  n_samples" class="gt_row gt_right">4</td>
<td headers="HST  pixy_pi" class="gt_row gt_right">0.000917</td>
<td headers="HST  weighted_mean_pi" class="gt_row gt_right">0.0418</td>
<td headers="HST  weighted_mean_ho" class="gt_row gt_right">0.0398</td>
<td headers="HST  weighted_mean_he" class="gt_row gt_right">0.0341</td>
<td headers="HST  weighted_mean_fis" class="gt_row gt_right">0.0598</td></tr>
    <tr><td headers="HST  subfacet" class="gt_row gt_left">North_Kimberley</td>
<td headers="HST  n_samples" class="gt_row gt_right">1</td>
<td headers="HST  pixy_pi" class="gt_row gt_right">0.000930</td>
<td headers="HST  weighted_mean_pi" class="gt_row gt_right">0.0467</td>
<td headers="HST  weighted_mean_ho" class="gt_row gt_right">0.0467</td>
<td headers="HST  weighted_mean_he" class="gt_row gt_right">0.0234</td>
<td headers="HST  weighted_mean_fis" class="gt_row gt_right">0</td></tr>
    <tr><td headers="HST  subfacet" class="gt_row gt_left">North_QLD</td>
<td headers="HST  n_samples" class="gt_row gt_right">2</td>
<td headers="HST  pixy_pi" class="gt_row gt_right">0.000884</td>
<td headers="HST  weighted_mean_pi" class="gt_row gt_right">0.0480</td>
<td headers="HST  weighted_mean_ho" class="gt_row gt_right">0.0378</td>
<td headers="HST  weighted_mean_he" class="gt_row gt_right">0.0354</td>
<td headers="HST  weighted_mean_fis" class="gt_row gt_right">0.288</td></tr>
    <tr><td headers="HST  subfacet" class="gt_row gt_left">Pilbara</td>
<td headers="HST  n_samples" class="gt_row gt_right">15</td>
<td headers="HST  pixy_pi" class="gt_row gt_right">0.00161</td>
<td headers="HST  weighted_mean_pi" class="gt_row gt_right">0.0484</td>
<td headers="HST  weighted_mean_ho" class="gt_row gt_right">0.0457</td>
<td headers="HST  weighted_mean_he" class="gt_row gt_right">0.0465</td>
<td headers="HST  weighted_mean_fis" class="gt_row gt_right">0.0570</td></tr>
  </tbody>
  
  
</table>
</div>

I’ve left the *π* calculation based on **just SNPs** in the table. We
can see that the values are much higher than the `pixy` estimations.

# Between population statistics

Lastly we have our F<sub>ST</sub> measures. I’ve separated out the upper
matricies for each of the three species.

``` r
# Fst
# fst_table <- basic_stats_list_tables |>
#     map(
#         \(i) 
#         i |> 
#             pluck("weighted.means") |> 
#             as_tibble() |> 
#             filter(! is.na(weighted_mean_fst)) |>
#             select(facet, subfacet, weighted_mean_fst, mean_fst)
#     )

fst_upper_matrix <- basic_stats_list_tables |>
    map(
        \(i) 
        i |> 
            pluck("fst.matrix") |>
            pluck("population") |>
            as_tibble()
    )
```

``` r
# Write FST matricies to file
iwalk(fst_upper_matrix, \(x, i) {
    x |>
        write_csv(
            file = here("results", "diversity-statistics", glue::glue("{i}-fst.csv")),
            na = ""
        )
})
```

### *Aipysurus laevis*

``` r
fst_upper_matrix$ALA |>
    gt() |>
    fmt_missing(missing_text = "") |>
    cols_label(p1 = "") |>
    fmt_number(n_sigfig = 3)
```

<div id="olsgxmcxes" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#olsgxmcxes table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#olsgxmcxes thead, #olsgxmcxes tbody, #olsgxmcxes tfoot, #olsgxmcxes tr, #olsgxmcxes td, #olsgxmcxes th {
  border-style: none;
}

#olsgxmcxes p {
  margin: 0;
  padding: 0;
}

#olsgxmcxes .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#olsgxmcxes .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#olsgxmcxes .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#olsgxmcxes .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#olsgxmcxes .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#olsgxmcxes .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#olsgxmcxes .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#olsgxmcxes .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#olsgxmcxes .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#olsgxmcxes .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#olsgxmcxes .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#olsgxmcxes .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#olsgxmcxes .gt_spanner_row {
  border-bottom-style: hidden;
}

#olsgxmcxes .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#olsgxmcxes .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#olsgxmcxes .gt_from_md > :first-child {
  margin-top: 0;
}

#olsgxmcxes .gt_from_md > :last-child {
  margin-bottom: 0;
}

#olsgxmcxes .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#olsgxmcxes .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#olsgxmcxes .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#olsgxmcxes .gt_row_group_first td {
  border-top-width: 2px;
}

#olsgxmcxes .gt_row_group_first th {
  border-top-width: 2px;
}

#olsgxmcxes .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#olsgxmcxes .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#olsgxmcxes .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#olsgxmcxes .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#olsgxmcxes .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#olsgxmcxes .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#olsgxmcxes .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#olsgxmcxes .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#olsgxmcxes .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#olsgxmcxes .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#olsgxmcxes .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#olsgxmcxes .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#olsgxmcxes .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#olsgxmcxes .gt_left {
  text-align: left;
}

#olsgxmcxes .gt_center {
  text-align: center;
}

#olsgxmcxes .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#olsgxmcxes .gt_font_normal {
  font-weight: normal;
}

#olsgxmcxes .gt_font_bold {
  font-weight: bold;
}

#olsgxmcxes .gt_font_italic {
  font-style: italic;
}

#olsgxmcxes .gt_super {
  font-size: 65%;
}

#olsgxmcxes .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#olsgxmcxes .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#olsgxmcxes .gt_indent_1 {
  text-indent: 5px;
}

#olsgxmcxes .gt_indent_2 {
  text-indent: 10px;
}

#olsgxmcxes .gt_indent_3 {
  text-indent: 15px;
}

#olsgxmcxes .gt_indent_4 {
  text-indent: 20px;
}

#olsgxmcxes .gt_indent_5 {
  text-indent: 25px;
}

#olsgxmcxes .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#olsgxmcxes div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id=""></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Broome">Broome</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Exmouth_Gulf">Exmouth_Gulf</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Gulf_of_Carpentaria">Gulf_of_Carpentaria</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Heywood_Shoal">Heywood_Shoal</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="New_Caledonia">New_Caledonia</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Pilbara">Pilbara</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Scott_Reef">Scott_Reef</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Shark_Bay">Shark_Bay</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Townsville">Townsville</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="p1" class="gt_row gt_left">Ashmore</td>
<td headers="Broome" class="gt_row gt_right">0.0307</td>
<td headers="Exmouth_Gulf" class="gt_row gt_right">0.0366</td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right">0.0973</td>
<td headers="Heywood_Shoal" class="gt_row gt_right">0.0255</td>
<td headers="New_Caledonia" class="gt_row gt_right">0.203</td>
<td headers="Pilbara" class="gt_row gt_right">0.0402</td>
<td headers="Scott_Reef" class="gt_row gt_right">0.0580</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.0554</td>
<td headers="Townsville" class="gt_row gt_right">0.0971</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Broome</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right">−0.000883</td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right">0.0777</td>
<td headers="Heywood_Shoal" class="gt_row gt_right">0.00339</td>
<td headers="New_Caledonia" class="gt_row gt_right">0.192</td>
<td headers="Pilbara" class="gt_row gt_right">0.00187</td>
<td headers="Scott_Reef" class="gt_row gt_right">0.00252</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.0104</td>
<td headers="Townsville" class="gt_row gt_right">0.0814</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Exmouth_Gulf</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right">0.0825</td>
<td headers="Heywood_Shoal" class="gt_row gt_right">0.00992</td>
<td headers="New_Caledonia" class="gt_row gt_right">0.211</td>
<td headers="Pilbara" class="gt_row gt_right">0.000713</td>
<td headers="Scott_Reef" class="gt_row gt_right">0.00902</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.0100</td>
<td headers="Townsville" class="gt_row gt_right">0.0867</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Gulf_of_Carpentaria</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right">0.0716</td>
<td headers="New_Caledonia" class="gt_row gt_right">0.168</td>
<td headers="Pilbara" class="gt_row gt_right">0.0824</td>
<td headers="Scott_Reef" class="gt_row gt_right">0.100</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.0991</td>
<td headers="Townsville" class="gt_row gt_right">0.0188</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Heywood_Shoal</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right"><br /></td>
<td headers="New_Caledonia" class="gt_row gt_right">0.267</td>
<td headers="Pilbara" class="gt_row gt_right">0.0102</td>
<td headers="Scott_Reef" class="gt_row gt_right">0.0436</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.0288</td>
<td headers="Townsville" class="gt_row gt_right">0.103</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">New_Caledonia</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right"><br /></td>
<td headers="New_Caledonia" class="gt_row gt_right"><br /></td>
<td headers="Pilbara" class="gt_row gt_right">0.200</td>
<td headers="Scott_Reef" class="gt_row gt_right">0.332</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.231</td>
<td headers="Townsville" class="gt_row gt_right">0.173</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Pilbara</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right"><br /></td>
<td headers="New_Caledonia" class="gt_row gt_right"><br /></td>
<td headers="Pilbara" class="gt_row gt_right"><br /></td>
<td headers="Scott_Reef" class="gt_row gt_right">0.0109</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.00833</td>
<td headers="Townsville" class="gt_row gt_right">0.0883</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Scott_Reef</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right"><br /></td>
<td headers="New_Caledonia" class="gt_row gt_right"><br /></td>
<td headers="Pilbara" class="gt_row gt_right"><br /></td>
<td headers="Scott_Reef" class="gt_row gt_right"><br /></td>
<td headers="Shark_Bay" class="gt_row gt_right">0.0264</td>
<td headers="Townsville" class="gt_row gt_right">0.151</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Shark_Bay</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right"><br /></td>
<td headers="New_Caledonia" class="gt_row gt_right"><br /></td>
<td headers="Pilbara" class="gt_row gt_right"><br /></td>
<td headers="Scott_Reef" class="gt_row gt_right"><br /></td>
<td headers="Shark_Bay" class="gt_row gt_right"><br /></td>
<td headers="Townsville" class="gt_row gt_right">0.109</td></tr>
  </tbody>
  
  
</table>
</div>

### *Hydrophis major*

``` r
fst_upper_matrix$HMA |>
    gt() |>
    fmt_missing(missing_text = "") |>
    cols_label(p1 = "") |>
    fmt_number(n_sigfig = 3)
```

<div id="ronmjhjqgr" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ronmjhjqgr table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#ronmjhjqgr thead, #ronmjhjqgr tbody, #ronmjhjqgr tfoot, #ronmjhjqgr tr, #ronmjhjqgr td, #ronmjhjqgr th {
  border-style: none;
}

#ronmjhjqgr p {
  margin: 0;
  padding: 0;
}

#ronmjhjqgr .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#ronmjhjqgr .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#ronmjhjqgr .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ronmjhjqgr .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ronmjhjqgr .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ronmjhjqgr .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ronmjhjqgr .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ronmjhjqgr .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#ronmjhjqgr .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#ronmjhjqgr .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ronmjhjqgr .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ronmjhjqgr .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#ronmjhjqgr .gt_spanner_row {
  border-bottom-style: hidden;
}

#ronmjhjqgr .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#ronmjhjqgr .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#ronmjhjqgr .gt_from_md > :first-child {
  margin-top: 0;
}

#ronmjhjqgr .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ronmjhjqgr .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#ronmjhjqgr .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#ronmjhjqgr .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#ronmjhjqgr .gt_row_group_first td {
  border-top-width: 2px;
}

#ronmjhjqgr .gt_row_group_first th {
  border-top-width: 2px;
}

#ronmjhjqgr .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ronmjhjqgr .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#ronmjhjqgr .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#ronmjhjqgr .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ronmjhjqgr .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ronmjhjqgr .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ronmjhjqgr .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#ronmjhjqgr .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ronmjhjqgr .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ronmjhjqgr .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ronmjhjqgr .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ronmjhjqgr .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ronmjhjqgr .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ronmjhjqgr .gt_left {
  text-align: left;
}

#ronmjhjqgr .gt_center {
  text-align: center;
}

#ronmjhjqgr .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ronmjhjqgr .gt_font_normal {
  font-weight: normal;
}

#ronmjhjqgr .gt_font_bold {
  font-weight: bold;
}

#ronmjhjqgr .gt_font_italic {
  font-style: italic;
}

#ronmjhjqgr .gt_super {
  font-size: 65%;
}

#ronmjhjqgr .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#ronmjhjqgr .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#ronmjhjqgr .gt_indent_1 {
  text-indent: 5px;
}

#ronmjhjqgr .gt_indent_2 {
  text-indent: 10px;
}

#ronmjhjqgr .gt_indent_3 {
  text-indent: 15px;
}

#ronmjhjqgr .gt_indent_4 {
  text-indent: 20px;
}

#ronmjhjqgr .gt_indent_5 {
  text-indent: 25px;
}

#ronmjhjqgr .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#ronmjhjqgr div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id=""></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Exmouth_Gulf">Exmouth_Gulf</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Fraser_Island">Fraser_Island</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Gulf_of_Carpentaria">Gulf_of_Carpentaria</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="New_Caledonia">New_Caledonia</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Pilbara">Pilbara</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Shark_Bay">Shark_Bay</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Weipa">Weipa</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="p1" class="gt_row gt_left">Broome</td>
<td headers="Exmouth_Gulf" class="gt_row gt_right">−0.00191</td>
<td headers="Fraser_Island" class="gt_row gt_right">0.00294</td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right">0.0238</td>
<td headers="New_Caledonia" class="gt_row gt_right">0.170</td>
<td headers="Pilbara" class="gt_row gt_right">0.00773</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.0126</td>
<td headers="Weipa" class="gt_row gt_right">0.0452</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Exmouth_Gulf</td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Fraser_Island" class="gt_row gt_right">−0.000300</td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right">0.0176</td>
<td headers="New_Caledonia" class="gt_row gt_right">0.174</td>
<td headers="Pilbara" class="gt_row gt_right">0.00586</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.0102</td>
<td headers="Weipa" class="gt_row gt_right">0.0367</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Fraser_Island</td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Fraser_Island" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right">0.0155</td>
<td headers="New_Caledonia" class="gt_row gt_right">0.206</td>
<td headers="Pilbara" class="gt_row gt_right">0.00681</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.00532</td>
<td headers="Weipa" class="gt_row gt_right">0.0481</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Gulf_of_Carpentaria</td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Fraser_Island" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="New_Caledonia" class="gt_row gt_right">0.146</td>
<td headers="Pilbara" class="gt_row gt_right">0.0269</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.0164</td>
<td headers="Weipa" class="gt_row gt_right">−0.00129</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">New_Caledonia</td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Fraser_Island" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="New_Caledonia" class="gt_row gt_right"><br /></td>
<td headers="Pilbara" class="gt_row gt_right">0.184</td>
<td headers="Shark_Bay" class="gt_row gt_right">0.122</td>
<td headers="Weipa" class="gt_row gt_right">0.251</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Pilbara</td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Fraser_Island" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="New_Caledonia" class="gt_row gt_right"><br /></td>
<td headers="Pilbara" class="gt_row gt_right"><br /></td>
<td headers="Shark_Bay" class="gt_row gt_right">0.00555</td>
<td headers="Weipa" class="gt_row gt_right">0.0203</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Shark_Bay</td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Fraser_Island" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="New_Caledonia" class="gt_row gt_right"><br /></td>
<td headers="Pilbara" class="gt_row gt_right"><br /></td>
<td headers="Shark_Bay" class="gt_row gt_right"><br /></td>
<td headers="Weipa" class="gt_row gt_right">−0.0188</td></tr>
  </tbody>
  
  
</table>
</div>

### *Hydrophis stokesii*

``` r
fst_upper_matrix$HST |>
    gt() |>
    fmt_missing(missing_text = "") |>
    cols_label(p1 = "") |>
    fmt_number(n_sigfig = 3)
```

<div id="ukqzffalgp" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ukqzffalgp table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#ukqzffalgp thead, #ukqzffalgp tbody, #ukqzffalgp tfoot, #ukqzffalgp tr, #ukqzffalgp td, #ukqzffalgp th {
  border-style: none;
}

#ukqzffalgp p {
  margin: 0;
  padding: 0;
}

#ukqzffalgp .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#ukqzffalgp .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#ukqzffalgp .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ukqzffalgp .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ukqzffalgp .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ukqzffalgp .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ukqzffalgp .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ukqzffalgp .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#ukqzffalgp .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#ukqzffalgp .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ukqzffalgp .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ukqzffalgp .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#ukqzffalgp .gt_spanner_row {
  border-bottom-style: hidden;
}

#ukqzffalgp .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#ukqzffalgp .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#ukqzffalgp .gt_from_md > :first-child {
  margin-top: 0;
}

#ukqzffalgp .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ukqzffalgp .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#ukqzffalgp .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#ukqzffalgp .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#ukqzffalgp .gt_row_group_first td {
  border-top-width: 2px;
}

#ukqzffalgp .gt_row_group_first th {
  border-top-width: 2px;
}

#ukqzffalgp .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ukqzffalgp .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#ukqzffalgp .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#ukqzffalgp .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ukqzffalgp .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ukqzffalgp .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ukqzffalgp .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#ukqzffalgp .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ukqzffalgp .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ukqzffalgp .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ukqzffalgp .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ukqzffalgp .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ukqzffalgp .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ukqzffalgp .gt_left {
  text-align: left;
}

#ukqzffalgp .gt_center {
  text-align: center;
}

#ukqzffalgp .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ukqzffalgp .gt_font_normal {
  font-weight: normal;
}

#ukqzffalgp .gt_font_bold {
  font-weight: bold;
}

#ukqzffalgp .gt_font_italic {
  font-style: italic;
}

#ukqzffalgp .gt_super {
  font-size: 65%;
}

#ukqzffalgp .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#ukqzffalgp .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#ukqzffalgp .gt_indent_1 {
  text-indent: 5px;
}

#ukqzffalgp .gt_indent_2 {
  text-indent: 10px;
}

#ukqzffalgp .gt_indent_3 {
  text-indent: 15px;
}

#ukqzffalgp .gt_indent_4 {
  text-indent: 20px;
}

#ukqzffalgp .gt_indent_5 {
  text-indent: 25px;
}

#ukqzffalgp .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#ukqzffalgp div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id=""></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Broome">Broome</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Exmouth_Gulf">Exmouth_Gulf</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Gulf_of_Carpentaria">Gulf_of_Carpentaria</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Heywood_Shoal">Heywood_Shoal</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="North_Kimberley">North_Kimberley</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="North_QLD">North_QLD</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Pilbara">Pilbara</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="p1" class="gt_row gt_left">Ashmore</td>
<td headers="Broome" class="gt_row gt_right">−0.0748</td>
<td headers="Exmouth_Gulf" class="gt_row gt_right">0.0324</td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right">0.102</td>
<td headers="Heywood_Shoal" class="gt_row gt_right">0.111</td>
<td headers="North_Kimberley" class="gt_row gt_right">0.0445</td>
<td headers="North_QLD" class="gt_row gt_right">0.00521</td>
<td headers="Pilbara" class="gt_row gt_right">0.0356</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Broome</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right">0.00219</td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right">0.00824</td>
<td headers="Heywood_Shoal" class="gt_row gt_right">0.00455</td>
<td headers="North_Kimberley" class="gt_row gt_right">−0.136</td>
<td headers="North_QLD" class="gt_row gt_right">−0.0493</td>
<td headers="Pilbara" class="gt_row gt_right">0.000309</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Exmouth_Gulf</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right">0.0845</td>
<td headers="Heywood_Shoal" class="gt_row gt_right">0.0888</td>
<td headers="North_Kimberley" class="gt_row gt_right">0.0525</td>
<td headers="North_QLD" class="gt_row gt_right">0.0573</td>
<td headers="Pilbara" class="gt_row gt_right">0.00373</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Gulf_of_Carpentaria</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right">0.145</td>
<td headers="North_Kimberley" class="gt_row gt_right">−0.00826</td>
<td headers="North_QLD" class="gt_row gt_right">−0.0200</td>
<td headers="Pilbara" class="gt_row gt_right">0.0961</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">Heywood_Shoal</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right"><br /></td>
<td headers="North_Kimberley" class="gt_row gt_right">0.105</td>
<td headers="North_QLD" class="gt_row gt_right">0.0925</td>
<td headers="Pilbara" class="gt_row gt_right">0.103</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">North_Kimberley</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right"><br /></td>
<td headers="North_Kimberley" class="gt_row gt_right"><br /></td>
<td headers="North_QLD" class="gt_row gt_right">−0.122</td>
<td headers="Pilbara" class="gt_row gt_right">0.0637</td></tr>
    <tr><td headers="p1" class="gt_row gt_left">North_QLD</td>
<td headers="Broome" class="gt_row gt_right"><br /></td>
<td headers="Exmouth_Gulf" class="gt_row gt_right"><br /></td>
<td headers="Gulf_of_Carpentaria" class="gt_row gt_right"><br /></td>
<td headers="Heywood_Shoal" class="gt_row gt_right"><br /></td>
<td headers="North_Kimberley" class="gt_row gt_right"><br /></td>
<td headers="North_QLD" class="gt_row gt_right"><br /></td>
<td headers="Pilbara" class="gt_row gt_right">0.0686</td></tr>
  </tbody>
  
  
</table>
</div>
