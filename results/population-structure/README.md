Population structure
================
2024-11-29

This directory contains results relating to diversity statistics, both
within and between populations. The script `population-structure.Rmd`
renders this `README` and is responsible for creating the relevant
outputs.

## Helper functions

The first helper function parses `STRUCTURE` output files into
dataframes which we can work with.

``` r
readStructure <- function(path) {
    k <- as.numeric(sub(".+-(.*)", "\\1", str_remove(basename(path), "\\..*") ))
    kpops <- paste0("k", 1:k)
    read_lines(path) |>
        str_remove("^\\s+") |>
        str_remove("\\s+$") |>
        str_remove(": ") |>
        str_replace_all("\\s+", "\t") |>
        as_tibble() |>
        separate(col = value, into = c("i1", "i2", "i3", "i4", kpops), sep = "\t") |>
        select(i1, starts_with("k")) |>
        pivot_longer(names_to = "k-value", values_to = "proportion", starts_with("k")) |>
        mutate(
            proportion = as.numeric(proportion),
            i1 = as.numeric(i1)
        )
}
```

The second helper function is simply for generating the PCAs. It handles
finding `x` and `y` maxima/minima and formatting.

``` r
# Function that takes the PCA list we'll generate below and plot the desired PCs
# for the specified population
plotPCA <- function(pca_list, species, pcA, pcB, xstep = 1, ystep = 1, pal = "Paired") {
    sp_pca <- pca_list |>
        pluck(species)
    
    # x-axis: max/min
    xmin <- sp_pca |>
        pull(pcA) |>
        min(na.rm = TRUE) |>
        floor()
    
    xmax <- sp_pca |>
        pull(pcA) |>
        max(na.rm = TRUE) |>
        ceiling()
    
    # y-axis: max/min
    ymin <- sp_pca |>
        pull(pcB) |>
        min(na.rm = TRUE) |>
        floor()
    
    ymax <- sp_pca |>
        pull(pcB) |>
        max(na.rm = TRUE) |>
        ceiling()
    
    sp_pca |>
        ggplot(aes(x = .data[[pcA]], y = .data[[pcB]], colour = population)) +
        geom_point(size = 2) +
        scale_color_brewer(palette = pal) +
        labs(title = glue::glue("{pcA} vs {pcB}")) +
        scale_y_continuous(
            limits = c(ymin, ymax),
            breaks = seq(ymin, ymax, by = ystep),
            labels = as.character(seq(ymin, ymax, by = ystep)),
            expand = c(0, 0)
        ) +
        scale_x_continuous(
            limits = c(xmin, xmax),
            breaks = seq(xmin, xmax, by = xstep),
            labels = as.character(seq(xmin, xmax, by = xstep)),
            expand = c(0, 0)
        ) +
        theme(
            legend.position = "bottom",
            legend.title = element_text(hjust = 0.5, face = "bold", size = 14),
            legend.title.position = "top",
            plot.title = element_text(hjust = 0.5),
            title = element_text(face = "bold")
        )
}
```

## PCA

Principal component analysis (PCA) was conducted using the [Ipyrad
analysis
toolkit](https://ipyrad.readthedocs.io/en/master/API-analysis/index.html).
The code responsible for the CSV files we’re loading in below can be
found
[here](https://github.com/a-lud/sea-snake-dart/blob/main/scripts/06-pca.ipynb).
Below, we’re simply creating custom figures.

Below we’re loading the CSV files as a list of dataframes and adapting
the column names slightly.

``` r
lst_pca <- fs::dir_ls(here("results", "population-structure", "pca"), glob = "*.csv") |>
    # For each vector element, we're setting the basename of the file as its name
    (\(x) set_names(x, str_remove(basename(x), "-.*")))() |>
    # Read the contents of each CSV and set the column names to PC... 
    imap(\(x, i) {
        tmp <- x |>
            read_csv(col_names = TRUE, col_types = cols())
        colnames(tmp) <- c("sample", paste0("PC", 1:(ncol(tmp) - 1)))
        pops |>
            filter(species == i) |>
            left_join(tmp)
    })
```

### *Aipysurus laevis*

Below is a patchworked plot of the first three PCs compared to each
other for *A. laevis*.

``` r
pc_1_2 <- plotPCA(pca_list = lst_pca, species = "ALA", pcA = "PC1", pcB = "PC2", xstep = 2)
pc_1_3 <- plotPCA(pca_list = lst_pca, species = "ALA", pcA = "PC1", pcB = "PC3", xstep = 2)
pc_2_3 <- plotPCA(pca_list = lst_pca, species = "ALA", pcA = "PC2", pcB = "PC3")


design <- "12\n34"

pc_1_2 + pc_1_3 + pc_2_3 + guide_area() +
    plot_layout(
        design = design,
        guides = "collect",
        axes = "collect"
    ) &
    guides(colour = guide_legend(ncol = 3, override.aes = list(size = 3)))
```

<img src="population-structure_files/figure-gfm/pca-ALA-1.png" width="100%" />

### *Hydrophis major*

Below is a patchworked plot of the first three PCs compared to each
other for *H. major*.

``` r
pc_1_2 <- plotPCA(pca_list = lst_pca, species = "HMA", pcA = "PC1", pcB = "PC2", xstep = 2, ystep = 2)
pc_1_3 <- plotPCA(pca_list = lst_pca, species = "HMA", pcA = "PC1", pcB = "PC3", xstep = 2, ystep = 2)
pc_2_3 <- plotPCA(pca_list = lst_pca, species = "HMA", pcA = "PC2", pcB = "PC3", ystep = 2)

pc_1_2 + pc_1_3 + pc_2_3 + guide_area() +
    plot_layout(
        design = design, 
        guides = "collect",
        axes = "collect"
    ) &
    guides(colour = guide_legend(ncol = 3, override.aes = list(size = 4)))
```

<img src="population-structure_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

### *Hydrophis stokesii*

Below is a patchworked plot of the first three PCs compared to each
other for *H. stokesii*.

``` r
pc_1_2 <- plotPCA(pca_list = lst_pca, species = "HST", pcA = "PC1", pcB = "PC2", xstep = 2, ystep = 2)
pc_1_3 <- plotPCA(pca_list = lst_pca, species = "HST", pcA = "PC1", pcB = "PC3", xstep = 2, ystep = 2)
pc_2_3 <- plotPCA(pca_list = lst_pca, species = "HST", pcA = "PC2", pcB = "PC3", xstep = 2, ystep = 2)

pc_1_2 + pc_1_3 + pc_2_3 + guide_area() +
    plot_layout(
        design = design, 
        guides = "collect",
        axes = "collect"
    ) &
    guides(colour = guide_legend(ncol = 3, override.aes = list(size = 4)))
```

<img src="population-structure_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

## STRUCTURE

Next, we generate the `STRUCTURE` plots for each species. We first read
in each output file as a dataframe, combining the results into a single
long-format tibble.

``` r
pops_structure <- pops |> mutate(i1 = row_number(), .by = species)

structure_df <- fs::dir_ls(
    here("results", "population-structure"), 
    glob = "*clumpp.outfile", 
    recurse = TRUE
) |>
    (\(x) set_names(x, str_remove(basename(x), "\\..*")))() |>
    map(readStructure) |>
    list_rbind(names_to = "species") |>
    separate_wider_delim(delim = "-", names = c("species", "K"), too_many = "merge", cols = species) |>
    left_join(pops_structure)
```

### Choosing the correct K (N-ancestral populations)

Choosing the best `K` value can be determined by plotting the estimated
log-probability for each `K` value, in addition to the delta-K value.
The values can be interpreted as follows:

- **Estimated log-probability**: Higher values equal better model fit.
- **Delta-K**: The rate of change between successive K-values. Larger
  values indicate reater change in the model fit.

The results for the three species are shown below.

``` r
etable <- fs::dir_ls(here("results"), glob = "*etable.csv", recurse = TRUE) |>
    read_csv(col_names = TRUE, col_types = cols(), id = "species") |>
    mutate(species = str_remove(basename(species), "-etable.csv")) |>
    rename(K = `...1`)

etable |>
    gt(groupname_col = "species") |>
    
    fmt_number(
        columns = c("K", "Nreps"),
        drop_trailing_zeros = TRUE
    ) |>
    fmt_number(
        columns = 3:7,
        n_sigfig = 4
    )
```

<div id="zprxuovbrt" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#zprxuovbrt table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#zprxuovbrt thead, #zprxuovbrt tbody, #zprxuovbrt tfoot, #zprxuovbrt tr, #zprxuovbrt td, #zprxuovbrt th {
  border-style: none;
}
&#10;#zprxuovbrt p {
  margin: 0;
  padding: 0;
}
&#10;#zprxuovbrt .gt_table {
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
&#10;#zprxuovbrt .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#zprxuovbrt .gt_title {
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
&#10;#zprxuovbrt .gt_subtitle {
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
&#10;#zprxuovbrt .gt_heading {
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
&#10;#zprxuovbrt .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#zprxuovbrt .gt_col_headings {
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
&#10;#zprxuovbrt .gt_col_heading {
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
&#10;#zprxuovbrt .gt_column_spanner_outer {
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
&#10;#zprxuovbrt .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#zprxuovbrt .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#zprxuovbrt .gt_column_spanner {
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
&#10;#zprxuovbrt .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#zprxuovbrt .gt_group_heading {
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
&#10;#zprxuovbrt .gt_empty_group_heading {
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
&#10;#zprxuovbrt .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#zprxuovbrt .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#zprxuovbrt .gt_row {
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
&#10;#zprxuovbrt .gt_stub {
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
&#10;#zprxuovbrt .gt_stub_row_group {
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
&#10;#zprxuovbrt .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#zprxuovbrt .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#zprxuovbrt .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#zprxuovbrt .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#zprxuovbrt .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#zprxuovbrt .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#zprxuovbrt .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#zprxuovbrt .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#zprxuovbrt .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#zprxuovbrt .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#zprxuovbrt .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#zprxuovbrt .gt_footnotes {
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
&#10;#zprxuovbrt .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#zprxuovbrt .gt_sourcenotes {
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
&#10;#zprxuovbrt .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#zprxuovbrt .gt_left {
  text-align: left;
}
&#10;#zprxuovbrt .gt_center {
  text-align: center;
}
&#10;#zprxuovbrt .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#zprxuovbrt .gt_font_normal {
  font-weight: normal;
}
&#10;#zprxuovbrt .gt_font_bold {
  font-weight: bold;
}
&#10;#zprxuovbrt .gt_font_italic {
  font-style: italic;
}
&#10;#zprxuovbrt .gt_super {
  font-size: 65%;
}
&#10;#zprxuovbrt .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#zprxuovbrt .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#zprxuovbrt .gt_indent_1 {
  text-indent: 5px;
}
&#10;#zprxuovbrt .gt_indent_2 {
  text-indent: 10px;
}
&#10;#zprxuovbrt .gt_indent_3 {
  text-indent: 15px;
}
&#10;#zprxuovbrt .gt_indent_4 {
  text-indent: 20px;
}
&#10;#zprxuovbrt .gt_indent_5 {
  text-indent: 25px;
}
&#10;#zprxuovbrt .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#zprxuovbrt div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="K">K</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Nreps">Nreps</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="lnPK">lnPK</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="lnPPK">lnPPK</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="deltaK">deltaK</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="estLnProbMean">estLnProbMean</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="estLnProbStdev">estLnProbStdev</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr class="gt_group_heading_row">
      <th colspan="7" class="gt_group_heading" scope="colgroup" id="ALA">ALA</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="ALA  K" class="gt_row gt_right">2</td>
<td headers="ALA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="ALA  lnPK" class="gt_row gt_right">0</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">0</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−102,900</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">2187.6857</td></tr>
    <tr><td headers="ALA  K" class="gt_row gt_right">3</td>
<td headers="ALA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="ALA  lnPK" class="gt_row gt_right">3,396</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">1,319</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.9791</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−99,550</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">1347.3262</td></tr>
    <tr><td headers="ALA  K" class="gt_row gt_right">4</td>
<td headers="ALA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="ALA  lnPK" class="gt_row gt_right">2,077</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">589.5</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.3606</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−97,470</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">1634.8270</td></tr>
    <tr><td headers="ALA  K" class="gt_row gt_right">5</td>
<td headers="ALA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="ALA  lnPK" class="gt_row gt_right">2,666</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">8,055</td>
<td headers="ALA  deltaK" class="gt_row gt_right">4.412</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−94,800</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">1825.9338</td></tr>
    <tr><td headers="ALA  K" class="gt_row gt_right">6</td>
<td headers="ALA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="ALA  lnPK" class="gt_row gt_right">−5,389</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">11,020</td>
<td headers="ALA  deltaK" class="gt_row gt_right">1.131</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−100,200</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">9743.5471</td></tr>
    <tr><td headers="ALA  K" class="gt_row gt_right">7</td>
<td headers="ALA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="ALA  lnPK" class="gt_row gt_right">5,631</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">8,450</td>
<td headers="ALA  deltaK" class="gt_row gt_right">4.447</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−94,560</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">1900.3825</td></tr>
    <tr><td headers="ALA  K" class="gt_row gt_right">8</td>
<td headers="ALA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="ALA  lnPK" class="gt_row gt_right">−2,819</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">2,072</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.5548</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−97,380</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">3734.2613</td></tr>
    <tr><td headers="ALA  K" class="gt_row gt_right">9</td>
<td headers="ALA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="ALA  lnPK" class="gt_row gt_right">−4,891</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">0</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−102,300</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">10567.6624</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="7" class="gt_group_heading" scope="colgroup" id="HMA">HMA</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="HMA  K" class="gt_row gt_right">2</td>
<td headers="HMA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HMA  lnPK" class="gt_row gt_right">0</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">0</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−53,630</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">618.1925</td></tr>
    <tr><td headers="HMA  K" class="gt_row gt_right">3</td>
<td headers="HMA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HMA  lnPK" class="gt_row gt_right">62.50</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">1,125</td>
<td headers="HMA  deltaK" class="gt_row gt_right">1.560</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−53,570</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">721.1716</td></tr>
    <tr><td headers="HMA  K" class="gt_row gt_right">4</td>
<td headers="HMA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HMA  lnPK" class="gt_row gt_right">1,188</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">1,639</td>
<td headers="HMA  deltaK" class="gt_row gt_right">2.614</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−52,380</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">626.9818</td></tr>
    <tr><td headers="HMA  K" class="gt_row gt_right">5</td>
<td headers="HMA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HMA  lnPK" class="gt_row gt_right">−451.4</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">1,949</td>
<td headers="HMA  deltaK" class="gt_row gt_right">5.538</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−52,830</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">351.9319</td></tr>
    <tr><td headers="HMA  K" class="gt_row gt_right">6</td>
<td headers="HMA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HMA  lnPK" class="gt_row gt_right">−2,400</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">2,959</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0.8013</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−55,230</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">3693.2316</td></tr>
    <tr><td headers="HMA  K" class="gt_row gt_right">7</td>
<td headers="HMA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HMA  lnPK" class="gt_row gt_right">558.9</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">5,665</td>
<td headers="HMA  deltaK" class="gt_row gt_right">2.614</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−54,670</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">2166.7466</td></tr>
    <tr><td headers="HMA  K" class="gt_row gt_right">8</td>
<td headers="HMA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HMA  lnPK" class="gt_row gt_right">−5,106</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">1,707</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0.2278</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−59,780</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">7493.2254</td></tr>
    <tr><td headers="HMA  K" class="gt_row gt_right">9</td>
<td headers="HMA  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HMA  lnPK" class="gt_row gt_right">−3,399</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">0</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−63,180</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">19178.9866</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="7" class="gt_group_heading" scope="colgroup" id="HST">HST</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="HST  K" class="gt_row gt_right">2</td>
<td headers="HST  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HST  lnPK" class="gt_row gt_right">0</td>
<td headers="HST  lnPPK" class="gt_row gt_right">0</td>
<td headers="HST  deltaK" class="gt_row gt_right">0</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−53,170</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">819.7494</td></tr>
    <tr><td headers="HST  K" class="gt_row gt_right">3</td>
<td headers="HST  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HST  lnPK" class="gt_row gt_right">4,169</td>
<td headers="HST  lnPPK" class="gt_row gt_right">1,036</td>
<td headers="HST  deltaK" class="gt_row gt_right">0.7794</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−49,000</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">1329.4763</td></tr>
    <tr><td headers="HST  K" class="gt_row gt_right">4</td>
<td headers="HST  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HST  lnPK" class="gt_row gt_right">3,132</td>
<td headers="HST  lnPPK" class="gt_row gt_right">10,350</td>
<td headers="HST  deltaK" class="gt_row gt_right">7.640</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−45,870</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">1354.5491</td></tr>
    <tr><td headers="HST  K" class="gt_row gt_right">5</td>
<td headers="HST  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HST  lnPK" class="gt_row gt_right">−7,217</td>
<td headers="HST  lnPPK" class="gt_row gt_right">9,282</td>
<td headers="HST  deltaK" class="gt_row gt_right">2.488</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−53,090</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">3730.5217</td></tr>
    <tr><td headers="HST  K" class="gt_row gt_right">6</td>
<td headers="HST  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HST  lnPK" class="gt_row gt_right">2,065</td>
<td headers="HST  lnPPK" class="gt_row gt_right">5,107</td>
<td headers="HST  deltaK" class="gt_row gt_right">1.738</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−51,020</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">2939.0177</td></tr>
    <tr><td headers="HST  K" class="gt_row gt_right">7</td>
<td headers="HST  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HST  lnPK" class="gt_row gt_right">−3,042</td>
<td headers="HST  lnPPK" class="gt_row gt_right">6,647</td>
<td headers="HST  deltaK" class="gt_row gt_right">1.128</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−54,070</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">5892.4193</td></tr>
    <tr><td headers="HST  K" class="gt_row gt_right">8</td>
<td headers="HST  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HST  lnPK" class="gt_row gt_right">3,605</td>
<td headers="HST  lnPPK" class="gt_row gt_right">4,159</td>
<td headers="HST  deltaK" class="gt_row gt_right">2.407</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−50,460</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">1727.7308</td></tr>
    <tr><td headers="HST  K" class="gt_row gt_right">9</td>
<td headers="HST  Nreps" class="gt_row gt_right">5.000</td>
<td headers="HST  lnPK" class="gt_row gt_right">−554.3</td>
<td headers="HST  lnPPK" class="gt_row gt_right">0</td>
<td headers="HST  deltaK" class="gt_row gt_right">0</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−51,010</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">3519.8323</td></tr>
  </tbody>
  &#10;  
</table>
</div>

Plotting the table above looks like the following.

``` r
prb <- etable |>
    ggplot(aes(x = K, y = estLnProbMean, colour = species)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(
        name = "Estimated mean log-probability",
        breaks = seq(0, -105000, -10000),
        limits = c(-105e3, -4e4),
        labels = scales::label_number(style_negative = "minus")
    ) +
    scale_x_continuous(breaks = seq(2, 9, 1))

dk <- etable |>
    ggplot(aes(x = K, y = deltaK, colour = species)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(2, 9, 1)) +
    scale_y_continuous(
        name = "Delta-K",
        breaks = seq(0, 8, 1)
    )

prb/dk +
    plot_layout(axes = "collect", guides = "collect")
```

<img src="population-structure_files/figure-gfm/unnamed-chunk-5-1.png" width="100%" />

### *Aipysurus laevis*

The two highest log-probability values are at `k = 5` and `K = 7`. The
Delta-K for both of these values is quite high, indicating that both
K-values improved the model fit substantially.

``` r
structure_df |>
    filter(species == "ALA", K == "K-5") |>
    ggplot(aes(x = sample, y = proportion, fill = `k-value`)) +
    geom_col(position = "stack") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL) + 
    facet_grid(
        cols = vars(population), scales = "free_x", space = "free_x"
        
    ) +
    guides(fill = guide_legend(title.position = "bottom")) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.98),
        strip.text = element_text(angle = 90, vjust = 0.5, hjust = 0.1),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5),
        panel.spacing.x = unit(0.2, "line")
    )
```

<img src="population-structure_files/figure-gfm/unnamed-chunk-6-1.png" width="100%" />

### *Hydrophis major*

`K = 5` appears to be the best choice for *H. major*. The estimated
log-probability is relatively stable for most values of `K`, before
dipping once reaching higher values. Delta-K peaks at `K = 5`, with
little improvement after.

``` r
structure_df |>
    filter(species == "HMA", K == "K-5") |>
    ggplot(aes(x = sample, y = proportion, fill = `k-value`)) +
    geom_col(position = "stack") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL) + 
    facet_grid(
        cols = vars(population), scales = "free_x", space = "free_x"
        
    ) +
    guides(fill = guide_legend(title.position = "bottom")) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.98),
        strip.text = element_text(angle = 90, vjust = 0.5, hjust = 0.1),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5),
        panel.spacing.x = unit(0.2, "line")
    )
```

<img src="population-structure_files/figure-gfm/unnamed-chunk-7-1.png" width="100%" />

### *Hydrophis stokesii*

The estimated log-probabilities for *stokesii* show a similar pattern to
*H. major*, being consistent across multiple `K` values. Delta-K spikes
at `K = 4`, indicating that this is likely the best value for the
current dataset.

``` r
structure_df |>
    filter(species == "HST", K == "K-4") |>
    ggplot(aes(x = sample, y = proportion, fill = `k-value`)) +
    geom_col(position = "stack") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL) +
    facet_grid(
        cols = vars(population), scales = "free_x", space = "free_x"
        
    ) +
    guides(fill = guide_legend(title.position = "bottom")) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.98),
        strip.text = element_text(angle = 90, vjust = 0.5, hjust = 0.1),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5),
        panel.spacing.x = unit(0.2, "line")
    )
```

<img src="population-structure_files/figure-gfm/unnamed-chunk-8-1.png" width="100%" />