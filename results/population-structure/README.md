Population structure
================
2025-07-18

- [PCA](#pca)
  - [All samples](#all-samples)
  - [WA Coast](#wa-coast)
- [STRUCTURE](#structure)
  - [Choosing the correct K (N-ancestral
    populations)](#choosing-the-correct-k-n-ancestral-populations)
  - [*Aipysurus laevis*](#aipysurus-laevis)
  - [*Hydrophis major*](#hydrophis-major)
  - [*Hydrophis stokesii*](#hydrophis-stokesii)
- [Isolation by Distance (IBD)](#isolation-by-distance-ibd)
  - [Mantel statistics and spatial
    auto-correlation](#mantel-statistics-and-spatial-auto-correlation)
  - [IBD visualisation](#ibd-visualisation)
  - [Geographic distance vs Genetic
    distance](#geographic-distance-vs-genetic-distance)
  - [Spatial auto-correlation](#spatial-auto-correlation)
  - [Without Shark Bay samples](#without-shark-bay-samples)

This directory contains population structure results. The script
`population-structure.Rmd` renders this `README` and is responsible for
creating the relevant outputs.

## PCA

Below we use `SNPRelate` to perform principal component analysis (PCA).
The VCF files have already been filtered for problem samples, on MAF and
LD. Below we simply make the PCA figures.

We’ll also process the WA *Hydrophis stokesii* results separately, as
these data are stored in a separate VCF file.

### All samples

<img src="population-structure_files/figure-gfm/pca-ALA-1.png" width="100%" />

### WA Coast

And the same for *just* the WA coastline samples

<img src="population-structure_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

We can see that the Shark Bay samples partition off in both *A. laevis*
and *H. major*.

## STRUCTURE

`STRUCTURE` analyses were run on the filtered VCF files for a range of
`K` values for both all populations ($2\leq\text{K}\leq10$), and WA
populations only ($2\leq\text{K}\leq5$).

### Choosing the correct K (N-ancestral populations)

Choosing the best `K` value can be determined by plotting the estimated
log-probability for each `K` value, in addition to the delta-K value.
The values can be interpreted as follows:

- **Estimated log-probability**: Higher values equal better model fit.
- **Delta-K**: The rate of change between successive K-values. Larger
  values indicate greater change in the model fit.

The results for the three species are shown below.

<div id="wgekwmfapj" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#wgekwmfapj table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#wgekwmfapj thead, #wgekwmfapj tbody, #wgekwmfapj tfoot, #wgekwmfapj tr, #wgekwmfapj td, #wgekwmfapj th {
  border-style: none;
}
&#10;#wgekwmfapj p {
  margin: 0;
  padding: 0;
}
&#10;#wgekwmfapj .gt_table {
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
&#10;#wgekwmfapj .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#wgekwmfapj .gt_title {
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
&#10;#wgekwmfapj .gt_subtitle {
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
&#10;#wgekwmfapj .gt_heading {
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
&#10;#wgekwmfapj .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#wgekwmfapj .gt_col_headings {
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
&#10;#wgekwmfapj .gt_col_heading {
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
&#10;#wgekwmfapj .gt_column_spanner_outer {
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
&#10;#wgekwmfapj .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#wgekwmfapj .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#wgekwmfapj .gt_column_spanner {
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
&#10;#wgekwmfapj .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#wgekwmfapj .gt_group_heading {
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
&#10;#wgekwmfapj .gt_empty_group_heading {
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
&#10;#wgekwmfapj .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#wgekwmfapj .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#wgekwmfapj .gt_row {
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
&#10;#wgekwmfapj .gt_stub {
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
&#10;#wgekwmfapj .gt_stub_row_group {
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
&#10;#wgekwmfapj .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#wgekwmfapj .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#wgekwmfapj .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#wgekwmfapj .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#wgekwmfapj .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#wgekwmfapj .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#wgekwmfapj .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#wgekwmfapj .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#wgekwmfapj .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#wgekwmfapj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#wgekwmfapj .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#wgekwmfapj .gt_footnotes {
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
&#10;#wgekwmfapj .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#wgekwmfapj .gt_sourcenotes {
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
&#10;#wgekwmfapj .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#wgekwmfapj .gt_left {
  text-align: left;
}
&#10;#wgekwmfapj .gt_center {
  text-align: center;
}
&#10;#wgekwmfapj .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#wgekwmfapj .gt_font_normal {
  font-weight: normal;
}
&#10;#wgekwmfapj .gt_font_bold {
  font-weight: bold;
}
&#10;#wgekwmfapj .gt_font_italic {
  font-style: italic;
}
&#10;#wgekwmfapj .gt_super {
  font-size: 65%;
}
&#10;#wgekwmfapj .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#wgekwmfapj .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#wgekwmfapj .gt_indent_1 {
  text-indent: 5px;
}
&#10;#wgekwmfapj .gt_indent_2 {
  text-indent: 10px;
}
&#10;#wgekwmfapj .gt_indent_3 {
  text-indent: 15px;
}
&#10;#wgekwmfapj .gt_indent_4 {
  text-indent: 20px;
}
&#10;#wgekwmfapj .gt_indent_5 {
  text-indent: 25px;
}
&#10;#wgekwmfapj .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#wgekwmfapj div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="condition">condition</th>
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
      <th colspan="8" class="gt_group_heading" scope="colgroup" id="ALA">ALA</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="ALA  condition" class="gt_row gt_left">all</td>
<td headers="ALA  K" class="gt_row gt_right">2</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">0</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">0</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−76,210</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">708.8</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">all</td>
<td headers="ALA  K" class="gt_row gt_right">3</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">2,221</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">1,374</td>
<td headers="ALA  deltaK" class="gt_row gt_right">1.657</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−73,990</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">829.2</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">all</td>
<td headers="ALA  K" class="gt_row gt_right">4</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">846.5</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">479.8</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.4466</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−73,140</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">1,074</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">all</td>
<td headers="ALA  K" class="gt_row gt_right">5</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">366.7</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">870.2</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.9576</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−72,780</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">908.7</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">all</td>
<td headers="ALA  K" class="gt_row gt_right">6</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">−503.4</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">989.3</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.8666</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−73,280</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">1,142</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">all</td>
<td headers="ALA  K" class="gt_row gt_right">7</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">485.9</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">630.0</td>
<td headers="ALA  deltaK" class="gt_row gt_right">1.127</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−72,790</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">558.9</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">all</td>
<td headers="ALA  K" class="gt_row gt_right">8</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">−144.1</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">504.6</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.8852</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−72,940</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">570.0</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">all</td>
<td headers="ALA  K" class="gt_row gt_right">9</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">−648.7</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">748.2</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.9090</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−73,590</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">823.1</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">all</td>
<td headers="ALA  K" class="gt_row gt_right">10</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">99.53</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">0</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−73,490</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">761.9</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">wa</td>
<td headers="ALA  K" class="gt_row gt_right">2</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">0</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">0</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−57,160</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">545.5</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">wa</td>
<td headers="ALA  K" class="gt_row gt_right">3</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">−127.0</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">520.6</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.7276</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−57,290</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">715.4</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">wa</td>
<td headers="ALA  K" class="gt_row gt_right">4</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">−647.6</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">112.6</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0.1194</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−57,930</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">943.1</td></tr>
    <tr><td headers="ALA  condition" class="gt_row gt_left">wa</td>
<td headers="ALA  K" class="gt_row gt_right">5</td>
<td headers="ALA  Nreps" class="gt_row gt_right">10</td>
<td headers="ALA  lnPK" class="gt_row gt_right">−534.9</td>
<td headers="ALA  lnPPK" class="gt_row gt_right">0</td>
<td headers="ALA  deltaK" class="gt_row gt_right">0</td>
<td headers="ALA  estLnProbMean" class="gt_row gt_right">−58,470</td>
<td headers="ALA  estLnProbStdev" class="gt_row gt_right">4,177</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="8" class="gt_group_heading" scope="colgroup" id="HMA">HMA</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="HMA  condition" class="gt_row gt_left">all</td>
<td headers="HMA  K" class="gt_row gt_right">2</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">0</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">0</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−28,130</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">364.0</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">all</td>
<td headers="HMA  K" class="gt_row gt_right">3</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">653.3</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">389.3</td>
<td headers="HMA  deltaK" class="gt_row gt_right">1.836</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−27,470</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">212.0</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">all</td>
<td headers="HMA  K" class="gt_row gt_right">4</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">264.0</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">257.3</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0.9707</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−27,210</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">265.1</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">all</td>
<td headers="HMA  K" class="gt_row gt_right">5</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">6.720</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">462.7</td>
<td headers="HMA  deltaK" class="gt_row gt_right">1.093</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−27,200</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">423.2</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">all</td>
<td headers="HMA  K" class="gt_row gt_right">6</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">469.4</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">347.0</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0.8346</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−26,730</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">415.8</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">all</td>
<td headers="HMA  K" class="gt_row gt_right">7</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">122.4</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">98.57</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0.2876</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−26,610</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">342.7</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">all</td>
<td headers="HMA  K" class="gt_row gt_right">8</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">23.83</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">54.32</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0.1942</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−26,590</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">279.7</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">all</td>
<td headers="HMA  K" class="gt_row gt_right">9</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">78.15</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">96.72</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0.2982</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−26,510</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">324.3</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">all</td>
<td headers="HMA  K" class="gt_row gt_right">10</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">174.9</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">0</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−26,330</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">285.0</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">wa</td>
<td headers="HMA  K" class="gt_row gt_right">2</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">0</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">0</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−15,680</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">148.2</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">wa</td>
<td headers="HMA  K" class="gt_row gt_right">3</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">18.90</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">231.4</td>
<td headers="HMA  deltaK" class="gt_row gt_right">1.443</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−15,670</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">160.4</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">wa</td>
<td headers="HMA  K" class="gt_row gt_right">4</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">−212.5</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">137.1</td>
<td headers="HMA  deltaK" class="gt_row gt_right">1.095</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−15,880</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">125.2</td></tr>
    <tr><td headers="HMA  condition" class="gt_row gt_left">wa</td>
<td headers="HMA  K" class="gt_row gt_right">5</td>
<td headers="HMA  Nreps" class="gt_row gt_right">10</td>
<td headers="HMA  lnPK" class="gt_row gt_right">−75.40</td>
<td headers="HMA  lnPPK" class="gt_row gt_right">0</td>
<td headers="HMA  deltaK" class="gt_row gt_right">0</td>
<td headers="HMA  estLnProbMean" class="gt_row gt_right">−15,950</td>
<td headers="HMA  estLnProbStdev" class="gt_row gt_right">283.4</td></tr>
    <tr class="gt_group_heading_row">
      <th colspan="8" class="gt_group_heading" scope="colgroup" id="HST">HST</th>
    </tr>
    <tr class="gt_row_group_first"><td headers="HST  condition" class="gt_row gt_left">all</td>
<td headers="HST  K" class="gt_row gt_right">2</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">0</td>
<td headers="HST  lnPPK" class="gt_row gt_right">0</td>
<td headers="HST  deltaK" class="gt_row gt_right">0</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−47,760</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">460.1</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">all</td>
<td headers="HST  K" class="gt_row gt_right">3</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">230.6</td>
<td headers="HST  lnPPK" class="gt_row gt_right">195.7</td>
<td headers="HST  deltaK" class="gt_row gt_right">0.6032</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−47,530</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">324.4</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">all</td>
<td headers="HST  K" class="gt_row gt_right">4</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">426.3</td>
<td headers="HST  lnPPK" class="gt_row gt_right">407.8</td>
<td headers="HST  deltaK" class="gt_row gt_right">0.8004</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−47,100</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">509.5</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">all</td>
<td headers="HST  K" class="gt_row gt_right">5</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">18.42</td>
<td headers="HST  lnPPK" class="gt_row gt_right">303.3</td>
<td headers="HST  deltaK" class="gt_row gt_right">0.5053</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−47,080</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">600.3</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">all</td>
<td headers="HST  K" class="gt_row gt_right">6</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">−284.9</td>
<td headers="HST  lnPPK" class="gt_row gt_right">99.44</td>
<td headers="HST  deltaK" class="gt_row gt_right">0.2100</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−47,370</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">473.5</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">all</td>
<td headers="HST  K" class="gt_row gt_right">7</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">−185.4</td>
<td headers="HST  lnPPK" class="gt_row gt_right">63.81</td>
<td headers="HST  deltaK" class="gt_row gt_right">0.1139</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−47,550</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">560.2</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">all</td>
<td headers="HST  K" class="gt_row gt_right">8</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">−249.3</td>
<td headers="HST  lnPPK" class="gt_row gt_right">128.5</td>
<td headers="HST  deltaK" class="gt_row gt_right">0.2926</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−47,800</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">439.1</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">all</td>
<td headers="HST  K" class="gt_row gt_right">9</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">−377.8</td>
<td headers="HST  lnPPK" class="gt_row gt_right">77.05</td>
<td headers="HST  deltaK" class="gt_row gt_right">0.1165</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−48,180</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">661.5</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">all</td>
<td headers="HST  K" class="gt_row gt_right">10</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">−300.7</td>
<td headers="HST  lnPPK" class="gt_row gt_right">0</td>
<td headers="HST  deltaK" class="gt_row gt_right">0</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−48,480</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">523.1</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">wa</td>
<td headers="HST  K" class="gt_row gt_right">2</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">0</td>
<td headers="HST  lnPPK" class="gt_row gt_right">0</td>
<td headers="HST  deltaK" class="gt_row gt_right">0</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−64,780</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">697.7</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">wa</td>
<td headers="HST  K" class="gt_row gt_right">3</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">101.5</td>
<td headers="HST  lnPPK" class="gt_row gt_right">603.4</td>
<td headers="HST  deltaK" class="gt_row gt_right">1.696</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−64,680</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">355.9</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">wa</td>
<td headers="HST  K" class="gt_row gt_right">4</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">−501.9</td>
<td headers="HST  lnPPK" class="gt_row gt_right">801.5</td>
<td headers="HST  deltaK" class="gt_row gt_right">1.279</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−65,180</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">626.4</td></tr>
    <tr><td headers="HST  condition" class="gt_row gt_left">wa</td>
<td headers="HST  K" class="gt_row gt_right">5</td>
<td headers="HST  Nreps" class="gt_row gt_right">10</td>
<td headers="HST  lnPK" class="gt_row gt_right">299.6</td>
<td headers="HST  lnPPK" class="gt_row gt_right">0</td>
<td headers="HST  deltaK" class="gt_row gt_right">0</td>
<td headers="HST  estLnProbMean" class="gt_row gt_right">−64,890</td>
<td headers="HST  estLnProbStdev" class="gt_row gt_right">525.1</td></tr>
  </tbody>
  &#10;  
</table>
</div>

Plotting the table above looks like the following.

<img src="population-structure_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

### *Aipysurus laevis*

#### All populations

The `estLnMeanProb` and `Delta-K` values don’t quite coincide. We know
from the PCA there are five broad groups (which matches the
`estLnMeanProb` value), even though the `Delta-K` value suggests
$\text{K}=3$ is where there is the biggest rate of change. Granted,
other `K` values show some significant peaks and troughs.

Plotting $\text{K}=3$

<img src="population-structure_files/figure-gfm/ALA-STRUCTURE-ALL-K3-1.png" width="100%" />

And plotting $\text{K}=5$ as in the PCA.

<img src="population-structure_files/figure-gfm/ALA-STRUCTURE-ALL-1.png" width="100%" />

#### WA populations

From the PCA, we can see that Shark Bay samples appear to cluster
separately to the rest of the WA populations. The $\text{K} = 3$ model
shows the highest `Delta-K`, however, it really looks like there are two
populations - Shark Bay and the rest.

When we plot $\text{K} = 2$, there’s not really any discernable
difference bewteen Shark Bay and the remaining populations.

<img src="population-structure_files/figure-gfm/ALA-STRUCTURE-WA-K2-1.png" width="100%" />

When we plot $\text{K} = 3$, we can see that the Shark Bay samples have
a larger portion of $\text{K} = 1$ than the other populations.

<img src="population-structure_files/figure-gfm/ALA-STRUCTURE-WA-K3-1.png" width="100%" />

This unique population portion is even more evident when plotting
$\text{K} = 4$

<img src="population-structure_files/figure-gfm/ALA-STRUCTURE-WA-K4-1.png" width="100%" />

### *Hydrophis major*

#### All populations

The PCA for *H. major* indicates that there could be five broad
populations of snake. The `Delta-K` peak for *H. major* is $\text{K}=3$,
though there is a second peak at $\text{K}=5$ with a significant
drop-off after this.

<img src="population-structure_files/figure-gfm/HMA-STRUCTURE-ALL-K3-1.png" width="100%" />

Increasing to $\text{K}=5$ highlights the unique population structure of
Shark Bay, in addition to samples in the Gulf of Carpentaria and South
Queensland.

<img src="population-structure_files/figure-gfm/HMA-STRUCTURE-ALL-K5-1.png" width="100%" />

#### WA populations

The peak `estLnMeanProb` peaks at $\text{K}=3$, while `Delta-K` shows
peaks at $\text{K}=3$ and $\text{K}=4$, with a significant fall-off
after the latter.

The $\text{K}=3$ plot clearly shows Shark Bay as comprising unique
population structure.

<img src="population-structure_files/figure-gfm/HMA-STRUCTURE-WA-K3-1.png" width="100%" />

Increasing $\text{K}=4$ results in two Shark Bay samples having a
significant portion of a unique population `K=3`.

<img src="population-structure_files/figure-gfm/HMA-STRUCTURE-WA-K4-1.png" width="100%" />

### *Hydrophis stokesii*

#### All populations

STRUCTURE models that there is likely between 3-5 genetic groups in *H.
stokesii*.

<img src="population-structure_files/figure-gfm/HST-STRUCTURE-ALL-K3-1.png" width="100%" />

<img src="population-structure_files/figure-gfm/HST-STRUCTURE-ALL-K4-1.png" width="100%" />

<img src="population-structure_files/figure-gfm/HST-STRUCTURE-ALL-K5-1.png" width="100%" />

#### WA populations

Interestingly, `Delta-K` and `estLnMeanProb` values suggest there are
three population groups in *H. stokesii* even though the PCA shows a
homogeneous cluster.

Plotting $\text{K}=2$ - $\text{K}=4$ illustrates that STRUCTURE seems to
simply be overfitting and making population groups when there are none.

<img src="population-structure_files/figure-gfm/HST-STRUCTURE-WA-K2-1.png" width="100%" />

## Isolation by Distance (IBD)

Lastly, we’ll perform Isolation-by-distance analysis. We’ll first load
the VCF data and estimate the genetic distance matrix (Nei’s genetic
distance).

Note: I also found this repo helpful -
<https://github.com/jdalapicolla/IBD_models.R> In addition to this
repository -
<https://anon0433.github.io/strong-individual-signatures/mantel_tests_and_spatial_autocorrelation#figure_4:_mantel-based_correlogram>

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7
    ##   header_line: 8
    ##   variant count: 2721
    ##   column count: 212
    ## Meta line 7 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 2721
    ##   Character matrix gt cols: 212
    ##   skip: 0
    ##   nrows: 2721
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant: 2721
    ## All variants processed

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7
    ##   header_line: 8
    ##   variant count: 1208
    ##   column count: 96
    ## Meta line 7 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 1208
    ##   Character matrix gt cols: 96
    ##   skip: 0
    ##   nrows: 1208
    ##   row_num: 0
    ## Processed variant 1000Processed variant: 1208
    ## All variants processed

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7
    ##   header_line: 8
    ##   variant count: 2177
    ##   column count: 96
    ## Meta line 7 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 2177
    ##   Character matrix gt cols: 96
    ##   skip: 0
    ##   nrows: 2177
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant: 2177
    ## All variants processed

### Mantel statistics and spatial auto-correlation

Next we’ll perform the Mantel tests and generate the autocorrelation
plots. We’re using `vegan` to perform Mantel tests as this package has
the autocorrelation function as well. Consequently, we’re not using the
in-built IBD function in `snpR`.

    ##         lon       lat      depth
    ## 1  117.2978 -20.44900 -17.857500
    ## 2  114.2231 -22.05562 -17.145000
    ## 3  114.2157 -22.06160 -16.465002
    ## 4  114.2180 -22.06772 -17.145000
    ## 5  114.2825 -22.11952 -11.245001
    ## 6  114.2825 -22.11952 -11.245001
    ## 7  114.2825 -22.11952 -11.245001
    ## 8  122.0936 -18.06773 -18.152500
    ## 9  121.8893 -17.96775 -28.327499
    ## 10 121.6570 -17.91592 -54.752499
    ## 11 121.6570 -17.91592 -54.752499
    ## 12 121.9214 -18.00228 -16.887501
    ## 13 121.9228 -17.99760 -16.887501
    ## 14 121.9701 -17.98095 -18.547501
    ## 15 121.9811 -17.98227 -18.547501
    ## 16 121.9867 -17.98167 -18.547501
    ## 17 121.9961 -17.98308 -18.547501
    ## 18 122.0952 -17.13131 -17.690001
    ## 19 122.6818 -16.47206 -28.812498
    ## 20 122.9477 -16.39226  -9.130001
    ## 21 122.3800 -16.53376 -33.902496
    ## 22 122.3188 -16.58201 -39.982498
    ## 23 114.1398 -21.99717 -10.345001
    ## 24 114.1163 -22.13187  -7.004999
    ## 25 114.1353 -22.11317 -13.525002
    ## 26 114.2189 -21.79992 -21.489998
    ## 27 114.2040 -22.15383 -12.495001
    ## 28 114.2825 -22.11952 -11.245001
    ## 29 114.2825 -22.11952 -11.245001
    ## 30 117.1041 -20.31729 -31.952499
    ## 31 117.3506 -20.33974 -23.782499
    ## 32 117.3880 -20.33490 -24.312500
    ## 33 117.3277 -20.32256 -29.072500
    ## 34 122.1694 -17.96593 -17.784998
    ## 35 122.1601 -17.97471 -12.107500
    ## 36 122.1450 -17.98799 -12.107500
    ## 37 122.1097 -18.07230 -18.152500
    ## 38 122.1593 -17.96580 -12.547500
    ## 39 114.7119 -21.50888 -32.342499
    ## 40 114.6918 -21.55711 -14.262499
    ## 41 114.7047 -21.47702 -63.919998
    ## 42 114.7221 -21.52123 -17.102499
    ## 43 114.7347 -21.68917  -1.822499
    ## 44 116.1159 -20.74362 -12.262502
    ## 45 116.2020 -20.81831  -7.042501
    ## 46 116.2483 -20.78315  -6.002501
    ## 47 116.2464 -20.78377  -6.002501
    ## 48 116.1621 -20.86957  -8.262502
    ## 49 119.7474 -19.72300 -27.109999
    ## 50 119.8719 -19.75854 -18.792501
    ## 51 119.9460 -19.79095 -12.032501
    ## 52 119.9420 -19.78478 -12.032501
    ## 53 121.0250 -18.71287 -45.762501
    ## 54 121.0250 -18.71287 -45.762501
    ## 55 117.9240 -20.15788 -31.982500
    ## 56 116.4369 -20.20950 -48.937500
    ## 57 116.0066 -20.50535 -35.562500
    ## 58 116.0534 -20.57125 -33.572498
    ## 59 115.7385 -20.34710 -48.075001
    ## 60 116.1428 -20.56015 -33.592499
    ## 61 122.1289 -17.74055 -17.307501
    ## 62 122.0457 -17.64753 -23.870001
    ## 63 122.0447 -17.64817 -23.870001
    ## 64 121.2675 -17.66360 -95.232498
    ## 65 121.7777 -16.67605 -33.590000
    ## 66 122.0958 -17.03932 -21.890001
    ## 67 122.0403 -17.01186 -21.750002
    ## 68 119.6252 -19.77451 -18.969999
    ## 69 119.6137 -19.77452 -18.969999
    ## 70 119.7481 -19.37013 -56.750000
    ## 71 114.0237 -26.37763  -1.017500
    ## 72 113.9903 -26.27869  -2.297500
    ## 73 114.0207 -26.37396  -1.017500
    ## 74 114.1876 -26.37625  15.945000
    ## 75 114.0127 -26.36556  -1.997500
    ## 76 114.0136 -26.36605  -1.997500
    ## 77 114.0504 -26.39904  -1.017500
    ## 78 114.0397 -26.39275  -1.017500
    ## 79 114.0446 -26.39628  -1.017500
    ## 80 114.0244 -26.37882  -1.017500
    ## 81 114.2825 -22.11952 -11.245001

    ##         lon       lat       depth
    ## 1  122.1067 -17.97833 -13.9575005
    ## 2  122.0600 -18.00167 -19.2374992
    ## 3  122.0400 -18.00500 -19.2374992
    ## 4  122.0450 -17.98000 -19.2374992
    ## 5  122.0400 -18.00333 -19.2374992
    ## 6  122.0850 -17.99333 -13.9575005
    ## 7  122.0983 -17.99167 -13.9575005
    ## 8  122.0867 -17.99167 -13.9575005
    ## 9  122.0400 -18.03833 -16.6725006
    ## 10 114.2222 -22.05493 -17.1450005
    ## 11 114.2251 -22.05482 -17.1450005
    ## 12 114.2617 -22.10583 -13.9550009
    ## 13 114.2825 -22.11952 -11.2450008
    ## 14 114.2825 -22.11952 -11.2450008
    ## 15 114.2825 -22.11952 -11.2450008
    ## 16 114.2825 -22.11952 -11.2450008
    ## 17 114.2825 -22.11952 -11.2450008
    ## 18 122.2465 -18.08618 -10.4324999
    ## 19 122.2465 -18.08618 -10.4324999
    ## 20 122.2465 -18.08618 -10.4324999
    ## 21 122.9350 -15.73817 -67.4124985
    ## 22 122.2069 -16.97113 -17.4700012
    ## 23 122.1577 -16.98824 -18.8700008
    ## 24 116.8626 -20.30130 -36.1549988
    ## 25 117.0407 -20.30280 -32.6124992
    ## 26 117.2243 -20.31342 -29.1124992
    ## 27 117.3285 -20.32714 -29.0725002
    ## 28 117.3584 -20.31939 -29.1825008
    ## 29 113.5445 -25.32753  -4.2450013
    ## 30 113.6564 -25.29146  -6.6375003
    ## 31 113.4100 -25.24705 -19.4350014
    ## 32 113.3827 -25.37616 -11.1250010
    ## 33 113.4071 -25.38893 -11.1250010
    ## 34 113.3268 -25.18902 -18.7350006
    ## 35 113.7737 -25.61593  -7.0125003
    ## 36 113.7000 -25.80000   0.7300003
    ## 37 114.7242 -21.52025 -17.1024990
    ## 38 119.7372 -19.74326 -19.1999989
    ## 39 116.9746 -20.30999 -35.6250000
    ## 40 119.7377 -19.64632 -36.1199989
    ## 41 119.7984 -19.64341 -35.3525009
    ## 42 113.6483 -25.73222  -2.0425003
    ## 43 113.7000 -25.80000   0.7300003
    ## 44 113.6483 -25.73222  -2.0425003

<div id="ujqaelmtge" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
  &#10;  <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false" style="-webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji'; display: table; border-collapse: collapse; line-height: normal; margin-left: auto; margin-right: auto; color: #333333; font-size: 16px; font-weight: normal; font-style: normal; background-color: #FFFFFF; width: auto; border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8; border-right-style: none; border-right-width: 2px; border-right-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #A8A8A8; border-left-style: none; border-left-width: 2px; border-left-color: #D3D3D3;" bgcolor="#FFFFFF">
  <thead style="border-style: none;">
    <tr class="gt_col_headings" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3;">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Species" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">Species</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Observed-correlation" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">Observed correlation</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Significance" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">Significance</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Replicates" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">Replicates</th>
    </tr>
  </thead>
  <tbody class="gt_table_body" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3;">
    <tr style="border-style: none;"><td headers="Species" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">ALA</td>
<td headers="Observed correlation" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0980</td>
<td headers="Significance" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0138</td>
<td headers="Replicates" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">10,000</td></tr>
    <tr style="border-style: none;"><td headers="Species" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">HMA</td>
<td headers="Observed correlation" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.1534</td>
<td headers="Significance" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0003</td>
<td headers="Replicates" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">10,000</td></tr>
    <tr style="border-style: none;"><td headers="Species" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">HST</td>
<td headers="Observed correlation" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0824</td>
<td headers="Significance" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0818</td>
<td headers="Replicates" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">10,000</td></tr>
  </tbody>
  &#10;  
</table>
</div>

### IBD visualisation

### Geographic distance vs Genetic distance

Lastly, let’s generate a scatter plot of distance vs genetic distance.
Add the correlation coefficients and p-values in Inkscape or something.

<img src="population-structure_files/figure-gfm/IBD-correlation-plot-1.png" width="100%" />

### Spatial auto-correlation

<img src="population-structure_files/figure-gfm/IBD-spatial-auto-correlation-1.png" width="100%" />

### Without Shark Bay samples

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7
    ##   header_line: 8
    ##   variant count: 2721
    ##   column count: 212
    ## Meta line 7 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 2721
    ##   Character matrix gt cols: 212
    ##   skip: 0
    ##   nrows: 2721
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant: 2721
    ## All variants processed

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7
    ##   header_line: 8
    ##   variant count: 1208
    ##   column count: 96
    ## Meta line 7 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 1208
    ##   Character matrix gt cols: 96
    ##   skip: 0
    ##   nrows: 1208
    ##   row_num: 0
    ## Processed variant 1000Processed variant: 1208
    ## All variants processed

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 7
    ##   header_line: 8
    ##   variant count: 2177
    ##   column count: 96
    ## Meta line 7 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 2177
    ##   Character matrix gt cols: 96
    ##   skip: 0
    ##   nrows: 2177
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant: 2177
    ## All variants processed

<div id="dflswbyhki" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
  &#10;  <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false" style="-webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji'; display: table; border-collapse: collapse; line-height: normal; margin-left: auto; margin-right: auto; color: #333333; font-size: 16px; font-weight: normal; font-style: normal; background-color: #FFFFFF; width: auto; border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8; border-right-style: none; border-right-width: 2px; border-right-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #A8A8A8; border-left-style: none; border-left-width: 2px; border-left-color: #D3D3D3;" bgcolor="#FFFFFF">
  <thead style="border-style: none;">
    <tr class="gt_col_headings" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3;">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Species" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: left;" bgcolor="#FFFFFF" valign="bottom" align="left">Species</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Observed-correlation" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">Observed correlation</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Significance" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">Significance</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Replicates" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">Replicates</th>
    </tr>
  </thead>
  <tbody class="gt_table_body" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3;">
    <tr style="border-style: none;"><td headers="Species" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">ALA</td>
<td headers="Observed correlation" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0356</td>
<td headers="Significance" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0936</td>
<td headers="Replicates" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">10,000</td></tr>
    <tr style="border-style: none;"><td headers="Species" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">HMA</td>
<td headers="Observed correlation" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0546</td>
<td headers="Significance" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.1336</td>
<td headers="Replicates" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">10,000</td></tr>
    <tr style="border-style: none;"><td headers="Species" class="gt_row gt_left" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">HST</td>
<td headers="Observed correlation" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0824</td>
<td headers="Significance" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">0.0789</td>
<td headers="Replicates" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">10,000</td></tr>
  </tbody>
  &#10;  
</table>
</div>

<img src="population-structure_files/figure-gfm/unnamed-chunk-6-1.png" width="100%" />

<img src="population-structure_files/figure-gfm/unnamed-chunk-7-1.png" width="100%" />
