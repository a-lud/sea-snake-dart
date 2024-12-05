Sequence QC
================
2024-12-05

- [Read counts before/after QC
  pipeline](#read-counts-beforeafter-qc-pipeline)

All raw data in this project was run through the same QC pipeline, the
steps of which are shown below:

1.  [BBduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/):
    Remove adapter content/barcodes.
2.  [Kraken2](https://github.com/DerrickWood/kraken2/tree/master):
    Remove contaminants (human + bacteria)
3.  [Fastp](https://github.com/OpenGene/fastp): Filter on average
    quality and minimum length
4.  [MultiQC](https://seqera.io/multiqc/): Aggregate QC results

Aggregated QC results can be found
[here](https://github.com/a-lud/sea-snake-dart/tree/main/results/qc/multiqc)
in the `MultiQC` report. For a comprehensive look at how samples fared
overall (after being processed through `Ipyrad`) see
[here](https://github.com/a-lud/sea-snake-dart/tree/main/results/ipyrad).

``` r
data <- fs::dir_ls(
    here("results", "qc", "stats"),
    glob = "*.tsv",
    recurse = TRUE
) |>
    read_tsv(id = "id", col_types = cols()) |>
    mutate(
        id = str_remove(basename(id), ".sequence.statistics.tsv"),
        id = forcats::fct_inorder(id),
        file = str_remove(file, "\\..*")
    ) |>
    left_join(meta, by = join_by(file == id_clean))
```

## Read counts before/after QC pipeline

Below is the code used to generate a simple comparison table showing the
amount of data before and after filtering. The table is too long to
include here, but you can view the output
[here](https://github.com/a-lud/sea-snake-dart/blob/main/results/qc/raw_trimmed_reads_bases.csv).

``` r
data |>
    unite(col = "Species", sep = " ", Genus, Species) |>
    mutate(across(c(id, Species), str_to_sentence)) |>
    select(id, Species, Region, file, num_seqs, sum_len) |>
    pivot_wider(names_from = id, values_from = c(num_seqs, sum_len)) |>
    arrange(Species, Region, -sum_len_Trimmed) |>
    write_csv(file = here("results", "qc", "raw_trimmed_reads_bases.csv"))
```
