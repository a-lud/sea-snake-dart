Effective population size
================
Alastair Ludington
2025-04-14

- [Introduction](#introduction)
- [Data processing](#data-processing)
- [Population sizes through time: N<sub>e</sub>
  curves](#population-sizes-through-time-ne-curves)

## Introduction

MSMC2 was used to model effective population size changes through time
from whole genome resequencing data.

## Data processing

*Hydrophis major* and *stokesii* samples were aligned to the *H. major*
reference genome, while *A. laevis* samples were mapped to the *A.
laevis* reference genome. As the *A. laevis* genome is not chromosome
scale, the sex chromosomes are unidentified. Contigs from the *A.
laevis* genome were first mapped to the *H. major* reference, and any
contig that had $\geq$ 50% of its length map to the Z-chromosome in *H.
major* was filtered out.

After filtering the A. laevis reference genome for Z-chromosome contigs,
mappability files were generated for each reference, before resequencing
data were mapped to their respective genomes using BWA-mem2. Input files
for `MSMC2` were built using the accompanying scripts, before `MSMC2`
was run on each sample individually, essentially running PSMC (often
referred to as `PSMC'`).

## Population sizes through time: N<sub>e</sub> curves

Below we plot effective population sizes through time. The data are
scaled using a mutation rate of
$2 \times 10^{-9}$/per-site/per-generation and a generation time of 10
years.

``` r
# ---------------------------------------------------------------------------- #
# Scaling parameters
# mu <- 2e-9
mu <- 4.71e-9 # https://academic.oup.com/mbe/article/37/6/1744/5741420#204169216
gen <- 10

# ---------------------------------------------------------------------------- #
# MSMC2 data
df_msmc <- fs::dir_ls(
    path = here("results", "demography", "psmc"),
    glob = '*.final.txt',
    recurse = TRUE
) %>%
    magrittr::extract(!str_detect(., "bootstrap")) |>
    read_tsv(col_names = TRUE, col_types = cols(), id = 'tmp') |>
    separate(col = tmp, into = c("tmp", 'clock'), sep = "__") |>
    mutate(
        tmp = basename(tmp),
        clock = str_remove(clock, '.final.txt'),
        clock = str_replace_all(clock, '_', "*"),
        clock = str_replace_all(clock, "-", '+'),
        # Scaling the data
        Years = (left_time_boundary/mu)*gen,
        Ne = (1/lambda)/(2*mu)
    ) |>
    separate(col = tmp, into = c("Species", 'sample'), sep = "-") |>
    arrange(Species) |>
    filter(Species %in% c("ALA", "HMA", "HST")) |>
    mutate(Species = forcats::fct_inorder(Species)) |>
    mutate(
        lt = as.numeric(factor(sample, levels = unique(sample))),
        .by = Species
    ) |>
    mutate(
        lt = forcats::fct_inorder(as.character(lt)),
        Species = case_when(
            str_starts(Species, "A") ~ "A. laevis",
            str_starts(Species, "HM") ~ "H. major",
            str_starts(Species, "HS") ~ "H. stokesii"
        )
    )

# Number of snakes in each species
df_msmc |>
    summarise(count = length(unique(sample)), .by = c(Species))
```

    ## # A tibble: 3 × 2
    ##   Species     count
    ##   <chr>       <int>
    ## 1 A. laevis       3
    ## 2 H. major        5
    ## 3 H. stokesii     1

``` r
# ---------------------------------------------------------------------------- #
# MSMC2 stairway plot
plot_stairway <- df_msmc |>
    # Removal of an obviously aberrant signal in distant past
    filter(clock == "2*2+25*1+2*3") |>
    ggplot(
        aes(
            x = Years, y = Ne,
            colour = Species,
            linetype = lt
        )
    ) +
    geom_step(linewidth = 1.4) +
    scale_x_log10(
        limits = c(3e4, 100e6),
        labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE)
    ) +
    annotation_logticks(sides = 'b', outside = TRUE, ) +
    coord_cartesian(clip = "off" ) +
    scale_y_continuous(
        limits = c(0, 600e3),
        breaks = seq(0, 600e3, 1e5),
        labels = scales::label_number(scale = 1e-3),
        expand = c(0, 0)
    ) +
    labs(
        y = bquote("Effective population size "~(N[e]~x~10^3)),
        x = bquote("Years ago (x10"^6*')')
    ) +
    scale_color_brewer(palette = "Set1") +
    guides(linetype = 'none') +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(vjust = -1),
        strip.text = element_text(size = 18),
        legend.position = "bottom",
        legend.text = element_text(size = 16, face = "italic"),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.7)
    )
ragg::agg_png(
    filename = here("results", "demography", "msmc.png"),
    width = 1200, height = 1200,
    units = "px",
    res = 150
)
plot_stairway
invisible(dev.off())

plot_stairway
```

<img src="demography_files/figure-gfm/msmc-data-1.png" width="100%" />
