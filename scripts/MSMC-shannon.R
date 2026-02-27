# ---------------------------------------------------------------------------- #
# MSMC: H. major, H. stokesii & Aipysurus laevis
suppressPackageStartupMessages({
    library(tidyverse)
    library(here)
    library(patchwork)
})

# ---------------------------------------------------------------------------- #
# Scaling parameters
mu <- 2e-9
gen <- 10

# ---------------------------------------------------------------------------- #
# MSMC2 data
df <- fs::dir_ls(
    path = 'results/psmc',
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
df |>
    summarise(count = length(unique(sample)), .by = c(Species))

# ---------------------------------------------------------------------------- #
# MSMC2 stairway plot
ragg::agg_png(
    filename = "figures/msmc-shannon.png",
    width = 1000, height = 1000,
    units = "px",
    res = 100
)
df |>
    filter(Ne <= 3e6, clock == "2*2+25*1+2*3") |>
    ggplot(
        aes(
            x = Years, y = Ne,
            colour = Species,
            linetype = lt
        )
    ) +
    geom_step(linewidth = 1.2) +
    scale_x_log10(labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE)) +
    annotation_logticks(sides = 'b', outside = TRUE, ) +
    coord_cartesian(clip = "off" ) +
    scale_y_continuous(labels = scales::label_number(scale = 1e-6)) +
    labs(
        y = bquote("Effective population size "~(N[e]~x~10^6)),
        x = bquote("Years ago (x10"^6*')')
    ) +
    scale_color_brewer(palette = "Set2") +
    guides(linetype = 'none') +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(vjust = 0.1),
        strip.text = element_text(size = 18),
        legend.position = "bottom",
        legend.text = element_text(size = 16, face = "italic"),
        legend.title = element_blank()
    )
invisible(dev.off())
