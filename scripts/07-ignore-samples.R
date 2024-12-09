# Create "drop-samples" list for each species: STRUCTURE
#
# Only want to run structure on certain populations. This script creates a
# "drop" list.

library(tidyverse)
library(here)

pops <- fs::dir_ls(here("data", "popmaps"), glob = "*.txt") |>
    read_delim(
        delim = " ",
        col_names = c("sample", "population"),
        col_types = cols(),
        id = "species"
    ) |>
    mutate(species = str_remove(basename(species), "-.*"))


# A. laevis:
pops |>
    filter(species == "ALA", ! population %in% c("Shark_Bay", "Pilbara", "Broome", "Exmouth_Gulf") ) |>
    select(population) |>
    distinct() |>
    pull(population) |>
    write_lines(here("results", "population-structure", "ALA-ignore.txt"))

# H. major
pops |>
    filter(species == "HMA", ! population %in% c("Shark_Bay", "Pilbara", "Broome", "Exmouth_Gulf") ) |>
    select(population) |>
    distinct() |>
    pull(population) |>
    write_lines(here("results", "population-structure", "HMA-ignore.txt"))

# H. stokesii
pops |>
    filter(species == "HST", ! population %in% c("North_Kimberley", "Pilbara", "Broome", "Exmouth_Gulf") ) |>
    select(population) |>
    distinct() |>
    pull(population) |>
    write_lines(here("results", "population-structure", "HST-ignore.txt"))
