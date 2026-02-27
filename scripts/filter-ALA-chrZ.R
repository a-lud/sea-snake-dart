# ---------------------------------------------------------------------------- #
# Remove Z-chromosome sequences from A. laevis assembly for PSMC'
suppressPackageStartupMessages({
    library(tidyverse)
})

df <- read_tsv(
    'results/z-chromosome-filter/ALA-to-HMA.stats',
    col_types = cols()
)

fai <- read_tsv(
    'data/genomes/aipysurus_laevis-garvin.fa.fai',
    col_names = FALSE,
    col_types = cols()
) |>
    select(query_name = X1, q_len = X2)

df <- df |>
    left_join(fai)

# ---------------------------------------------------------------------------- #
# Identify Z-chromosome sequences in A. laevis
# query_name    q_len       cum_query_aln_width prop_aln
# contig_106    2373556     2348180.            0.989
# contig_60     7375929     6812746.            0.924
# contig_38     14481695    12991699.           0.897
# contig_121    1793549     1589810.            0.886
# contig_12     41696726    36503943.           0.875
# contig_44     12375658    10834198.           0.875
# contig_28     21725459    18511163.           0.852
# contig_80     3901514     3251567.            0.833
# contig_20     27652581    22492369.           0.813
# contig_129    1627560     1275529.            0.784
# contig_33     19045678    14358875.           0.754
# contig_66     5828272     3870666.            0.664
# contig_139    1386786     770716.             0.556
df.z <- df |>
    mutate(query_aln_width = abs(query_start - query_end)) |>
    select(
        ref_name = `#reference_name`, strand, query_name, query_aln_width,
        q_len, starts_with('per')
    ) |>
    filter(ref_name == 'chrZ') |>
    group_by(query_name, q_len) |>
    summarise(cum_query_aln_width = sum(query_aln_width)) |>
    arrange(q_len) |>
    mutate(prop_aln = cum_query_aln_width/q_len) |>
    filter(prop_aln >= 0.5) |>
    arrange(-prop_aln)

# Total len: 161,264,963bp
df.z |>
    pull(q_len) |>
    sum()

# Contigs to exclude from PSMC' analysis
fai |>
    filter(! query_name %in% df.z$query_name, q_len > 1e6) |>
    pull(query_name) |>
    write_lines('results/z-chromosome-filter/ALA-contigs-non-chrz.txt')

