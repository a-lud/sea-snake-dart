#!/usr/bin/env python3

# Import libraries
import itertools

from pandas import read_csv
import ipyrad.analysis as ipa
import ipyparallel as ipp
import toyplot
import toyplot.png

# Connect to ipp to control threads
ipyclient = ipp.Client(cluster_id="ipyrad_ALA")

# Files/paths
data = "ALA-stringent.highQ.filtered.LD50k.snps.hdf5"
popmap = "ALA-popmap.txt"

# Burnin/nreps
burnin = 100000
numreps = 200000

# Population-map file
populations = read_csv(popmap, comment="#", sep=" ", names=["sample", "grouping"])

# Convert to dictionary
imap = dict(populations.groupby("grouping")["sample"].apply(list))

# 50% missing data per population group
minmap = {i: 0.5 for i in imap}

# Initiate the structure object
struct = ipa.structure(name="ALA", data=data, imap=imap, minmap=minmap, mincov=0.75)

# Burning and
struct.mainparams.burnin = burnin
struct.mainparams.numreps = numreps
struct.run(
    nreps=5,
    kpop=[2, 3, 4, 5, 6, 7, 8, 9],
    auto=True,
    show_cluster=True,
    ipyclient=ipyclient,
)

etable = struct.get_evanno_table([2, 3, 4, 5, 6, 7, 8, 9])

# Generate plot
canvas = toyplot.Canvas(width=1500, height=800)

# plot the mean log probability of the models in red
axes = canvas.cartesian(ylabel="estLnProbMean")
axes.plot(etable.estLnProbMean * -1, color="darkred", marker="o")
axes.y.spine.style = {"stroke": "darkred"}

# plot delta K with its own scale bar of left side and in blue
axes = axes.share(
    "x", ylabel="deltaK", ymax=etable.deltaK.max() + etable.deltaK.max() * 0.25
)
axes.plot(etable.deltaK, color="steelblue", marker="o")
axes.y.spine.style = {"stroke": "steelblue"}

# set x labels
axes.x.ticks.locator = toyplot.locator.Explicit(range(len(etable.index)), etable.index)
axes.x.label.text = "K (N ancestral populations)"

# Save figure
toyplot.png.render(canvas, f"/ALA-deltaK-burnin_{burnin}-numreps_{numreps}.png")

for i in range(2, 9):
    table = struct.get_clumpp_table(i)

    # sort list by columns
    table.sort_values(by=list(range(i)), inplace=True)

    # or, sort by a list of names (here taken from imap)
    onames = list(itertools.chain(*imap.values()))
    table = table.loc[onames]

    # build barplot
    canvas = toyplot.Canvas(width=1500, height=800)
    axes = canvas.cartesian(bounds=("10%", "90%", "10%", "45%"))
    axes.bars(table)

    # add labels to x-axis
    ticklabels = [i for i in table.index.tolist()]
    axes.x.ticks.locator = toyplot.locator.Explicit(labels=ticklabels)
    axes.x.ticks.labels.angle = -60
    axes.x.ticks.show = True
    axes.x.ticks.labels.offset = 10
    axes.x.ticks.labels.style = {"font-size": "12px"}

    toyplot.png.render(
        canvas, f"/ALA-structure-k{i}-burnin_{burnin}-numreps_{numreps}.png"
    )
