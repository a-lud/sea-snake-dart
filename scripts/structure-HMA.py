#!/usr/bin/env python3

# Import libraries
import itertools

import ipyparallel as ipp
import ipyrad.analysis as ipa
import toyplot
import toyplot.png
from pandas import read_csv

# Connect to ipp to control threads
ipyclient = ipp.Client(cluster_id="ipyrad_HMA")

# Files/paths
probs = "dropped-samples.tsv"
data = "HMA-stringent.highQ.filtered.LD50k.snps.hdf5"
popmap = "HMA-popmap.txt"
ign = "HMA-ignore.txt"

# Problem samples
problems = read_csv(probs, sep="\t", comment="#", names=["sample", "population"])
problems = dict(problems.groupby("population")["sample"].apply(list))

with open(ign) as f:
    ignore = f.read().splitlines()


# Accessory function for dropping samples
def drop_samples(probs: dict, pop_dict: dict) -> dict:
    for k, v in probs.items():
        if k in pop_dict.keys():
            for sample in v:
                if sample in pop_dict[k]:
                    print(f"Dropping: {sample}")
                    pop_dict[k].remove(sample)
    return pop_dict


def drop_pop(drop_list: list, pop_dict: dict) -> dict:
    for p in drop_list:
        if p in pop_dict.keys():
            print(f"Dropping: {p}")
            pop_dict.pop(p)
    return pop_dict


# Burnin/nreps
burnin = 1000000
numreps = 100000

# Population-map file
populations = read_csv(popmap, comment="#", sep=" ", names=["sample", "grouping"])

# Convert to dictionary
imap = dict(populations.groupby("grouping")["sample"].apply(list))
imap = drop_samples(probs=problems, pop_dict=imap)
imap = drop_pop(drop_list=ignore, pop_dict=imap)

# 50% missing data per population group
minmap = {i: 0.8 for i in imap}

# Initiate the structure object
struct = ipa.structure(name="HMA", data=data, imap=imap, minmap=minmap, mincov=0.9)

# Burning and
struct.mainparams.burnin = burnin
struct.mainparams.numreps = numreps
struct.extraparams.printqhat = 1
struct.extraparams.printlikes = 1

print("Main parameters:")
print(struct.mainparams)

print("Extra parameters:")
print(struct.extraparams)

struct.run(
    nreps=10,
    kpop=[2, 3, 4, 5],
    seed=12345,
    auto=True,
    show_cluster=True,
    ipyclient=ipyclient,
    force=True
)

etable = struct.get_evanno_table([2, 3, 4, 5])
etable.to_csv(f"HMA-burnin_{burnin}-numreps_{numreps}.csv")

# Generate plot
canvas = toyplot.Canvas(width=1500, height=800)
canvas.style.update({"background-color": "white"})

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
toyplot.png.render(canvas, f"HMA-deltaK-burnin_{burnin}-numreps_{numreps}.png")

for i in range(2, 6):
    table = struct.get_clumpp_table(i)

    # sort list by columns
    table.sort_values(by=list(range(i)), inplace=True)

    # or, sort by a list of names (here taken from imap)
    onames = list(itertools.chain(*imap.values()))
    table = table.loc[onames]
    table.to_csv(f"HMA-structure-k{i}-burnin_{burnin}-numreps_{numreps}.plot.csv")

    # build barplot
    canvas = toyplot.Canvas(width=1500, height=800)
    canvas.style.update({"background-color": "white"})
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
        canvas, f"HMA-structure-k{i}-burnin_{burnin}-numreps_{numreps}.png"
    )
