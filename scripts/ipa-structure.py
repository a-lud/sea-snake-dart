#!/usr/bin/env python3

import time
import itertools
import argparse 
from pathlib import Path
from sys import exit, stderr
from ipyparallel import Client
from pandas import read_csv
from ipyrad.analysis import structure
from loguru import logger


def setup_logging(log_file: str, level: str = "INFO"):
    logger.remove()

    # Print to screen
    logger.add(
        stderr,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
        level=level,
    )

    # To file
    logger.add(
        log_file,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=level,
    )


def parse_arguments():
    """Parse command line arguments for the analysis script."""
    parser = argparse.ArgumentParser(
        description="Ipyrad analusis toolkit: STRUCTURE",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Mandatory arguments
    parser.add_argument("name", type=str, help="Name of the analysis")
    parser.add_argument("hdf5", type=str, help="Input HDF5 file")
    parser.add_argument("popmap", type=str, help="Population map file")

    parser.add_argument(
        "cluster_id",
        type=str,
        help="Name of the `ipyparallel` cluster to use (start before running this script)",
    )

    # Optional arguments
    parser.add_argument(
        "-p", "--problem_samples", type=str, help="File containing problem samples"
    )

    parser.add_argument(
        "-x",
        "--exclude-populations",
        type=str,
        help="File containing populations to exclude",
    )

    parser.add_argument(
        "-b",
        "--burnin",
        type=int,
        default=1000000,
        help="Burn-in value (default: %(default)s)",
    )

    parser.add_argument(
        "-n",
        "--numreps",
        type=int,
        default=100000,
        help="Number of repetitions (default: %(default)s)",
    )

    parser.add_argument(
        "-m",
        "--min-missing",
        type=float,
        default=0.75,
        help="SNPs must be present in at least this proportion of samples within a population (0-1, default: %(default)s)",
    )

    parser.add_argument(
        "-c",
        "--min-cov",
        type=float,
        default=0.75,
        help="SNPs must be present in at least this proportion of all samples (0-1, default: %(default)s)",
    )

    parser.add_argument(
        "-r",
        "--nreps",
        type=int,
        default=10,
        help="Number of STRUCTURE repetitions (default: %(default)s)",
    )

    parser.add_argument(
        "-k",
        "--k_pops",
        type=int,
        default=(2, 5),
        nargs=2,
        metavar=("k_min", "k_max"),
        help="Range of K-population values to test (default: %(default)s)",
    )

    args = parser.parse_args()

    # Validate min_missing range
    if not 0 <= args.min_missing <= 1:
        parser.error("--min-missing must be between 0 and 1")

    if not 0 <= args.min_cov <= 1:
        parser.error("--min-cov must be between 0 and 1")

    # Validate file existence for mandatory files
    for file_arg in [args.hdf5, args.popmap]:
        if not Path(file_arg).exists():
            parser.error(f"File does not exist: {file_arg}")

    # Validate optional file if provided
    if args.problem_samples and not Path(args.problem_samples).exists():
        parser.error(f"File does not exist: {args.problem_samples}")

    if args.exclude_populations and not Path(args.exclude_populations).exists():
        parser.error(f"File does not exist: {args.exclude_populations}")

    return args


# Accessory function for dropping samples
def drop_samples(probs: dict, pop_dict: dict) -> dict:
    for k, v in probs.items():
        if k in pop_dict.keys():
            for sample in v:
                if sample in pop_dict[k]:
                    logger.info(f"Dropping: {sample}")
                    pop_dict[k].remove(sample)
    return pop_dict


def drop_pop(drop_list: list, pop_dict: dict) -> dict:
    for p in drop_list:
        if p in pop_dict.keys():
            logger.info(f"Dropping: {p}")
            pop_dict.pop(p)
    return pop_dict


def main():
    """Run structure"""
    try:
        args = parse_arguments()
        start_time = time.time()
        time_str = time.strftime("%Y%m%d-%H%M%S")
        log_file = f"{args.name}_STRUCTURE_{time_str}.log"
        setup_logging(log_file)

        logger.info(f"Analysis: started at{time.strftime('%Y-%m-%d %H:%M:%S')}\n")

        # Your analysis logic here
        logger.info("-" * 10)
        logger.info("Arguments:")
        logger.info(f"Analysis name: {args.name}")
        logger.info(f"hdf5 file: {args.hdf5}")
        logger.info(f"Cluster ID: {args.cluster_id}")
        logger.info(f"Problem samples file: {args.problem_samples}")
        logger.info(f"Popmap file: {args.popmap}")
        logger.info(f"Exclude populations file: {args.exclude_populations}")
        logger.info(f"Burn-in: {args.burnin}")
        logger.info(f"Number of repetitions: {args.numreps}")
        logger.info(f"Minimum missing threshold: {args.min_missing}")
        logger.info(f"Minimum coverage threshold: {args.min_cov}")
        logger.info(f"Number of STRUCTURE repetitions: {args.nreps}")
        logger.info(f"K-population range: {args.k_pops}")
        logger.info("-" * 10 + "\n")

        # Connect to ipp to control threads
        logger.info(f"Connecting to ipyparallel cluster: {args.cluster_id}")
        ipyclient = Client(cluster_id=args.cluster_id)

        # Load population map file
        logger.info(f"Loading population map file: {args.popmap}")
        populations = read_csv(
            args.popmap, comment="#", sep="\t", names=["sample", "grouping"]
        )
        imap = dict(populations.groupby("grouping")["sample"].apply(list))

        # Drop samples if provided
        if args.problem_samples:
            logger.info(
                f"Dropping problem samples from population map: {args.problem_samples}"
            )
            problems = read_csv(
                args.problem_samples,
                sep="\t",
                comment="#",
                names=["sample", "population"],
            )
            problems = dict(problems.groupby("population")["sample"].apply(list))
            imap = drop_samples(probs=problems, pop_dict=imap)

        # Drop populations if provided
        if args.exclude_populations:
            logger.info(
                f"Dropping populations from population map: {args.exclude_populations}"
            )
            with open(args.exclude_populations) as f:
                remove = f.read().splitlines()
            imap = drop_pop(drop_list=remove, pop_dict=imap)

        # Required proportion of samples present at each SNP
        minmap = {i: args.min_missing for i in imap}

        # Initialise STRUCTURE object
        logger.info(f"Initialising STRUCTURE object")
        struct = structure(
            name=args.name,
            data=args.hdf5,
            imap=imap,
            minmap=minmap,
            mincov=args.min_cov,
            workdir=f"workdir-{args.name}",
        )

        # Set some parameters
        struct.mainparams.burnin = args.burnin
        struct.mainparams.numreps = args.numreps
        struct.extraparams.printqhat = 1
        struct.extraparams.printlikes = 1

        k_pops = list(range(args.k_pops[0], args.k_pops[1] + 1))
        logger.info(f"Running STRUCTURE with K-populations: {k_pops}")

        # Run STRUCTURE
        struct.run(
            nreps=args.nreps,
            kpop=k_pops,
            seed=12345,
            auto=True,
            show_cluster=True,
            ipyclient=ipyclient,
            force=True,
        )

        logger.info(f"Exporting delta-k table")
        # export delta-k table
        etable = struct.get_evanno_table(k_pops)
        etable.to_csv(
            f"{args.name}-burnin_{args.burnin}-numreps_{args.numreps}-etable.csv"
        )

        # export data for structure plots (by K)
        for i in k_pops:
            logger.info(f"Exporting plotting data for K-population: {i}")
            table = struct.get_clumpp_table(i)

            # sort list by columns
            table.sort_values(by=list(range(i)), inplace=True)

            # or, sort by a list of names (here taken from imap)
            onames = list(itertools.chain(*imap.values()))
            table = table.loc[onames]
            table.to_csv(
                f"{args.name}-burnin_{args.burnin}-numreps_{args.numreps}-k{i}.plot.csv"
            )
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(f"Pipeline: completed in {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")

    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        exit(1)
    except Exception as e:
        print(f"Error: {e}")
        exit(1)


if __name__ == "__main__":
    main()
