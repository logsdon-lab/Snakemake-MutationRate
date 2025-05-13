import os
import json
import argparse
import matplotlib
import matplotlib.axes
import numpy as np
import polars as pl
import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap, rgb2hex

IN_COLS = ["ref", "ref_st", "ref_end", "qry", "qry_sm", "qry_sm_divergence_time", "mu"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--infile",
        help=f"BED file of mutation rates. Format: {IN_COLS}",
        required=True,
        type=str,
    )
    ap.add_argument(
        "-c",
        "--sample_colors",
        help="Sample colors JSON as {'sm': 'color'}",
        default=None,
        type=str,
    )
    ap.add_argument(
        "-s",
        "--sample_shapes",
        help="Sample shapes JSON as {'sm': 'shape'}",
        default=None,
        type=str,
    )
    ap.add_argument("-o", "--outplot", help="Output plot", default=None, type=str)
    args = ap.parse_args()

    if args.outplot:
        outplot = args.outplot
    else:
        fname, _ = os.path.splitext(os.path.basename(args.infile))
        outplot = f"{fname}.pdf"

    if args.sample_colors:
        sample_colors = json.loads(args.sample_colors)
    else:
        sample_colors = {}

    if args.sample_shapes:
        sample_shapes = json.loads(args.sample_shapes)
    else:
        sample_shapes = {}

    df = (
        pl.read_csv(args.infile, has_header=False, separator="\t", new_columns=IN_COLS)
        .with_columns(
            # TODO: Should be in mutation rate.
            hap=pl.col("qry").str.extract(
                r"(mat|pat|haplotype1|haplotype2|hap1|hap2|h1|h2)"
            )
        )
        .with_columns(
            hap=pl.when(
                (pl.col("hap") == "pat")
                | (pl.col("hap") == "haplotype1")
                | (pl.col("hap") == "hap1")
                | (pl.col("hap") == "h1")
            )
            .then(pl.lit("H1"))
            .when(
                (pl.col("hap") == "mat")
                | (pl.col("hap") == "haplotype2")
                | (pl.col("hap") == "hap2")
                | (pl.col("hap") == "h2")
            )
            .then(pl.lit("H2"))
            .otherwise(pl.col("hap").fill_null(pl.lit("H1")))
        )
        .with_columns(
            pos=((pl.col("ref_end") - pl.col("ref_st")) / 2 + pl.col("ref_st"))
        )
        .sort(by=["qry_sm"])
    )

    unique_vals = df["qry_sm"].unique(maintain_order=True)
    # Generate random number of colors.
    colors = np.random.rand(len(unique_vals), 3)
    cmap = ListedColormap(colors)
    sm_color_mapping = {
        val: rgb2hex(color)
        for val, color in zip(unique_vals, cmap(np.linspace(0, 1, len(unique_vals))))
    }
    fig, axes = plt.subplots(figsize=(16, 8), nrows=2, height_ratios=[0.8, 0.2])
    ax: matplotlib.axes.Axes = axes[0]
    legend_ax: matplotlib.axes.Axes = axes[1]

    for sm, df_sm in df.group_by(["qry_sm"], maintain_order=True):
        sm = sm[0]
        ax.plot(
            df_sm["pos"],
            df_sm["mu"],
            sample_shapes.get(sm, "o"),
            markersize=5,
            markeredgecolor="black",
            markeredgewidth=0.1,
            color=sample_colors.get(sm, sm_color_mapping[sm]),
            label=sm,
        )

    ax.set_ylim(10e-10, 10e-4)
    ax.set_yscale("log")
    ax.set_ylabel("# of mutations per site per generation")
    ax.set_xlabel("Position (bp)")
    for side in ["right", "top"]:
        ax.spines[side].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()

    # Minimize legend axis.
    legend_ax.set_xticks([], [])
    legend_ax.set_yticks([], [])
    for side in ["left", "right", "top", "bottom"]:
        legend_ax.spines[side].set_visible(False)

    legend_ax.legend(
        frameon=False, loc="center", handles=handles, labels=labels, ncols=10
    )
    fig.savefig(outplot, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
