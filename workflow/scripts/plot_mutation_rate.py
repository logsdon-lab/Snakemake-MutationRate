import os
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
    ap.add_argument("-o", "--outplot", help="Output plot", default=None, type=str)
    args = ap.parse_args()

    if args.outplot:
        outplot = args.outplot 
    else:
        fname, _ = os.path.splitext(os.path.basename(args.infile))
        outplot = f"{fname}.pdf"

    df = (
        pl.read_csv(
            args.infile,
            has_header=False,
            separator="\t",
            new_columns=IN_COLS
        )
        .with_columns(
            # TODO: Should be in mutation rate.
            hap=pl.col("qry").str.extract(r"(mat|pat|haplotype1|haplotype2|hap1|hap2)")
        )
        .with_columns(
            hap=pl.when((pl.col("hap") == "pat") | (pl.col("hap") == "haplotype1") | (pl.col("hap") == "hap1"))
            .then(pl.lit("H1"))
            .when((pl.col("hap") == "mat") | (pl.col("hap") == "haplotype2") | (pl.col("hap") == "hap2"))
            .then(pl.lit("H2"))
            .otherwise(pl.col("hap"))
        )
        .with_columns(
            pos=((pl.col("ref_end") - pl.col("ref_st")) / 2 + pl.col("ref_st")) 
        )
        .sort(by=["qry_sm"])
    )

    unique_vals = df.select(["qry_sm", "hap"]).unique(maintain_order=True)
    # Generate random number of colors.
    colors = np.random.rand(len(unique_vals), 3)
    cmap = ListedColormap(colors)
    val_color_mapping = {
        val: rgb2hex(color)
        for val, color in zip(
            unique_vals.iter_rows(), cmap(np.linspace(0, 1, len(unique_vals)))
        )
    }
    fig, ax = plt.subplots(figsize=(16, 8))
    ax: matplotlib.axes.Axes
    
    for sm_hap, df_sm in df.group_by(["qry_sm", "hap"], maintain_order=True):
        ax.plot(
            df_sm["pos"],
            df_sm["mu"],
            "o",
            markersize=5,
            markeredgecolor="black",
            markeredgewidth=0.1,
            color=val_color_mapping[sm_hap],
            label="_".join(sm_hap)
        )

    ax.set_ylabel("# of mutations per site per generation")
    ax.set_xlabel("Position (bp)")
    
    fig.legend(loc="center right")
    fig.savefig(outplot, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())