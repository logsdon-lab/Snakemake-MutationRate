import os
import sys
import ast
import json
import argparse
import statistics
import matplotlib
import matplotlib.axes
import numpy as np
import polars as pl
import intervaltree as it
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap, rgb2hex

IN_COLS = ["ref", "ref_st", "ref_end", "qry", "qry_sm", "qry_sm_divergence_time", "mu"]
OUT_COLS = ["mu_ref", "mu_other", "rel_rate"]


def minimize_ax(ax: matplotlib.axes.Axes):
    ax.set_xticks([], [])
    ax.set_yticks([], [])
    for side in ["left", "right", "top", "bottom"]:
        ax.spines[side].set_visible(False)


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
        "--outfile",
        help=f"BED file of mutation rates. Format: {OUT_COLS}",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "-r",
        "--ref_annot",
        help="Reference annotations as BED file. Does not check mutation rate chrom names. Must limit to range of infile.",
        default=None,
        type=argparse.FileType("rb"),
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
    ap.add_argument(
        "--ref_cmp_bed",
        help="BED file to calculate relative mutation rate against other regions.",
        default=None,
        type=str,
    )
    ap.add_argument("-o", "--outplot", help="Output plot", default=None, type=str)
    ap.add_argument(
        "--legend_ncols",
        help="Number of columns for legend elements.",
        default=2,
        type=str,
    )
    args = ap.parse_args()

    if args.outplot:
        outplot = args.outplot
    else:
        fname, _ = os.path.splitext(os.path.basename(args.infile))
        outplot = f"{fname}.png"

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
    if args.ref_annot:
        fig, axes = plt.subplots(
            figsize=(16, 8),
            nrows=4,
            height_ratios=[0.1, 0.8, 0.2, 0.01],
            layout="tight",
        )

        ax: matplotlib.axes.Axes = axes[1]
        legend_ax: matplotlib.axes.Axes = axes[2]

        ref_ax: matplotlib.axes.Axes = axes[0]
        minimize_ax(ref_ax)

        legend_ref_ax: matplotlib.axes.Axes = axes[3]
        minimize_ax(legend_ref_ax)
    else:
        fig, axes = plt.subplots(
            figsize=(16, 8), nrows=2, height_ratios=[0.8, 0.2], layout="tight"
        )
        ax: matplotlib.axes.Axes = axes[0]
        legend_ax: matplotlib.axes.Axes = axes[1]

    if args.ref_cmp_bed:
        df_ref_mut = pl.read_csv(
            args.ref_cmp_bed,
            separator="\t",
            has_header=False,
            columns=[0, 1, 2],
            new_columns=["chrom", "st", "end"],
        )
        itree_ref_mut = it.IntervalTree(
            it.Interval(itv["st"], itv["end"])
            for itv in df_ref_mut.iter_rows(named=True)
        )
    else:
        itree_ref_mut = it.IntervalTree()

    xmin, xmax = df["pos"].min(), df["pos"].max()
    ymin, ymax = 10e-11, 10e-5

    ax.set_ylim(ymin, ymax)
    ax.set_yscale("log")
    ax.set_ylabel("# of mutations per site per generation")

    for side in ["right", "top"]:
        ax.spines[side].set_visible(False)

    ref_muts = []
    other_muts = []
    for sm, df_sm in df.group_by(["qry_sm"], maintain_order=True):
        sm = sm[0]
        for st, end, mu in df_sm.select("ref_st", "ref_end", "mu").iter_rows():
            if itree_ref_mut.overlaps(st, end):
                ref_muts.append(mu)
            else:
                other_muts.append(mu)

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

    if ref_muts:
        median_mut = statistics.median(ref_muts)
        ax.axhline(median_mut, linestyle="dotted", color="red")
    else:
        median_mut = 0.0

    if other_muts:
        other_mut = statistics.median(other_muts)
        if other_mut == 0:
            other_mut = statistics.mean(other_muts)
        ax.axhline(other_mut, linestyle="dotted", color="black")
    else:
        other_mut = 1.0

    rel_mut_rate = round(median_mut / other_mut)
    print("\t".join(OUT_COLS), file=args.outfile)
    print(f"{median_mut}\t{other_mut}\t{rel_mut_rate}", file=args.outfile)

    ax.xaxis.set_major_formatter("plain")
    new_xtick_labels = []
    xticks = [xmin, *ax.get_xticks(), xmax]
    for tk in xticks:
        # Convert units and round.
        new_x_txt = round(tk / 1_000_000, 1)
        new_xtick_labels.append(str(new_x_txt))

    ax.set_xticks(xticks, new_xtick_labels)
    ax.set_xlabel("Position (Mbp)")
    ax.set_xlim(xmin, xmax)

    # Minimize legend axis.
    minimize_ax(legend_ax)

    if args.ref_annot and ref_ax:
        ref_ax.set_xlim(xmin, xmax)

        df_ref_annot = pl.read_csv(
            args.ref_annot,
            separator="\t",
            has_header=False,
            new_columns=[
                "chrom",
                "st",
                "end",
                "name",
                "score",
                "strand",
                "thick_st",
                "thick_end",
                "item_rgb",
            ],
        )

        df_ref_annot = df_ref_annot.filter(
            (pl.col("st") >= xmin) & (pl.col("end") <= xmax)
        )
        ylim = ref_ax.get_ylim()
        ylim_ht = ylim[1] - ylim[0]
        for row in df_ref_annot.iter_rows(named=True):
            st, end, item_rgb, lbl = row["st"], row["end"], row["item_rgb"], row["name"]
            try:
                item_rgb = [elem / 255 for elem in ast.literal_eval(f"({item_rgb})")]
            except Exception:
                pass

            ref_ax.add_patch(
                Rectangle((st, 0), end - st, ylim_ht, color=item_rgb, label=lbl)
            )

        # Add border
        rect = Rectangle(
            (xmin, 0),
            xmax - xmin,
            ylim_ht,
            edgecolor="black",
            linewidth=1.0,
            fill=None,
        )
        ref_ax.add_patch(rect)

        handles, labels = ref_ax.get_legend_handles_labels()
        legend_elems: dict[str, Rectangle] = dict(zip(labels, handles))
        legend_ref_ax.legend(
            frameon=False,
            loc="center",
            labels=legend_elems.keys(),
            handles=legend_elems.values(),
            ncols=args.legend_ncols,
        )

    handles, labels = ax.get_legend_handles_labels()
    legend_ax.legend(
        frameon=False, loc="center", handles=handles, labels=labels, ncols=10
    )
    fig.savefig(outplot, bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
