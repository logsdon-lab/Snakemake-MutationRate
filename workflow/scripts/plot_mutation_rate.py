import os
import sys
import json
import argparse
import matplotlib
import matplotlib.axes
import polars as pl

from cenplot import (
    LegendTrackSettings,
    Track,
    TrackPosition,
    TrackType,
    LineTrackSettings,
    PlotSettings,
    read_tracks,
    plot_tracks,
)


OUT_COLS = ["mu_ref", "mu_other", "rel_rate"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--infile",
        help="BED file of mutation rates. Has header",
        required=True,
        type=str,
    )
    ap.add_argument(
        "--outfile",
        help=f"BED file of summary mutation rates. Format: {OUT_COLS}",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "-t",
        "--toml_base",
        help="Base cenplot tracks.",
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-s",
        "--sample_opts",
        help="Sample plot opts JSON as {'sm': {'color': ..., 'display_name': ..., 'shape': ...}}",
        default=None,
        type=str,
    )
    ap.add_argument(
        "--display_order",
        help="Order by display_name. Only used if sample_opts provided.",
        type=str,
        nargs="*",
    )
    ap.add_argument(
        "-o",
        "--outplot_prefix",
        help="Output plot prefix. Creates pdf and png.",
        default=None,
        type=str,
    )
    args = ap.parse_args()

    if args.outplot_prefix:
        outplot_prefix = args.outplot_prefix
    else:
        dr = os.path.dirname(args.infile)
        fname, _ = os.path.splitext(os.path.basename(args.infile))
        outplot_prefix = os.path.join(dr, fname)

    df = pl.read_csv(args.infile, has_header=True, separator="\t").unique()

    if True:
        # Remove primates.
        df = df.filter(~pl.col("qry_sm").str.contains("^m"))

    display_order = []
    if args.sample_opts:
        sample_opts = json.loads(args.sample_opts)
        samples, colors, display_names, shapes = [], [], [], []
        for sample, opts in sample_opts.items():
            samples.append(sample)
            colors.append(opts.get("color", "gray"))
            display_names.append(opts.get("display_name", sample))
            shapes.append(opts.get("shape", "o"))
        sample_opts = {
            "sample": samples,
            "color": colors,
            "display_name": display_names,
            "shape": shapes,
        }
        df_sample_opts = pl.DataFrame(sample_opts)

        df = df.join(
            df_sample_opts, left_on="qry_sm", right_on="sample", how="left"
        ).with_columns(display_name=pl.col("display_name").fill_null(pl.col("qry_sm")))
        if args.display_order:
            df = df.join(
                pl.DataFrame(
                    {
                        # Sort in reverse.
                        "display_name": list(reversed(args.display_order)),
                        "order": range(len(args.display_order)),
                    }
                ),
                on="display_name",
                how="left",
            ).sort(by="order")
            display_order = df["display_name"].unique(maintain_order=True).to_list()
    else:
        df = df.with_columns(
            color=pl.lit("gray"), display_name=pl.col("qry_sm"), shape=pl.lit("o")
        )

    # category should be str col of comma delimited values {ref, other, ignore}
    ref_muts = df.filter(pl.col("category").str.contains("ref"))["mu"]
    other_muts = df.filter(pl.col("category").str.contains("other"))["mu"]

    try:
        tracks, settings = read_tracks(args.toml_base)
        chrom = next(iter(tracks.chroms))
        tracks = tracks.tracks
    except (ValueError, AttributeError) as err:
        print(f"Using defaults due to error: {err}", file=sys.stderr)
        tracks = []
        chrom = df["ref_chrom"][0]
        settings = PlotSettings(title=chrom, dim=(18.0, 8.0))

    ax_mut_idx = max(0, len(tracks) - 1)
    for i, df_part in enumerate(
        df.partition_by(["color", "display_name", "shape"], maintain_order=True)
    ):
        if i != 0:
            pos = TrackPosition.Overlap
            title = None
        else:
            pos = TrackPosition.Relative
            title = "Mutation rate\n(mutations per bp\nper generation)"

        color = df_part["color"][0]
        shape = df_part["shape"][0]
        display_name = df_part["display_name"][0]
        df_part = df_part.select(
            chrom=pl.col("ref_chrom"),
            chrom_st=pl.col("ref_st"),
            chrom_end=pl.col("ref_end"),
            name=pl.col("mu"),
            score=pl.lit(0),
            strand=pl.lit("."),
            thick_st=pl.col("ref_st"),
            thick_end=pl.col("ref_end"),
            item_rgb=pl.col("color"),
        ).sort(by=["chrom", "chrom_st"])

        new_track = Track(
            title=title,
            pos=pos,
            opt=TrackType.Line,
            prop=0.5,
            data=df_part,
            options=LineTrackSettings(
                legend=False,
                marker=shape,
                markersize=3,
                color=color,
                ymin=10e-11,
                ymax=10e-5,
                hide_x=False,
                linestyle="none",
                label=display_name,
                log_scale=True,
            ),
        )
        tracks.append(new_track)

    # Add legend tracks.
    legend_tracks = []
    idx = 0
    for trk in tracks:
        # If overlaps, index unaffected.
        if trk.pos == TrackPosition.Overlap:
            continue
        curr_idx = idx
        idx += 1
        if trk.data.is_empty():
            continue

        n_unique = trk.data["name"].n_unique()
        # One element. Don't add.
        if n_unique == 1:
            print(f"Skipping single track with element: {trk}", file=sys.stderr)
            continue

        if trk.title == "Mutation rate\n(mutations per bp\nper generation)":
            legend_label_order = (
                trk.options.legend_label_order
                if not display_order
                else args.display_order
            )
        else:
            legend_label_order = trk.options.legend_label_order

        legend_tracks.append(
            Track(
                title=None,
                pos=TrackPosition.Relative,
                prop=0.05,
                data=trk.data,
                opt=TrackType.Legend,
                options=LegendTrackSettings(
                    index=curr_idx,
                    legend_ncols=10,
                    legend_label_order=legend_label_order,
                ),
            )
        )

    tracks.extend(legend_tracks)
    # Create fig and axes
    fig, axes, _ = plot_tracks(tracks=tracks, chrom=chrom, settings=settings)
    # Get axes.
    # Create secondary axis with shared x
    ax_mut: matplotlib.axes.Axes = axes[ax_mut_idx, 0]
    ax_mut_yaxis: matplotlib.axes.Axes = ax_mut.twinx()

    # Log scale.
    ax_mut_yaxis.set_yscale("log")

    # Add lines
    mean_mut = ref_muts.mean()
    if not mean_mut:
        mean_mut = sys.float_info.min

    other_mut = other_muts.mean()
    if not other_mut or other_mut == 0.0:
        other_mut = 1.0

    assert isinstance(mean_mut, float) and isinstance(other_mut, float), (
        "Invalid type for mutation rate."
    )
    ax_mut.axhline(mean_mut, linestyle="dotted", color="red")
    ax_mut.axhline(other_mut, linestyle="dotted", color="black")

    # Add median ticks.
    second_yticks, second_yticklabels = [], []
    second_yticks.append(mean_mut)
    second_yticks.append(other_mut)
    second_yticklabels.append(f"{mean_mut:.2e}")
    second_yticklabels.append(f"{other_mut:.2e}")
    ax_mut_yaxis.set_yticks(second_yticks, second_yticklabels)
    ax_mut_yaxis.get_yticklabels()[0].set_color("red")
    # Reset ylim
    ax_mut_yaxis.set_ylim(ax_mut.get_ylim())

    rel_mut_rate = round(mean_mut / other_mut)
    print(f"{mean_mut}\t{other_mut}\t{rel_mut_rate}", file=args.outfile)

    fig.savefig(f"{outplot_prefix}.png", bbox_inches="tight", dpi=600)
    fig.savefig(f"{outplot_prefix}.pdf", bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
