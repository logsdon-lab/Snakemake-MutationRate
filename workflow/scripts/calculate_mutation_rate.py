import os
import sys
import json
import argparse
import polars as pl

IN_COLS = ["ref", "qry", "gc_ref", "gc_qry", "gc", "tn93_div"]
OUT_COLS = ["ref", "ref_st", "ref_end", "qry", "qry_sm", "qry_sm_divergence_time", "mu"]


def compute_mu_div(col_tn93_div: str, col_div_time: str, Ne: int, gen: int) -> pl.Expr:
    return pl.col(col_tn93_div) / (2 * (pl.col(col_div_time) / gen) + 4 * Ne)


def compute_mu_poly(col_tn93_div: str, Ne: int) -> pl.Expr:
    return pl.col(col_tn93_div) / pl.lit(4 * Ne)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--infile",
        help=f"TSV file of divergence time estimates. Samples should be delimited by '~'. Format: {IN_COLS}",
        required=True,
        type=str,
    )
    ap.add_argument(
        "-d",
        "--divergence_times",
        help="Divergence times by sample. Either a json string or a JSON file.",
        required=True,
        type=str,
    )
    ap.add_argument(
        "--ref_ancestral_pop_size",
        default=10_000,
        help="Ancestral population size of reference.",
        type=int,
    )
    ap.add_argument(
        "--ref_generation_time",
        default=20,
        help="Generation time of reference.",
        type=int,
    )
    ap.add_argument(
        "-o",
        "--outfile",
        help=f"Bedfile with mutation rate per ref-qry pair. Format: {OUT_COLS}",
        type=argparse.FileType("wt"),
        default=sys.stdout,
    )

    args = ap.parse_args()
    Ne = args.ref_ancestral_pop_size
    gen = args.ref_generation_time

    if os.path.exists(args.divergence_times):
        with open(args.divergence_times) as fh:
            divergence_times = json.load(fh)
    else:
        divergence_times = json.loads(args.divergence_times)

    df_divergence_times = pl.DataFrame(divergence_times).unpivot(
        variable_name="sm", value_name="divergence_time"
    )

    df = (
        pl.scan_csv(args.infile, separator="\t", has_header=False, new_columns=IN_COLS)
        .with_columns(
            # Expects ^{sm}~{chrom}:{st}-{end}$
            mtch_ref=pl.col("ref").str.extract_groups(
                r"^(?<ref_sm>.*?)~(?<ref_chrom>.*?):(?<ref_st>\d+)-(?<ref_end>\d+)$"
            ),
            qry_sm=pl.col("qry").str.extract(r"^(.*?)~"),
        )
        .unnest("mtch_ref")
        .join(
            df_divergence_times.lazy(), left_on=["qry_sm"], right_on=["sm"], how="left"
        )
        .rename({"divergence_time": "qry_sm_divergence_time"})
        .collect()
    )
    df = df.with_columns(
        mu=pl.when(pl.col("qry_sm_divergence_time") == 0)
        .then(compute_mu_poly("tn93_div", Ne=Ne))
        .otherwise(compute_mu_div("tn93_div", "qry_sm_divergence_time", Ne=Ne, gen=gen))
    ).select(OUT_COLS)
    df.write_csv(args.outfile, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
