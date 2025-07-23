import os
import sys
import json
import argparse
import polars as pl

IN_COLS = ["ref", "qry", "gc_ref", "gc_qry", "gc", "tn93_div"]
OUT_COLS = [
    "ref_chrom",
    "ref_st",
    "ref_end",
    "mu",
    "qry_chrom",
    "qry_st",
    "qry_end",
    "ref_sm",
    "qry_sm",
]
REF_SM_DIV_TIMES_COLS = ["ref_chrom", "qry_chrom", "repl_time"]


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
        "--ref_sm_ctg_div_times",
        default=None,
        help=f"TSV file of generation times of sample relative to reference for each contig. Overrides divergence times. Format: {REF_SM_DIV_TIMES_COLS}",
        type=str,
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

    if args.ref_sm_ctg_div_times:
        df_ref_sm_ctg_div_times = pl.read_csv(
            args.ref_sm_ctg_div_times,
            has_header=False,
            separator="\t",
            new_columns=REF_SM_DIV_TIMES_COLS,
        )
    else:
        df_ref_sm_ctg_div_times = pl.DataFrame(schema=REF_SM_DIV_TIMES_COLS)

    df_divergence_times = pl.DataFrame(divergence_times).unpivot(
        variable_name="sm", value_name="divergence_time"
    )

    df = (
        pl.scan_csv(args.infile, separator="\t", has_header=False, new_columns=IN_COLS)
        .with_columns(
            # Expects ^{sm}~{chrom}:{st}-{end}$
            mtch_ref=pl.col("ref").str.extract_groups(
                r"^(?<ref_sm>.*?)~(?<ref_chrom>.*?):.*?(?<ref_st>\d+)-(?<ref_end>\d+)$"
            ),
            mtch_qry=pl.col("qry").str.extract_groups(
                r"^(?<qry_sm>.*?)~(?<qry_chrom>.*?):.*?(?<qry_st>\d+)-(?<qry_end>\d+)$"
            ),
        )
        .unnest("mtch_ref")
        .unnest("mtch_qry")
        .join(
            df_divergence_times.lazy(), left_on=["qry_sm"], right_on=["sm"], how="left"
        )
        .rename({"divergence_time": "qry_sm_divergence_time"})
        .collect()
    )
    # Replace divergence times if provided.
    if args.ref_sm_ctg_div_times:
        df = (
            df.join(df_ref_sm_ctg_div_times, on=["ref_chrom", "qry_chrom"], how="left")
            .with_columns(
                qry_sm_divergence_time=pl.when(~pl.col("repl_time").is_null())
                .then(pl.col("repl_time"))
                .otherwise(pl.col("qry_sm_divergence_time"))
            )
            .drop_nulls()
        )

    df = df.with_columns(
        mu=pl.when(pl.col("qry_sm_divergence_time") == 0)
        .then(compute_mu_poly("tn93_div", Ne=Ne))
        .otherwise(compute_mu_div("tn93_div", "qry_sm_divergence_time", Ne=Ne, gen=gen))
    ).select(OUT_COLS)

    df.write_csv(args.outfile, separator="\t", include_header=True)


if __name__ == "__main__":
    raise SystemExit(main())
