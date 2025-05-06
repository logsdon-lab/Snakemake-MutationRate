import argparse
import ast
import sys
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--fai", type=str, required=True)
    ap.add_argument(
        "--rgx_species", type=str, default=r"(.*?)~", help="Species name regex."
    )
    ap.add_argument(
        "--rgx_hap",
        type=str,
        default=r"(mat|pat|haplotype1|haplotype2|hap1|hap2|h1|h2)",
        help="Regex pattern to extract haplotype.",
    )
    ap.add_argument(
        "--rgx_chrom",
        type=str,
        default=r"_(hsa(.*?)|chr(.*?)|cen(.*?))[:_]",
        help="Regex pattern to extract chromosome.",
    )
    ap.add_argument(
        "-c",
        "--rchrom",
        type=str,
        required=True,
        help="Reference chromsome to filter for. Should be [0-9XYM]{1,2} only.",
    )
    ap.add_argument(
        "-l",
        "--length",
        type=lambda x: ast.literal_eval(x),
        default=(5_000, 20_000),
        help="Length range of sequence.",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("wt"),
        default=sys.stdout,
        help="Output list of contigs to keep.",
    )

    args = ap.parse_args()
    fai = args.fai
    rgx_hap = args.rgx_hap
    rgx_chrom = args.rgx_chrom
    rgx_species = args.rgx_species
    # 10, X, etc.
    rchrom = args.rchrom

    df = pl.read_csv(
        fai,
        has_header=False,
        separator="\t",
        columns=[0, 1],
        new_columns=["chrom", "length"],
    )
    df = df.with_columns(
        # Get chrom in name ex. chr5 -> 5
        # For cases like chr5x17 only first 5 captured.
        rchrom=pl.col("chrom").str.extract(rgx_chrom).str.extract(r"([0-9XY]+)"),
        hap=pl.col("chrom").str.extract(rgx_hap),
        species=pl.col("chrom").str.extract(rgx_species),
    )
    # Filter for reference chrom and if multiple mapping to a hap, remove.
    (
        df.filter(
            (pl.col("rchrom") == rchrom) & pl.col("length").is_between(*args.length)
        )
        .filter((pl.col("chrom").count() <= 1).over(["species", "hap", "rchrom"]))
        .select("chrom")
        .write_csv(args.output, include_header=False)
    )


if __name__ == "__main__":
    raise SystemExit(main())
