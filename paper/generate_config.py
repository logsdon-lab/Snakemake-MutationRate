import os
import re
import sys
import glob
import json
import yaml
import tempfile
import tomllib
import subprocess
import argparse
import polars as pl
from concurrent.futures import Future, ProcessPoolExecutor
from typing import Any

NHP_DIV_TIMES = {
    "mGorGor1": 10700000,
    "mPanPan1": 6200000,
    "mPanTro3": 6200000,
    "mPonAbe1": 19100000,
    "mPonPyg2": 19100000,
}

RGX_CHROM = re.compile(r"chr([0-9XY]+)")
RGX_CHROM_NHP = re.compile(r"hsa([0-9XY]+)")
RGX_HAP_NHP = re.compile(r"(hap1|hap2|mat|pat)")
RGX_COORDS = r"^(?<ctg>.+):(?<ctg_st>\d+)-(?<ctg_end>\d+)$"

SUPER_POP_COLORS = {
    "#7FA771": "EAS",
    "#E49C1C": "AFR",
    "#855BC0": "SAS",
    "#CA5150": "AMR",
    "#548CB8": "EUR",
}
CT_HEX = "#DDDDDD"
DEF_PATTERNS = {
    "asat": {"pattern": r"ALR", "color": "#58245B"},
    "bsat": {"pattern": r"BSR", "color": "#3A3A3A"},
    "gsat": {"pattern": r"GSAT", "color": "#3A3A3A"},
    "hsat1A": {"pattern": r"SAR", "color": "#3A3A3A"},
    "hsat2": {"pattern": r"^HSATII$", "color": "#3A3A3A"},
    "hsat1B": {"pattern": r"^HSATI$", "color": "#3A3A3A"},
    "hsat3": {"pattern": r"(CATTC)n|(GAATG)n", "color": "#3A3A3A"},
    "ct": {"pattern": None, "color": CT_HEX},
}
OUTPUT_COLS = [
    "chrom",
    "chrom_st",
    "chrom_end",
    "name",
    "score",
    "strand",
    "thick_st",
    "thick_end",
    "item_rgb",
]


def format_rm(df: pl.DataFrame):
    df = (
        df.with_columns(
            ctg_st=pl.col("chrom").str.extract(r":(\d+)-").fill_null(0).cast(pl.Int64),
            strand=pl.when(pl.col("strand") == "C")
            .then(pl.lit("-"))
            .otherwise(pl.lit("+")),
        )
        .with_columns(
            chrom_st=pl.col("chrom_st") + pl.col("ctg_st"),
            chrom_end=pl.col("chrom_end") + pl.col("ctg_st"),
        )
        .with_columns(
            thick_st=pl.col("chrom_st"),
            thick_end=pl.col("chrom_end"),
            score=pl.lit(0),
        )
    )

    patterns = DEF_PATTERNS

    # Remove repeats not matched.
    # Add ct to fill in the rest.
    expr_name = pl.when(pl.lit(False)).then(None)
    expr_item_rgb = pl.when(pl.lit(False)).then(None)
    for sat, pat in patterns.items():
        rgx_pattern = pat.get("pattern")
        color = pat.get("color")

        if not rgx_pattern or not color:
            continue

        expr_name = expr_name.when(pl.col("rtype").str.contains(rgx_pattern)).then(
            pl.lit(sat)
        )
        expr_item_rgb = expr_item_rgb.when(
            pl.col("rtype").str.contains(rgx_pattern)
        ).then(pl.lit(color))
    expr_name = expr_name.otherwise(pl.col("rtype"))
    expr_item_rgb = expr_item_rgb.otherwise(None)

    df_sat = (
        df.with_columns(name=expr_name, item_rgb=expr_item_rgb)
        .select(OUTPUT_COLS)
        .drop_nulls()
    )
    df_ct = (
        df_sat.get_column("chrom")
        .unique()
        .str.extract_groups(RGX_COORDS)
        .to_frame(name="mtch_ctg")
        .unnest("mtch_ctg")
        .with_columns(
            chrom=pl.col("ctg") + ":" + pl.col("ctg_st") + "-" + pl.col("ctg_end"),
        )
        .cast({"ctg_st": pl.Int64, "ctg_end": pl.Int64})
        .with_columns(
            chrom_st=pl.col("ctg_st"),
            chrom_end=pl.col("ctg_end"),
            name=pl.lit("ct"),
            score=pl.lit(0),
            strand=pl.lit("."),
            thick_st=pl.col("ctg_st"),
            thick_end=pl.col("ctg_end"),
            item_rgb=pl.lit(CT_HEX),
        )
        .select(OUTPUT_COLS)
    )

    df = pl.concat((df_sat, df_ct)).sort(by=["chrom", "chrom_st"])
    return df, df_sat


def generate_cfg(
    args,
    ref_row: dict[str, Any],
    df_ref_arr_len: pl.DataFrame,
    df_ref_cdr: pl.DataFrame,
    output_ref_bed_dir: str,
    output_ref_fa_dir: str,
    output_qry_fa_dir: str,
    output_qry_bed_dir: str,
    output_cfg_dir: str,
    df_qry_clade: pl.DataFrame,
    sample_colors: dict[str, str],
    # This is not good code.
    ref_data_sources: dict[str, str],
):
    with open(args.cfg_template, "rt") as fh:
        config = yaml.safe_load(fh)
    chrom = ref_row["chrom"]
    chrom_name = ref_row["chrom_name"]
    clade_n = ref_row["clade_n"]
    clade_name = f"{chrom_name}@{clade_n}"
    print(f"Running: {clade_name}", file=sys.stderr)

    sample = ref_row["sm"]
    ref_chrom_coords = f"{ref_row['chrom']}:{ref_row['st']}-{ref_row['end']}"
    df_ref_cdr = df_ref_cdr.select(
        "chrom",
        "st",
        "end",
        name=pl.lit("cdr"),
        score=pl.lit("0"),
        strand=pl.lit("."),
        tst=pl.col("st"),
        tend=pl.col("end"),
        item_rgb=pl.lit("0,0,0"),
    )

    # RM
    df_rm = pl.read_csv(
        ref_data_sources["rm"].format(sm=sample, chrom_coords=ref_chrom_coords),
        columns=[0, 4, 5, 6, 8, 9, 10],
        separator="\t",
        has_header=False,
        new_columns=[
            "score",
            "chrom",
            "chrom_st",
            "chrom_end",
            "strand",
            "rtype",
            "rclass",
        ],
        truncate_ragged_lines=True,
    )
    df_rm, df_sat = format_rm(df_rm)
    df_sat = df_sat.with_columns(pl.col("chrom").str.extract(r"(.*?):"))

    # HOR
    df_hor = pl.read_csv(
        ref_data_sources["hor"].format(sm=sample, chrom_coords=ref_chrom_coords),
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom",
            "st",
            "end",
            "name",
            "score",
            "strand",
            "tst",
            "tend",
            "item_rgb",
        ],
    )

    # ModDotPlot
    df_mdp = pl.read_csv(
        ref_data_sources["mdp"].format(sm=sample, chrom_coords=ref_chrom_coords),
        separator="\t",
        has_header=True,
    )
    print(f"Reference data filtered: {clade_name}", file=sys.stderr)

    ref_bedfiles = {
        "rm": df_rm,
        "hor": df_hor,
        "cdr": df_ref_cdr,
        "mdp": df_mdp,
        "hor_array": df_ref_arr_len,
    }
    ref_bedfile_paths = {}
    for dtype, df_dtype in ref_bedfiles.items():
        bed_path = os.path.join(output_ref_bed_dir, f"{clade_name}_{dtype}.bed")
        # Adjust st/end coordinates. Remove coodinates in name.
        if dtype == "mdp":
            df_dtype = (
                df_dtype.with_columns(
                    ctg_st=pl.col("#query_name").str.extract(r":(\d+)-").cast(pl.Int64)
                )
                .with_columns(
                    pl.col("#query_name").str.extract(r"(.*?):"),
                    pl.col("reference_name").str.extract(r"(.*?):"),
                    pl.col("query_start") + pl.col("ctg_st"),
                    pl.col("query_end") + pl.col("ctg_st"),
                    pl.col("reference_start") + pl.col("ctg_st"),
                    pl.col("reference_end") + pl.col("ctg_st"),
                )
                .drop("ctg_st")
            )
        elif dtype == "hor":
            df_dtype = (
                df_dtype.with_columns(
                    ctg_st=pl.col("chrom").str.extract(r":(\d+)-").cast(pl.Int64)
                )
                .with_columns(
                    pl.col("chrom").str.extract(r"(.*?):"),
                    pl.col("st") + pl.col("ctg_st"),
                    pl.col("end") + pl.col("ctg_st"),
                )
                .drop("ctg_st")
            )
        else:
            df_dtype = df_dtype.with_columns(pl.col("chrom").str.extract(r"(.*?):"))

        df_dtype.write_csv(
            bed_path,
            separator="\t",
            include_header=False,
        )
        ref_bedfile_paths[dtype] = os.path.abspath(bed_path)

    time = ref_row["time"] * 1_000_000

    # Write ref bed.
    ref_bed = os.path.join(output_ref_bed_dir, f"{clade_name}.bed")
    with open(ref_bed, "wt") as fh:
        values = list(ref_row.values())
        # Sub 1 so match data.
        values[1] -= 1
        print("\t".join(str(v) for v in values), file=fh)

    asm_fa = glob.glob(os.path.join(args.fasta_dir, f"{sample}*.fa"))[0]

    ref_fa = os.path.join(output_ref_fa_dir, f"{clade_name}.fa")
    if not os.path.exists(ref_fa) or os.path.getsize(ref_fa) == 0:
        with open(ref_fa, "w") as fh, tempfile.NamedTemporaryFile("wt") as tfh:
            # Get whole contig.
            tfh.write(ref_row["chrom"] + "\n")
            tfh.flush()
            subprocess.run(("seqtk", "subseq", asm_fa, tfh.name), stdout=fh, check=True)

    # Write satellite bed
    bed_other_ignore = os.path.join(
        output_ref_bed_dir, f"{clade_name}_other_ignore.bed"
    )
    df_sat.write_csv(bed_other_ignore, separator="\t", include_header=False)

    config["output_dir"] = os.path.join("results", clade_name)
    config["log_dir"] = os.path.join("logs", clade_name)
    config["benchmark_dir"] = os.path.join("benchmarks", clade_name)
    config["reference"][0]["name"] = chrom
    config["reference"][0]["path"] = os.path.abspath(ref_fa)
    config["reference"][0]["bed"] = os.path.abspath(ref_bed)
    config["reference"][0]["bed_comparison"] = os.path.abspath(
        ref_bedfile_paths["hor_array"]
    )
    config["reference"][0]["bed_other_ignore"] = os.path.abspath(bed_other_ignore)

    # Write to toml file here.
    output_plot_format = os.path.join(output_ref_bed_dir, f"{clade_name}.yaml")
    with (
        open(args.plot_template, "rb") as fh,
        open(output_plot_format, "wt") as ofh,
    ):
        cfg = tomllib.load(fh)
        cfg["settings"]["title"] = ref_chrom_coords
        # Update paths
        for trk in cfg["tracks"]:
            path = trk["path"]
            trk["path"] = os.path.abspath(ref_bedfile_paths[path])

        yaml.safe_dump(cfg, ofh)

    config["reference"][0]["plot_format"] = os.path.abspath(output_plot_format)

    qry_fasta = []
    for row in df_qry_clade.iter_rows(named=True):
        qry_full_chrom = f"{row['chrom']}:{row['st']}-{row['end']}"
        qry_chrom, qry_st, qry_end = row["chrom"], row["st"], row["end"]
        qry_fa = os.path.join(output_qry_fa_dir, f"{qry_full_chrom}.fa")
        qry_bed = os.path.join(output_qry_bed_dir, f"{qry_full_chrom}.bed")

        with open(qry_bed, "wt") as fh:
            fh.write(f"{qry_chrom}\t{qry_st}\t{qry_end}\n")

        qry_sample = qry_chrom.split("_")[0]
        qry_asm_fa = glob.glob(os.path.join(args.fasta_dir, f"{qry_sample}*.fa"))[0]

        if not os.path.exists(qry_fa) or os.path.getsize(qry_fa) == 0:
            with open(qry_fa, "w") as fh:
                subprocess.run(
                    ["seqtk", "subseq", qry_asm_fa, qry_bed], stdout=fh, check=True
                )

        qry_fasta.append(qry_fa)

    sample_config = []
    for fasta in qry_fasta:
        name, _ = os.path.splitext(os.path.basename(fasta))
        fasta = os.path.abspath(fasta)
        sm, _, _ = name.split("_")
        color = sample_colors[sm]
        display_name = SUPER_POP_COLORS[color]
        sm_config = {
            "name": name,
            "path": fasta,
            "mm2_opts": "-x asm20 -K 8G --eqx -Y -t 8 -r 500000,500000",
            "divergence_time": time,
            "color": color,
            "display_name": display_name,
            "shape": "o",
        }
        sample_config.append(sm_config)

    print(f"Writing config: {clade_name}", file=sys.stderr)

    config["samples"] = sample_config
    with open(os.path.join(output_cfg_dir, f"{clade_name}.yaml"), "wt") as fh:
        yaml.safe_dump(config, fh)


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument(
        "-t",
        "--cfg_template",
        help="Config template for mutation rate workflow.",
        default="/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/config/config_template.yaml",
    )
    ap.add_argument(
        "-r",
        "--ref_bed",
        help="All reference bed coordinates.",
        default="/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/data/all_ref.bed",
    )
    ap.add_argument(
        "-q",
        "--qry_bed",
        help="All query bed coordinates.",
        default="/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/data/all_query.bed",
    )
    ap.add_argument(
        "-f",
        "--fasta_dir",
        help="Fasta dir with both reference and query assemblies.",
        default="/project/logsdon_shared/projects/HGSVC3/new_65_asms_renamed/",
    )
    ap.add_argument(
        "--ref_data_sources",
        help="Reference annotation data sources.",
        default="/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/config/ref_sources.json",
    )
    ap.add_argument(
        "--color_code",
        help="Color codes for samples.",
        default="/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/config/sample_colors.tsv",
    )
    ap.add_argument(
        "--plot_template",
        help="Cenplot template.",
        default="/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/config/base.toml",
    )
    ap.add_argument("-o", "--outdir", default="input")
    args = ap.parse_args()

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    sample_colors = dict(
        pl.read_csv(
            args.color_code,
            separator="\t",
            has_header=False,
            new_columns=["sm", "color"],
        ).iter_rows()
    )

    # HG02953_chr1_haplotype1-0000013	16773102	25725928	0.07	HG02953	chr1	3	HG02953_chr1_haplotype1-0000013:17273102-25225928
    df_ref = pl.read_csv(
        args.ref_bed,
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom",
            "st",
            "end",
            "time",
            "sm",
            "chrom_name",
            "clade_n",
            "full_chrom",
        ],
    )
    df_qry = pl.read_csv(
        args.qry_bed,
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom",
            "st",
            "end",
            "original_chrom",
            "clade",
        ],
    )
    with open(args.ref_data_sources, "rt") as fh:
        ref_data_sources = json.load(fh)

    df_arr_len = pl.read_csv(
        ref_data_sources["hor_arr"],
        separator="\t",
        has_header=False,
        new_columns=["chrom", "st", "end", "length"],
    )
    df_cdr = pl.read_csv(
        ref_data_sources["cdr"],
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom",
            "st",
            "end",
        ],
    )
    output_ref_fa_dir = os.path.join(outdir, "fa", "ref_new")
    output_ref_bed_dir = os.path.join(outdir, "bed", "ref_new")
    output_qry_fa_dir = os.path.join(outdir, "fa", "qry_new")
    output_qry_bed_dir = os.path.join(outdir, "bed", "qry_new")
    output_cfg_dir = os.path.join(outdir, "cfg")
    for dr in (
        output_ref_fa_dir,
        output_ref_bed_dir,
        output_qry_fa_dir,
        output_qry_bed_dir,
        output_cfg_dir,
    ):
        os.makedirs(dr, exist_ok=True)

    with ProcessPoolExecutor(24, max_tasks_per_child=1) as pool:
        futures: list[Future] = []
        for ref_row in df_ref.iter_rows(named=True):
            chrom_name = ref_row["chrom_name"]
            clade_n = ref_row["clade_n"]
            clade_name = f"{chrom_name}@{clade_n}"

            # generate_cfg(
            #     args,
            #     ref_row,
            #     df_arr_len.filter(pl.col("chrom").str.contains(ref_row['chrom'])),
            #     df_cdr.filter(pl.col("chrom").str.contains(ref_row['chrom'])),
            #     output_ref_bed_dir,
            #     output_ref_fa_dir,
            #     output_qry_fa_dir,
            #     output_qry_bed_dir,
            #     output_cfg_dir,
            #     df_qry.filter(pl.col("clade") == clade_name),
            #     sample_colors,
            #     ref_data_sources
            # )
            future = pool.submit(
                generate_cfg,
                args,
                ref_row,
                df_arr_len.filter(pl.col("chrom").str.contains(ref_row["chrom"])),
                df_cdr.filter(pl.col("chrom").str.contains(ref_row["chrom"])),
                output_ref_bed_dir,
                output_ref_fa_dir,
                output_qry_fa_dir,
                output_qry_bed_dir,
                output_cfg_dir,
                df_qry.filter(pl.col("clade") == clade_name),
                sample_colors,
                ref_data_sources,
            )
            futures.append(future)

    for future in futures:
        if future.cancelled():
            print(future.exception(), fh=sys.stderr)
        else:
            _ = future.result()


if __name__ == "__main__":
    raise SystemExit(main())
