import os
import re
import sys
import glob
import yaml
import tomllib
import polars as pl

NHP_DIV_TIMES = {
    "mGorGor1": 10700000,
    "mPanPan1": 6200000,
    "mPanTro3": 6200000,
    "mPonAbe1": 19100000,
    "mPonPyg2": 19100000,
}
DATA_SOURCES = {
    "hor_arr": lambda chrom_name: "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/all_array_length.bed",
    "rm": lambda chrom_name: f"/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/RM/all_cens_{chrom_name}.annotation.fa.out",
    "hor": lambda chrom_name: f"/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/HOR/{chrom_name}_AS-HOR_stv_row.all.bed",
    "cdr": lambda chrom_name: "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/all_cdrs.bed",
    "mdp": lambda chrom_name: "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/1DModplot.bed",
}

RGX_CHROM = re.compile(r"chr([0-9XY]+)")
RGX_CHROM_NHP = re.compile(r"hsa([0-9XY]+)")

COLOR_SCALE = [
    ((0.0, 90.0), "#4b3991"),
    ((90.0, 97.5), "#2974af"),
    ((97.5, 97.75), "#4a9da8"),
    ((97.75, 98.0), "#57b894"),
    ((98.0, 98.25), "#9dd893"),
    ((98.25, 98.5), "#e1f686"),
    ((98.5, 98.75), "#ffffb2"),
    ((98.75, 99.0), "#fdda79"),
    ((99.0, 99.25), "#fb9e4f"),
    ((99.25, 99.5), "#ee5634"),
    ((99.5, 99.75), "#c9273e"),
    ((99.75, 100.0), "#8a0033"),
]
super_population_color = {
    "#7FA771": "EAS",
    "#E49C1C": "AFR",
    "#855BC0": "SAS",
    "#CA5150": "AMR",
    "#548CB8": "EUR",
}


def convert_ident_to_color(
    df: pl.DataFrame, color_scale: list[tuple[tuple[float, float], str]]
) -> pl.DataFrame:
    color_expr = None
    rng_expr = None
    for rng, color in color_scale:
        if not isinstance(color_expr, pl.Expr):
            color_expr = pl.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(color))
            rng_expr = pl.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(f"{rng[0]}-{rng[1]}"))
        else:
            color_expr = color_expr.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(color))
            rng_expr = rng_expr.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(f"{rng[0]}-{rng[1]}"))

    if isinstance(color_expr, pl.Expr):
        color_expr = color_expr.otherwise(None)
    else:
        color_expr = pl.lit(None)
    if isinstance(rng_expr, pl.Expr):
        rng_expr = rng_expr.otherwise(None)
    else:
        rng_expr = pl.lit(None)

    return df.with_columns(
        name=rng_expr,
        color=color_expr,
    )


def main():
    # "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/exp/ref_clade_based_new/config_template.yaml"
    template = sys.argv[1]
    ref_fa = sys.argv[2]
    # "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/exp/ref_clade_based_new/ref_clades/chr1@3.bed"
    ref_bed = sys.argv[3]
    # "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/exp/ref_clade_based_new/others/chr1@3"
    qry_dir = sys.argv[4]
    # "/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/data/cens_primates"
    nhp_dir = sys.argv[5]
    # /project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/exp/ref_clade_based_new/sample_super_color.xls
    sample_color_code = sys.argv[6]
    # Plot format.
    base_plot_format = sys.argv[7]

    outfile = sys.stdout

    with open(template, "rt") as fh:
        config = yaml.safe_load(fh)

    sample_colors = dict(
        pl.read_csv(
            sample_color_code,
            separator="\t",
            has_header=False,
            new_columns=["sm", "color"],
        ).iter_rows()
    )

    # HG02953_chr1_haplotype1-0000013	16773102	25725928	0.07	HG02953	chr1	3	HG02953_chr1_haplotype1-0000013:17273102-25225928
    df_ref_row = pl.read_csv(
        ref_bed,
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
    chrom, st, end, time, sm, chrom_name, clade_n, full_chrom = df_ref_row.row(0)

    ref_bed_dir = os.path.dirname(ref_bed)
    ref_bed_dir = os.path.abspath(ref_bed_dir)

    clade_name = f"{chrom_name}@{clade_n}"

    # HOR ort
    arr_len_bed = DATA_SOURCES["hor_arr"](chrom_name)
    df_arr_len = pl.read_csv(
        arr_len_bed,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "st", "end", "length"],
    ).filter(pl.col("chrom").str.contains(full_chrom))

    # RM
    rm_bed = DATA_SOURCES["rm"](chrom_name)
    df_rm = pl.read_csv(
        rm_bed,
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
    ).filter(pl.col("chrom").str.contains(full_chrom))

    # HOR
    hor_bed = DATA_SOURCES["hor"](chrom_name)
    df_hor = pl.read_csv(
        hor_bed,
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
    ).filter(pl.col("chrom").str.contains(full_chrom))

    # CDR
    cdr_bed = DATA_SOURCES["cdr"](chrom_name)
    df_cdr = (
        pl.read_csv(
            cdr_bed,
            separator="\t",
            has_header=False,
            new_columns=[
                "chrom",
                "st",
                "end",
            ],
        )
        .filter(pl.col("chrom").str.contains(full_chrom))
        .select(
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
    )

    # ModDotPlot
    mdp_bed = DATA_SOURCES["mdp"](chrom_name)
    df_mdp = pl.read_csv(
        mdp_bed,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "st", "end", "percent_identity_by_events"],
    ).filter(pl.col("chrom").str.contains(full_chrom))
    df_mdp = convert_ident_to_color(df_mdp, COLOR_SCALE).select(
        "chrom",
        "st",
        "end",
        "name",
        score=pl.lit("0"),
        strand=pl.lit("."),
        tst=pl.col("st"),
        tend=pl.col("end"),
        item_rgb=pl.col("color"),
    )

    ref_bedfiles = {
        "rm": df_rm,
        "hor": df_hor,
        "cdr": df_cdr,
        "mdp": df_mdp,
        "hor_array": df_arr_len,
    }
    ref_bedfile_paths = {}
    for dtype, df_dtype in ref_bedfiles.items():
        bed_path = os.path.join(ref_bed_dir, f"{clade_name}_{dtype}.bed")
        df_dtype = df_dtype.with_columns(pl.col("chrom").str.extract(r"(.*?):"))
        df_dtype.write_csv(
            bed_path,
            separator="\t",
            include_header=False,
        )
        ref_bedfile_paths[dtype] = bed_path

    time = time * 1_000_000

    config["output_dir"] = os.path.join("results", clade_name)
    config["log_dir"] = os.path.join("logs", clade_name)
    config["benchmark_dir"] = os.path.join("benchmarks", clade_name)
    config["reference"][0]["name"] = chrom
    config["reference"][0]["path"] = os.path.abspath(ref_fa)
    config["reference"][0]["bed"] = os.path.abspath(ref_bed)
    config["reference"][0]["bed_comparison"] = ref_bedfile_paths["hor_array"]

    # Write to toml file here.
    output_plot_format = os.path.join(ref_bed_dir, f"{clade_name}.yaml")
    with (
        open(base_plot_format, "rb") as fh,
        open(output_plot_format, "wt") as ofh,
    ):
        cfg = tomllib.load(fh)
        cfg["settings"]["title"] = full_chrom
        # Update paths
        for trk in cfg["tracks"]:
            path = trk["path"]
            # Set colorscale order.
            if trk["path"] == "mdp":
                trk["options"]["legend_label_order"] = [
                    f"{cs[0]}-{cs[1]}" for cs, _ in COLOR_SCALE
                ]
            trk["path"] = ref_bedfile_paths[path]

        yaml.safe_dump(cfg, ofh)

    config["reference"][0]["plot_format"] = output_plot_format

    sample_config = []
    for fasta in glob.glob(os.path.join(qry_dir, "*.fa")):
        name, _ = os.path.splitext(os.path.basename(fasta))
        fasta = os.path.abspath(fasta)
        sm, _, _ = name.split("_")
        color = sample_colors[sm]
        display_name = super_population_color[color]
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

    chrom_name = re.search(RGX_CHROM, chrom).group(1)
    for fa in glob.glob(os.path.join(nhp_dir, "*", "*.fa")):
        species = os.path.basename(os.path.dirname(fa))
        species_chrom = os.path.basename(fa)
        species_chrom_name = re.search(RGX_CHROM_NHP, species_chrom).group(1)
        color = sample_colors[species]

        if species_chrom_name != chrom_name:
            continue
        species_sm_config = {
            "name": species,
            "path": fa,
            "mm2_opts": "-x asm20 -K 8G --eqx -Y -t 8 -r 500000,500000",
            "divergence_time": NHP_DIV_TIMES[species],
            "color": color,
            "shape": "^",
        }
        sample_config.append(species_sm_config)

    config["samples"] = sample_config

    yaml.safe_dump(config, outfile)


if __name__ == "__main__":
    raise SystemExit(main())
