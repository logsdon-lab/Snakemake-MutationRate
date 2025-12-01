DIVERGENCE_TIMES = json.dumps(
    {sample: cfg.get("divergence_time") for sample, cfg in SAMPLES.items()}
)

PLOT_OPTS = json.dumps(
    {
        sample: {
            "color": cfg.get("color"),
            "shape": cfg.get("shape"),
            "display_name": cfg.get("display_name"),
        }
        for sample, cfg in SAMPLES.items()
    }
)


rule calculate_mutation_rate:
    input:
        script="workflow/scripts/calculate_mutation_rate.py",
        div=rules.aggregate_divergence_times.output,
        ref_sm_ctg_div_times=(
            config["reference_sample_ctg_divergence_times"]
            if config.get("reference_sample_ctg_divergence_times")
            else []
        ),
    output:
        mut_bedpe=temp(join(OUTPUT_DIR, "mutation_rate", "{ref}", "{rgn}.bedpe")),
    params:
        sample_divergence_times=DIVERGENCE_TIMES,
        ref_sm_ctg_div_times=lambda wc, input: (
            f"--ref_sm_ctg_div_times {input.ref_sm_ctg_div_times}"
            if input.ref_sm_ctg_div_times
            else ""
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "calculate_mutation_rate_{ref}_{rgn}.log"),
    shell:
        """
        python {input.script} -i {input.div} -d '{params.sample_divergence_times}' \
        {params.ref_sm_ctg_div_times} > {output} 2> {log}
        """


rule annotate_mutation_rates:
    input:
        mut_bedpe=rules.calculate_mutation_rate.output,
        ref_cmp_bed=lambda wc: REFERENCES[wc.ref].get("bed_comparison", []),
        other_ign_bed=lambda wc: REFERENCES[wc.ref].get("bed_other_ignore", []),
    output:
        mut_bedpe=join(OUTPUT_DIR, "mutation_rate", "{ref}", "{rgn}_annot.bedpe"),
    conda:
        "../envs/tools.yaml"
    params:
        ref_cmp_bed=lambda wc, input: (
            f"{input.ref_cmp_bed}" if input.ref_cmp_bed else "<(echo '')"
        ),
        other_ign_bed=lambda wc, input: (
            f"{input.other_ign_bed}" if input.other_ign_bed else "<(echo '')"
        ),
    log:
        join(LOGS_DIR, "calculate_mutation_rate_{ref}_{rgn}.log"),
    shell:
        """
        # Add header
        head -n1 {input.mut_bedpe} | awk -v OFS="\\t" '{{ print $0, "category"}}' > {output.mut_bedpe}
        # Intersect with headerless BEDPE and add overlap category (reference/ignore), if any.
        # No overlap is other.
        # Then merge the results into a 10th comma-delimited column.
        {{ bedtools intersect \
            -a <(tail -n+2 {input.mut_bedpe}) \
            -b <(cat \
                <(awk -v OFS="\\t" '{{ print $1, $2, $3, "ref" }}' {params.ref_cmp_bed}) \
                <(awk -v OFS="\\t" '{{ print $1, $2, $3, "ignore" }}' {params.other_ign_bed}) \
            ) \
            -loj | \
        cut -f 1-9,13 | \
        awk -v OFS="\\t" '{{ $10=($10 == ".") ? "other" : $10; print }}' | \
        bedtools groupby -i - -g 1-9 -c 10 -o distinct ;}} >> {output.mut_bedpe} 2> {log}
        """


rule plot_mutation_rate:
    input:
        script="workflow/scripts/plot_mutation_rate.py",
        mut_bedpe=rules.annotate_mutation_rates.output,
        plot_format=lambda wc: REFERENCES[wc.ref].get("plot_format", []),
    output:
        mut_plt=join(OUTPUT_DIR, "mutation_rate", "{ref}", "{rgn}.png"),
        mut_plt_pdf=join(OUTPUT_DIR, "mutation_rate", "{ref}", "{rgn}.pdf"),
        mut_rate=join(OUTPUT_DIR, "mutation_rate", "{ref}", "{rgn}_mut.tsv"),
    params:
        plot_opts=PLOT_OPTS,
        plot_format=lambda wc, input: (
            f"-t {input.plot_format}" if input.plot_format else ""
        ),
        output_prefix=lambda wc, output: splitext(output.mut_plt)[0],
        display_order=lambda wc: (
            f"--display_order {" ".join(config["display_order"])}"
            if config.get("display_order")
            else ""
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "plot_mutation_rate_{ref}_{rgn}.log"),
    shell:
        """
        {{ python {input.script} \
        -i {input.mut_bedpe} \
        -o {params.output_prefix} \
        -s '{params.plot_opts}' \
        {params.plot_format} {params.display_order} | \
        awk -v OFS="\\t" '{{ print "{wildcards.ref}", "{wildcards.rgn}", $0 }}';}} \
        > {output.mut_rate} 2> {log}
        """


rule mutation_rate_all:
    input:
        [
            (
                expand(
                    rules.annotate_mutation_rates.output,
                    ref=ref,
                    rgn=REFERENCE_REGIONS[ref],
                ),
                expand(
                    rules.plot_mutation_rate.output,
                    ref=ref,
                    rgn=REFERENCE_REGIONS[ref],
                ),
            )
            for ref in REFERENCE_REGIONS.keys()
        ],
