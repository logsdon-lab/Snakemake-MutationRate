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
        mut_tsv=join(OUTPUT_DIR, "mutation_rate", "{ref}", "{rgn}.tsv"),
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


rule plot_mutation_rate:
    input:
        script="workflow/scripts/plot_mutation_rate.py",
        mut_tsv=rules.calculate_mutation_rate.output,
        ref_cmp_bed=lambda wc: REFERENCES[wc.ref].get("bed_comparison", []),
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
        ref_cmp_bed=lambda wc, input: (
            f"--ref_cmp_bed {input.ref_cmp_bed}" if input.ref_cmp_bed else ""
        ),
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "plot_mutation_rate_{ref}_{rgn}.log"),
    shell:
        """
        python {input.script} \
        -i {input.mut_tsv} \
        -o {params.output_prefix} \
        -s '{params.plot_opts}' \
        {params.plot_format} \
        {params.ref_cmp_bed} > {output.mut_rate} 2> {log}
        """


rule aggregate_mutation_rate:
    input:
        lambda wc: expand(
            rules.plot_mutation_rate.output.mut_rate,
            ref=wc.ref,
            rgn=REFERENCE_REGIONS[wc.ref],
        ),
    output:
        mut_rate=join(OUTPUT_DIR, "mutation_rate", "{ref}", "all_mut.tsv"),
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "aggregate_mutation_rate_{ref}.log"),
    shell:
        """
        {{ awk -v OFS="\\t" '{{
            if ($1 != "mu_ref") {{
                sub(".*/", "", FILENAME);
                sub("_mut.tsv", "", FILENAME);
                print FILENAME, $0
            }}
        }}' {input} | \
        sort -k 4,4n ;}} > {output}
        """


rule mutation_rate_all:
    input:
        [
            (
                expand(
                    rules.calculate_mutation_rate.output,
                    ref=ref,
                    rgn=REFERENCE_REGIONS[ref],
                ),
                expand(
                    rules.plot_mutation_rate.output,
                    ref=ref,
                    rgn=REFERENCE_REGIONS[ref],
                ),
                expand(
                    rules.aggregate_mutation_rate.output,
                    ref=ref,
                ),
            )
            for ref in REFERENCE_REGIONS.keys()
        ],
