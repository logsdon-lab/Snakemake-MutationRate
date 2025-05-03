
DIVERGENCE_TIMES = json.dumps({
    sample: cfg["divergence_time"]
    for sample, cfg in SAMPLES.items()
})

rule calculate_mutation_rate:
    input:
        script="workflow/scripts/calculate_mutation_rate.py",
        div=rules.aggregate_divergence_times.output
    output:
        mut_tsv=join(OUTPUT_DIR, "mutation_rate", "{ref}", "{rgn}.tsv"),
    params:
        sample_divergence_times=DIVERGENCE_TIMES
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "calculate_mutation_rate_{ref}_{rgn}.log"),
    shell:
        """
        python {input.script} -i {input.div} -d '{params.sample_divergence_times}' > {output} 2> {log}
        """

rule plot_mutation_rate:
    input:
        script="workflow/scripts/plot_mutation_rate.py",
        mut_tsv=rules.calculate_mutation_rate.output,
    output:
        mut_plt=join(OUTPUT_DIR, "mutation_rate", "{ref}", "{rgn}.pdf"),
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "plot_mutation_rate_{ref}_{rgn}.log"),    
    shell:
        """
        python {input.script} -i {input.mut_tsv} -o {output} 2> {log}
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
            )
            for ref in REFERENCE_REGIONS.keys()
        ],
        