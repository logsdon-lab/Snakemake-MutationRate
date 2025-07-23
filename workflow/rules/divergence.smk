
rule compile_tn93:
    input:
        src="workflow/scripts/libtn93/python",
    output:
        touch(join(OUTPUT_DIR, "divergence", "compiled_tn93.done")),
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "compile_tn93.log"),
    shell:
        """
        {{ pip install {input.src} || true ;}}  &> {log}
        """


rule split_msa_into_pairs:
    input:
        script="workflow/scripts/split_msa_into_pairs.py",
        libtn93_chkpt=rules.compile_tn93.output,
        fa_msa=rules.msa.output,
    output:
        # Outputs "{p1}+{p2}.fa"
        fa_pair=directory(join(OUTPUT_DIR, "divergence", "{ref}", "{window}")),
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "split_msa_into_pairs_{ref}_{window}.log"),
    shell:
        """
        python {input.script} {input.fa_msa} {wildcards.ref} {output} 2> {log}
        """


rule compute_div_tn93:
    input:
        script="workflow/scripts/computePairwiseDivTN93.py",
        fa_pair=join(OUTPUT_DIR, "divergence", "{ref}", "{window}"),
    output:
        div_tsv=join(OUTPUT_DIR, "divergence", "{ref}", "{window}", "divergence.tsv"),
    log:
        join(LOGS_DIR, "compute_div_tn93_{ref}_{window}.log"),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        python {input.script} -f $(find {input.fa_pair} -name "*.fa") > {output.div_tsv} 2> {log}
        touch {output}
        """


def divergence_times(wc):
    outputs = []
    chrom = get_chrom(wc.rgn)
    if not chrom:
        raise ValueError(f"Region {wc.rgn} must contain chromosome.")

    chkpts = [checkpoints.get_aligned_query_regions.get(**wc, sm=sm) for sm in SAMPLES]
    windows = []
    for chkpt in chkpts:
        wcs = glob_wildcards(
            join(
                dirname(chkpt.output[0]),
                "{window}.fa",
            )
        )

        for window in wcs.window:
            window_chrom = get_chrom(window)
            if window_chrom and window_chrom == chrom:
                windows.append(window)
    results = expand(
        rules.compute_div_tn93.output,
        ref=wc.ref,
        window=windows,
    )
    if not results:
        raise ValueError(
            "No extracted regions for checkpoint get_aligned_query_regions."
        )

    return results


rule aggregate_divergence_times:
    input:
        divergence_times,
    output:
        fa=join(OUTPUT_DIR, "divergence", "{ref}", "{rgn}.tsv"),
    shell:
        """
        for file in {input}; do
            cat $file >> {output}
        done
        """


rule divergence_all:
    input:
        [
            expand(
                rules.aggregate_divergence_times.output,
                ref=ref,
                rgn=REFERENCE_REGIONS[ref],
            )
            for ref in REFERENCE_REGIONS.keys()
        ],
