def extract_sm_seq_cmd(wc, input, output):
    if not input.bed:
        return f"ln -s {input.fa} {output.fa}"

    if isinstance(input.bed, str):
        bed_cmd = f"{input.bed}"
    elif isinstance(SAMPLES[wc.sm], list):
        bed_cmd = f"<(cat {input.bed})"
    else:
        raise ValueError(f"Sample {wc.sm} bedfile is not valid type.")

    return f"seqtk subseq {input.fa} {bed_cmd} > {output.fa}"


rule extract_sm_seq:
    input:
        fa=lambda wc: SAMPLES[wc.sm]["path"],
        bed=lambda wc: SAMPLES[wc.sm].get("bed", []),
    output:
        fa=join(OUTPUT_DIR, "fastas", "{sm}", "subset.fa"),
    params:
        cmd=lambda wc, input, output: extract_sm_seq_cmd(wc, input, output),
    conda:
        "../envs/tools.yaml"
    shell:
        """
        {params.cmd}
        """


ALN_CFG = {
    "ref": {ref: cfg_ref["path"] for ref, cfg_ref in REFERENCES.items()},
    "sm": {qry: rules.extract_sm_seq.output for qry, cfg_qry in SAMPLES.items()},
    "temp_dir": join(OUTPUT_DIR, "temp"),
    "output_dir": OUTPUT_DIR,
    "logs_dir": LOGS_DIR,
    "benchmarks_dir": BENCHMARK_DIR,
    "aln_threads": config["threads"],
    "aln_mem": config["mem"],
    "aln_filter_flag": 260,
    # https://github.com/koisland/asm-to-reference-alignment/blob/remove_sm_num_index/config/clint.yaml
    "mm2_opts": {
        qry: cfg_qry.get("mm2_opts", "-x asm20 -Y --eqx -K 8G")
        for qry, cfg_qry in SAMPLES.items()
    },
}


# Align assemblies to reference.
module align_asm_to_ref:
    snakefile:
        github(
            "koisland/asm-to-reference-alignment",
            path="workflow/Snakefile",
            branch="feature/mm2-sample-params",
        )
    config:
        ALN_CFG


use rule * from align_asm_to_ref as asm_ref_*


rule alignment_idx:
    input:
        bam=rules.asm_ref_alignment.output,
    output:
        bam=join(OUTPUT_DIR, "{ref}", "bam", "{sm}.bam"),
        bam_idx=join(OUTPUT_DIR, "{ref}", "bam", "{sm}.bam.bai"),
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "asm_{sm}_{ref}_alignment_idx.log"),
    shell:
        """
        samtools sort {input.bam} -o {output.bam} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule align_all:
    input:
        expand(rules.asm_ref_all.input, ref=REFERENCES.keys(), sm=SAMPLES.keys()),
        expand(
            rules.alignment_idx.output,
            ref=REFERENCES.keys(),
            sm=SAMPLES.keys(),
        ),
