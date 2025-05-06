
checkpoint get_aligned_query_regions:
    input:
        bam=rules.alignment_idx.output.bam,
        bam_idx=rules.alignment_idx.output.bam_idx,
        bed=lambda wc: REFERENCES[wc.ref]["bed"][wc.rgn],
    output:
        outdir=touch(join(OUTPUT_DIR, "fastas", "msa", "{ref}", "{rgn}_{sm}.done")),
    singularity:
        "docker://eichlerlab/subseqfa:1.0"
    params:
        outdir=lambda wc, output: dirname(output[0]),
    log:
        join(LOGS_DIR, "get_aligned_query_regions_{ref}_{rgn}_{sm}.log"),
    shell:
        """
        mkdir -p "{params.outdir}"
        while IFS='' read -r line; do
            window=$(echo $line | awk '{{print $1":"$2"-"$3}}')
            output_fa="{params.outdir}/${{window}}.fa"
            {{ subseqfa -b -v -r "${{window}}" {input.bam} | \
                awk '{{
                    if (substr($0, 1, 1)==">") {{
                        fname=substr($0, 2)
                        # tilde since we need to reserve _
                        $0=">{wildcards.sm}~"fname
                    }}
                    print $0
                }}' >> "${{output_fa}}" ;}} 2>> {log}
            touch "${{output_fa}}"
        done < {input.bed}
        """


rule filter_regions:
    input:
        script="workflow/scripts/filter_regions.py",
        all_rgns_done=lambda wc: expand(
            rules.get_aligned_query_regions.output,
            ref=wc.ref,
            rgn=REFERENCE_REGIONS[wc.ref],
            sm=SAMPLES,
        ),
        ref_fa=lambda wc: REFERENCES[wc.ref]["path"],
        fa=join(OUTPUT_DIR, "fastas", "msa", "{ref}", "{window}.fa"),
    output:
        fai=join(OUTPUT_DIR, "fastas", "msa", "{ref}", "{window}.fa.fai"),
        filtered_fa=join(OUTPUT_DIR, "fastas", "msa", "{ref}_filtered", "{window}.fa"),
        filtered_fai=join(
            OUTPUT_DIR, "fastas", "msa", "{ref}_filtered", "{window}.fa.fai"
        ),
    params:
        # Expects: ^chr(?):.*?$ in window name from reference bed file
        chrom=lambda wc: get_chrom(wc.window),
        window_bed_record=lambda wc: get_window_bed_record(wc.window),
        rgx_hap=config["regex_sm_hap"],
        rgx_chrom=config["regex_sm_chrom"],
        rgx_species=r"'(.*?)~'",
        length_range="'(5_000,20_000)'",
    conda:
        "../envs/tools.yaml"
    log:
        join(LOGS_DIR, "filter_regions_{ref}_{window}.log"),
    shell:
        """
        {{ samtools faidx {input.fa} || true ;}} &> /dev/null
        seqtk subseq {input.fa} \
            <(python {input.script} -i {output.fai} --rgx_hap {params.rgx_hap} --rgx_chrom {params.rgx_chrom} --rgx_species {params.rgx_species} -c {params.chrom} -l {params.length_range}) > {output.filtered_fa} 2> {log}
        seqtk subseq {input.ref_fa} <(printf "{params.window_bed_record}") | \
            awk '{{
                if (substr($0, 1, 1)==">") {{
                    fname=substr($0, 2)
                    $0=">{wildcards.ref}~"fname
                }}
                print $0
            }}'>> {output.filtered_fa} 2> {log}
        {{ samtools faidx {output.filtered_fa} || true ;}} &> /dev/null
        touch {output}
        """


rule msa:
    input:
        fa=rules.filter_regions.output.filtered_fa,
    output:
        fa=join(OUTPUT_DIR, "msa", "{ref}", "{window}.fa"),
    conda:
        "../envs/tools.yaml"
    threads: 4
    resources:
        mem="60GB",
    log:
        join(LOGS_DIR, "msa_{ref}_{window}.log"),
    shell:
        """
        mafft --maxiterate 1000 --localpair --thread {threads} {input.fa} > {output.fa} 2> {log}
        """
