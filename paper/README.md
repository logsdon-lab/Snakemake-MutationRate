
Requires the following files:
1. `cfg_template`
    * Config template for mutation rate workflow.
    * `/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/config/config_template.yaml`
    ```yaml
    mem: 100GB
    threads: 8
    output_dir: "results"
    log_dir: "logs"
    benchmark_dir: "benchmarks"
    reference:
    - name: ""
        path: ""
        bed: []
        window_size: 10000
    regex_sm_hap: "'(mat|pat|haplotype1|haplotype2|hap1|hap2|h1|h2)'"
    regex_sm_chrom: "'_(hsa(.*?)|chr(.*?)|cen(.*?))[:_]'"
    regex_ref_chrom: "(?:chr|hsa|cen)([0-9XY]{1,2})"
    length_range: [8000, 12000]
    display_order: [SAS, AFR, AMR, EAS, EUR]
    # https://www.nature.com/articles/s41586-025-08816-3
    samples:
    - name: ""
        path: ""
        mm2_opts: -x asm20 -K 8G --eqx -Y -t 8 -r 500000,500000
        divergence_time: 0
    ```
2. `ref_bed`
    * All reference bed coordinates.
    * `/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/data/all_ref.bed`
    ```
    #chrom	chromStart	chromEnd	contigNameCoords	divergenceTimeEstimate	sample	chromName	cladeNumber	contigNameCoords
    HG01352_chr10_haplotype2-0000163	37045901	43675849	0.016	HG01352	chr10	10	HG01352_chr10_haplotype2-0000163:38545900-42175849
    HG01352_chr10_haplotype1-0000016	37095519	44252088	0.068	HG01352	chr10	11	HG01352_chr10_haplotype1-0000016:38595518-42752088
    ```
3. `qry_bed`
    * All query bed coordinates.
    * `/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/data/all_query.bed`
    ```
    #chrom	chromStart	chromEnd	contigNameCoords	cladeName
    NA19650_chr10_haplotype2-0000081	37113889	44294588	NA19650_chr10_haplotype2-0000081:38613889-42794588	chr10@10
    HG00358_chr10_haplotype1-0000023	36420134	45162695	HG00358_chr10_haplotype1-0000023:37920134-43662695	chr10@10
    ```
4. `fasta_dir`
    * Fasta dir with both reference and query assemblies.
    * `/project/logsdon_shared/projects/HGSVC3/new_65_asms_renamed/`
5. `ref_data_sources`
    * Reference annotation data sources.
    * JSON file map with keys: []
    * `/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/config/ref_sources.json`
    ```json
    {
        "hor_arr": "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/all_array_length.bed",
        "rm": "/project/logsdon_shared/projects/HGSVC3/Snakemake-Repeat-Annotation/results/repeatmasker/repeats/{sm}/{chrom_coords}.fa.out",
        "hor": "/project/logsdon_shared/projects/HGSVC3/Snakemake-Repeat-Annotation/results/humas_annot/{sm}_{chrom_coords}/stv_row.bed",
        "cdr": "/project/logsdon_shared/projects/HGSVC3/HGSVC_centromere_annotation/all_cdrs_revised.bed",
        "mdp": "/project/logsdon_shared/projects/HGSVC3/Snakemake-Repeat-Annotation/results/moddotplot/{sm}/{chrom_coords}/{chrom_coords}.bed"
    }
    ```
6. `color_code`
    * Color codes for samples.
    * `/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/config/sample_colors.tsv`
    ```
    HG01890	#E49C1C
    HG02011	#E49C1C
    HG02282	#E49C1C
    ```
7. `plot_template`
    * Cenplot template.
    * `/project/logsdon_shared/projects/HGSVC3/Snakemake-MutationRate/paper/config/base.toml`
    ```toml
    [settings]
    format = ["png", "pdf"]
    transparent = true
    dim = [18.0, 8.0]
    axis_h_pad = 0.1
    dpi = 600

    [[tracks]]
    title = "CDR"
    position = "relative"
    type = "label"
    proportion = 0.01
    path = "cdr"
    options = { legend = false, hide_x = true }
    ```

Generate config.
```bash
# Requires seqtk in PATH
python paper/generate_config.py
```

Run mutation rate workflow per clade.
```bash
module load singularity
snakemake -np -s paper/run_multiple_mut.smk -c 12 --config input_dir=paper/input/cfg
```
