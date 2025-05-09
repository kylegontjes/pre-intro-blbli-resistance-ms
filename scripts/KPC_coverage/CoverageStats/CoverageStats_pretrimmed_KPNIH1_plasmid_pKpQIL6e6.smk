# Author: Kyle Gontjes

configfile: "config/config_pretrimmed_KPNIH1_plasmid_pKpQIL6e6.yaml"

import pandas as pd
import os
import re

PREFIX = config["prefix"]

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

if not os.path.exists("results/" + PREFIX):
    try:
        os.makedirs("results/" + PREFIX)
    except OSError as e:
        print(f"Error creating directory: {e}")

# Rule all information
rule all:
    input:
        reference_size_file = expand(str(config["reference_genome"]) +".size"),
        reference_window_file = expand(str(config["reference_genome"]) +".bed"), 
        aligned_sam_file = expand("results/{prefix}/{sample}/align_reads/{sample}_aln.sam",prefix=PREFIX,sample=SAMPLE),
        sorted_bam_file = expand("results/{prefix}/{sample}/post_align/{sample}_sorted_aln.bam",prefix=PREFIX,sample=SAMPLE),
        final_bam = expand("results/{prefix}/{sample}/post_align/remove_duplicates/{sample}_final.bam",prefix=PREFIX,sample=SAMPLE),
        summary = expand("results/{prefix}/{sample}/gatk/{sample}_summary",prefix=PREFIX,sample=SAMPLE)
 
# Step 0.1: Curate reference file for gatk
rule curate_reference_size:
    input:
        reference_genome = config["reference_genome"]
    output:
        reference_size_file = expand(str(config["reference_genome"]) +".size")
    singularity:
        "docker://quay.io/biocontainers/bioawk:1.0--he4a0461_12"
    shell:
        "bioawk -c fastx '{{ print $name, length($seq) }}' < {input.reference_genome} > {output.reference_size_file}"

# Step 0.2: Curate reference bed file for gatk
rule curate_reference_bed:
    input: 
       reference_size_file = expand(str(config["reference_genome"]) +".size")
    output: 
        reference_window_file = expand(str(config["reference_genome"]) +".bed")
    singularity:
        "docker://staphb/bedtools"
    shell:
        "bedtools makewindows -g {input.reference_size_file} -w 1000 > {output.reference_window_file}"
 
# Step 1: Align trimmed reads using bwa
rule align_reads:
    input:
        r1_paired = lambda wildcards: expand(str(config["input_reads"] +"/" f"{wildcards.sample}_R1_trim_paired.fastq.gz")),
        r2_paired = lambda wildcards: expand(str(config["input_reads"] +"/" f"{wildcards.sample}_R2_trim_paired.fastq.gz"))
    output:
        aligned_sam_out = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_aln.sam"
    params:
        reference_genome = config["reference_genome"],
        num_cores = config["num_cores"]
    log:
        bwa_log = "logs/{prefix}/{sample}/bwa/{sample}_bwa.log"
    singularity:
        "docker://staphb/bwa:0.7.17"
    shell:
        "bash ./bash_scripts/bwa.sh {input.r1_paired} {input.r2_paired} {params.reference_genome} {params.num_cores} {output.aligned_sam_out} &> {log.bwa_log}"

# Step 2: Convert sam to bam files
rule sam_to_bam:
    input:
        aligned_sam_file = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_aln.sam"
    output:
        # Bam output files
        bam_file = f"results/{{prefix}}/{{sample}}/post_align/{{sample}}_aln.bam",
        # Bam sorted files
        sorted_bam_file = f"results/{{prefix}}/{{sample}}/post_align/{{sample}}_sorted_aln.bam",
    params:
        outdir_temp = "results/{prefix}/{sample}/post_align/{sample}_sorted_aln_temp"
    log:
        sam_to_bam_log= "logs/{prefix}/{sample}/post_align/{sample}_sam_to_bam.log"
    singularity:
        "docker://staphb/samtools:1.19"
    shell:
        """
        samtools view -Sb {input.aligned_sam_file} > {output.bam_file}  
        samtools sort {output.bam_file} -m 500M -@ 0 -o {output.sorted_bam_file} -T {params.outdir_temp} &> {log.sam_to_bam_log}
        """
        
# Step 3: Remove duplicates
rule remove_pcr_duplicates:
    input:
        sorted_bam = f"results/{{prefix}}/{{sample}}/post_align/{{sample}}_sorted_aln.bam"
    output:
        # Bam dupliactes removed
        bam_duplicates_removed = f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_aln_marked.bam",
        # Picard duplicates_removed_name
        picard_metrics = f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_picard_metrics.txt"
    log:
        picard_log = "logs/{prefix}/{sample}/post_align/remove_duplicates/{sample}_picard.log"
    singularity:
        "docker://broadinstitute/picard:latest"
    shell:
        "java -jar /usr/picard/picard.jar MarkDuplicates -REMOVE_DUPLICATES true -INPUT {input.sorted_bam} -OUTPUT {output.bam_duplicates_removed} -METRICS_FILE {output.picard_metrics} -CREATE_INDEX false -VALIDATION_STRINGENCY LENIENT &> {log.picard_log}"

# Step 4: Sort bam file w/ removed duplicates
rule bam_sort:
    input:
        bam_duplicates_removed = f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_aln_marked.bam"
    output: 
        # duplicates removed sorted bam
        final_bam = f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_final.bam"
    params:
        samtools_temp = "results/{prefix}/{sample}/post_align/remove_duplicates/{sample}_final_temp"
    log:
        post_picard_sort = "logs/{prefix}/{sample}/post_align/remove_duplicates/{sample}_post_picard_sort.log"
    singularity:
        "docker://staphb/samtools:1.19"
    shell:
        """
        samtools sort {input.bam_duplicates_removed} -m 500M -@ 0 -o {output.final_bam} -T {params.samtools_temp} &> {log.post_picard_sort}
        samtools index {output.final_bam}
        """

# Step 5: Calculate coverage statistics
rule gatk_coverage:
    input:
        final_bam = f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_final.bam",
        reference_window_file = expand(str(config["reference_genome"]) +".bed")
    output:
        summary = f"results/{{prefix}}/{{sample}}/gatk/{{sample}}_summary"
    params: 
        reference_genome = config["reference_genome"]
    log:
        gatk_log = "logs/{prefix}/{sample}/gatk/{sample}_gatk.log"
    singularity:
        "docker://broadinstitute/gatk" 
    shell:
        "gatk DepthOfCoverage -R {params.reference_genome} -O {output.summary} -I {input.final_bam}  --summary-coverage-threshold 1 --summary-coverage-threshold 5 --summary-coverage-threshold 9 --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25 --ignore-deletion-sites --intervals {input.reference_window_file} &> {log.gatk_log}"