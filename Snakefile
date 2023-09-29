import os


def get_fast5(wildcards):
    return config["fast5"][wildcards.sample]


def get_fastq(wildcards):
    return config["fastq"][wildcards.sample]


def get_bam(wildcards):
    return config["bam"][wildcards.sample]


def get_region(wildcards):
    return config["region"][wildcards.region]


rule all:
    input:
        expand("{prefix}methylation_frequency/{region}_{sample}.tsv",
               sample=config["fastq"],
               region=config["region"],
               prefix=config.get("prefix", ""))

rule nanopolish_index:
    input:
        f5 = get_fast5,
        fq = get_fastq,
    output:
        "{prefix}indices/index_done_{sample}"
    threads: 10  # Just to ensure that this is not ran in parallel for too many samples
    log:
        "{prefix}logs/nanopolish_index/index_{sample}.log"
    shell:
        """
        nanopolish index -d {input.f5}/ {input.fq} 2> {log}
        touch {output}
        """

rule call_methylation:
    input:
        idx = "indices/index_done_{sample}",
        fq = get_fastq,
        bam = get_bam,
        genome = config["genome"],
    output:
        "{prefix}methylation_calls/{region}_{sample}.tsv"
    threads: 8
    params:
        region = get_region,
    log:
        "{prefix}logs/methylation/methcall_{region}_{sample}.log"
    shell:
        """
        nanopolish call-methylation \
         -t {threads} \
         -r {input.fq} \
         -b {input.bam} \
         -g {input.genome} \
         -w {params.region} > {output} 2> {log}
        """

rule methylation_frequency:
    input:
        "{prefix}methylation_calls/{region}_{sample}.tsv"
    output:
        "{prefix}methylation_frequency/{region}_{sample}.tsv"
    log:
        "{prefix}logs/methylation/methfreq_{region}_{sample}.log"
    params:
        script = os.path.join(workflow.basedir, "scripts/calculate_methylation_frequency.py")
    shell:
        """
        python {params.script} -i {input} > {output} 2> {log}
        """
