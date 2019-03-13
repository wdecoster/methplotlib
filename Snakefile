configfile: "config.yaml"


def get_fast5(wildcards):
    return config["fast5"][wildcards.sample]


def get_fastq(wildcards):
    return config["fastq"][wildcards.sample]


def get_bam(wildcards):
    return config["bam"][wildcards.sample]


def get_summary(wildcards):
    return config["summary"][wildcards.sample]


rule all:
    pass

rule nanopolish_index:
    input:
        f5 = get_fast5,
        fq = get_fastq,
        sm = get_summary,
    output:
        "indices/index_done_{sample}"
    shell:
        """
        nanopolish index -f {input.sm}-d {input.f5}/ {input.fq}
        touch {output}
        """

rule call_methylation:
    input:
        idx = "indices/index_done_{sample}",
        fq = get_fastq,
        bam = get_bam,
        genome = config["genome"],
        region = config["region"],
    output:
        "methylation_calls/{input.region}_{sample}.tsv"
    threads: 8
    shell:
        """
        nanopolish call-methylation \
         -t {threads} \
         -r {input.fq} \
         -b {input.bam} \
         -g {input.genome} \
         -w {input.region} > {output}
        """

rule methylation_frequency:
    input:
        "methylation_calls/{input.region}_{sample}.tsv"
    output:
        "methylation_frequency/{input.region}_{sample}.tsv"
    shell:
        """
        scripts/calculate_methylation_frequency.py -i {input} > {output}
        """
