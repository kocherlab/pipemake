rule all:
    input:
        expand("FASTQ/Filtered/{sample}.fcsfilt.fastq.gz", sample=config["samples"]),


ruleorder: pacbio_bam_to_fasta > hifi_fastq_to_fastq


rule pacbio_bam_to_fasta:
    input:
        pacbio_bam="BAM/PacBio/{sample}.bam",
    output:
        temp("FASTQ/Unfiltered/{sample}.fasta"),
    singularity:
        "docker://aewebb/bamtools:v2.5.2"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "bamtools convert -format fasta -in {input} -out {output}"


rule hifi_fastq_to_fastq:
    input:
        "FASTQ/Unfiltered/{sample}.fastq.gz",
    output:
        temp("FASTQ/Unfiltered/{sample}.fasta"),
    singularity:
        "docker://aewebb/seqkit:v2.10.0"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "seqkit fq2fa {input} -o {output}"


rule hifi_reads_screen_fcs_adaptor:
    input:
        "FASTQ/Unfiltered/{sample}.fasta",
    output:
        "FASTQ/Unfiltered/fcs-adaptor/{sample}/fcs_adaptor_report.txt",
    params:
        out_prefix="FASTQ/Unfiltered/fcs-adaptor/{sample}",
        prok="--prok" if config["prok"] else "",
        euk="--euk" if config["euk"] else "",
    singularity:
        "docker://ncbi/fcs-adaptor:0.5.4"
    resources:
        mem_mb=32000,
        shell_exec="sh",
    threads: 1
    shell:
        "av_screen_x {input} --output {params.out_prefix} {params.prok} {params.euk}"


rule link_fastq_hifiadapterfiltFCS:
    input:
        "FASTQ/Unfiltered/{sample}.fastq.gz",
    output:
        temp("FASTQ/Filtered/{sample}.fastq.gz"),
    run:
        import os

        os.symlink(os.path.abspath(input[0]), output[0])


rule hifi_screen_hifiadapterfiltFCS:
    input:
        reads="FASTQ/Filtered/{sample}.fastq.gz",
        fcs_adaptor_report="FASTQ/Unfiltered/fcs-adaptor/{sample}/fcs_adaptor_report.txt",
    output:
        "FASTQ/Filtered/{sample}.fcsfilt.fastq.gz",
        temp("FASTQ/Filtered/{sample}.blocklist"),
        temp("FASTQ/Filtered/{sample}.fastq"),
    params:
        out_dir="FASTQ/Filtered/",
    singularity:
        "docker://aewebb/hifiadapterfilt:v3.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "hifiadapterfiltFCS.sh -r {input.reads} -f {input.fcs_adaptor_report} -o {params.out_dir} -t {threads}"
