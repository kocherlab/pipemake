rule all:
    input:
        expand("FASTQ/Processed_Wildcard/{sample}_R1.fq.gz", sample=config["samples"]),


rule fastq_Unprocessed_Wildcards:
    input:
        config["fastq_file_input"],
    output:
        temp(
            "FASTQ/Unprocessed_Wildcards/"
            + os.path.basename(config["fastq_file_input"]).replace(".fastq", ".fq")
        ),
    run:
        import os

        os.symlink(os.path.abspath(input[0]), output[0])


ruleorder: fastq_process_paired_end > fastq_process_single_end > fastq_rename_basic_paired_end > fastq_rename_reads_paired_end > fastq_rename_basic_single_end > fastq_rename_reads_single_end


rule fastq_rename_basic_single_end:
    input:
        r1_reads="FASTQ/Unprocessed_Wildcards/{sample}_1.fq.gz",
    output:
        r1_reads=temp("FASTQ/Unprocessed_Wildcards/{sample}_R1.fq.gz"),
    shell:
        "mv {input.r1_reads} {output.r1_reads}"


rule fastq_rename_basic_paired_end:
    input:
        r1_reads="FASTQ/Unprocessed_Wildcards/{sample}_1.fq.gz",
        r2_reads="FASTQ/Unprocessed_Wildcards/{sample}_2.fq.gz",
    output:
        r1_reads=temp("FASTQ/Unprocessed_Wildcards/{sample}_R1.fq.gz"),
        r2_reads=temp("FASTQ/Unprocessed_Wildcards/{sample}_R2.fq.gz"),
    shell:
        """
        mv {input.r1_reads} {output.r1_reads}
        mv {input.r2_reads} {output.r2_reads}
        """


rule fastq_rename_reads_single_end:
    input:
        r1_reads="FASTQ/Unprocessed_Wildcards/{sample}_read-1.fq.gz",
    output:
        r1_reads=temp("FASTQ/Unprocessed_Wildcards/{sample}_R1.fq.gz"),
    shell:
        "mv {input.r1_reads} {output.r1_reads}"


rule fastq_rename_reads_paired_end:
    input:
        r1_reads="FASTQ/Unprocessed_Wildcards/{sample}_read-1.fq.gz",
        r2_reads="FASTQ/Unprocessed_Wildcards/{sample}_read-4.fq.gz",
    output:
        r1_reads=temp("FASTQ/Unprocessed_Wildcards/{sample}_R1.fq.gz"),
        r2_reads=temp("FASTQ/Unprocessed_Wildcards/{sample}_R2.fq.gz"),
    shell:
        """
        mv {input.r1_reads} {output.r1_reads}
        mv {input.r2_reads} {output.r2_reads}
        """


rule fastq_process_single_end:
    input:
        r1_reads="FASTQ/Unprocessed_Wildcards/{sample}_R1.fq.gz",
    output:
        r1_reads=temp("FASTQ/Processed_Wildcard/{sample}_R1.fq.gz"),
    shell:
        "mv {input.r1_reads} {output.r1_reads}"


rule fastq_process_paired_end:
    input:
        r1_reads="FASTQ/Unprocessed_Wildcards/{sample}_R1.fq.gz",
        r2_reads="FASTQ/Unprocessed_Wildcards/{sample}_R2.fq.gz",
    output:
        r1_reads=temp("FASTQ/Processed_Wildcard/{sample}_R1.fq.gz"),
        r2_reads=temp("FASTQ/Processed_Wildcard/{sample}_R2.fq.gz"),
    shell:
        """
        mv {input.r1_reads} {output.r1_reads}
        mv {input.r2_reads} {output.r2_reads}
        """
