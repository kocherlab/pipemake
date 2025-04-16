rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["assembly_dir"],
                "{sample}_{ec_type}_converted.fasta",
            ),
            sample=config["samples"],
            ec_type=["ec", "noec"],
        ),

rule index_hifi_assembly_samtools:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa.fai",
        ),
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "samtools faidx {input}"

rule index_hifi_assembly_bwa:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa.amb",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa.ann",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa.bwt",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa.pac",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa.sa",
        ),
    singularity:
        "docker://aewebb/bwa:v0.7.18"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "bwa index {input}"


rule align_hic_reads_bwa:
    input:
        read_fastq=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_fastq_dir"],
            "{sample}_{read}.fq.gz
        ),
        index_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["hic_aligned_bam_dir"],
                "{sample}_{read}.aligned.bam",
            )
        ),
    singularity:
        "docker://aewebb/bwa:v0.7.18"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "bwa mem -t {threads} {input.index_fasta} {input.read_fastq} | samtools view --threads {threads} -bh -o {output}"

rule filter_hic_reads:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_aligned_bam_dir"],
            "{sample}_{read}.aligned.bam",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["hic_aligned_bam_dir"],
                "{sample}_{read}.filtered.bam",
            )
        ),
    singularity:
        "docker://aewebb/arima_mapping:05222024"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "samtools view -h {input} | filter_five_end.pl | samtools view {threads} -bh -o {output}"

rule combine_hic_reads:
    input:
        r1_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_aligned_bam_dir"],
            "{sample}_R1.filtered.bam",
        ),
        r2_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_aligned_bam_dir"],
            "{sample}_R2.filtered.bam",
        ),
        index_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["hic_sorted_bam_dir"],
                "{sample}.sorted.bam",
            )
        ),
    params:
        mapq_filter="10",
    singularity:
        "docker://aewebb/arima_mapping:05222024"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "two_read_bam_combiner.pl {input.r1_bam} {input.r2_bam} samtools {params.mapq_filter} | samtools view -bh -t {input.index_fasta} | samtools sort -@ {threads} -o {output}"

rule add_read_groups:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_sorted_bam_dir"],
            "{sample}.sorted.bam",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["hic_sorted_bam_dir"],
                "{sample}.sorted_groups.bam",
            )
        ),
    params:
        mapq_filter="10",
    singularity:
        "docker://aewebb/gatk4:v4.6.1.0"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gatk AddOrReplaceReadGroups INPUT={input} OUTPUT={output} ID={wildcards.sample} LB={wildcards.sample} SM={wildcards.sample} PL=ILLUMINA PU=none"

rule mark_duplicates:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_sorted_bam_dir"],
            "{sample}.sorted_groups.bam",
        ),
    output: 
        bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_sorted_bam_dir"],
            "{sample}.dedup.bam",
        ),
        metrics=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_sorted_bam_dir"],
            "{sample}.dedup.metrics.txt",
        ),
    params:
        tmp_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            ".tmp",
        ),
    singularity:
        "docker://aewebb/gatk4:v4.6.1.0"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        """
        gatk MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} TMP_DIR={params.tmp_dir} ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
        samtools index {output.bam}
        """

rule dedup_bam_stats:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_sorted_bam_dir"],
            "{sample}.dedup.bam",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_sorted_bam_dir"],
            "{sample}.dedup.stats.txt",
        ),
    singularity:
        "docker://aewebb/arima_mapping:05222024"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "get_stats.pl {input} > {output}"

rule yahs_ec:
    input:
        dedup_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_sorted_bam_dir"],
            "{sample}.dedup.bam",
        ),
        index_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['assembly_sample']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_yahs_dir"],
            "{sample}_yahsout_ec.bin",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_yahs_dir"],
            "{sample}_yahsout_ec_scaffolds_final.agp",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_yahs_dir"],
            "{sample}_yahsout_ec",
        ),
    singularity:
        "docker://aewebb/yahs:v1.2.2"
    resources:
        mem_mb=32000,
    threads: 1
    shell:
        "yahs {input.index_fasta} {input.dedup_bam} -o {params.out_prefix}"


rule yahs_noec:
    input:
        dedup_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_sorted_bam_dir"],
            "{sample}.dedup.bam",
        ),
        index_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_yahs_dir"],
            "{sample}_yahsout_noec.bin",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_yahs_dir"],
            "{sample}_yahsout_noec_scaffolds_final.agp",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_yahs_dir"],
            "{sample}_yahsout_noec",
        ),
    singularity:
        "docker://aewebb/yahs:v1.2.2"
    resources:
        mem_mb=32000,
    threads: 1
    shell:
        "yahs --no-contig-ec {input.index_fasta} {input.dedup_bam} -o {params.out_prefix}"

rule yahs_juicer_pre_noec:
    input:
        reads_bin=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_yahs_dir"],
            "{sample}_yahsout_noec.bin",
        ),
        assembly_agp=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hic_yahs_dir"],
            "{sample}_yahsout_{ec_type}_scaffolds_final.agp",
        ),
        assembly_index=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa.fai",
        ),
    output:
        txt=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'juicebox',
            "{sample}_{ec_type}.txt",
        ),
        log=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'juicebox',
            "{sample}_{ec_type}_juicer_pre.log",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'juicebox',
            "{sample}_{ec_type}",
        ),
    singularity:
        "docker://aewebb/yahs:v1.2.2"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "juicer pre -a -o {params.out_prefix} {input.reads_bin} {input.assembly_agp} {input.assembly_index} 2> {output.log}"

rule juicer_tools_pre:
    input:
        txt=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'juicebox',
            "{sample}_{ec_type}.txt",
        ),
        log=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'juicebox',
            "{sample}_{ec_type}_juicer_pre.log",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'juicebox',
            "{sample}_{ec_type}.hic",
        ),
    singularity:
        "docker://aewebb/juicer_tools:v1.19.02"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        """
        assembly_size=$(grep 'PRE_C_SIZE' {input.log} | awk '{{print $3}}')
        java -jar /opt/juicer_tools.jar pre {input.txt} {output} <(echo "assembly ${assembly_size}")
        """

rule create_liftover_fasta:
    input:
        agp=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'juicebox',
            "{sample}_{ec_type}.liftover.agp",
        ),
        fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'hifiasm',
            f"{config['assembly_sample']}.fa",
        ),
    output:
        bed=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["assembly_dir"],
                'juicebox',
                "{sample}_{ec_type}.bed",
            )
        ),
        fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'juicebox',
            "{sample}_{ec_type}_converted.fasta",
        ),
    singularity:
        "docker://aewebb/bedtools:v2.31.1"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        """
        awk -v OFS='\t' '{{print $6, $7, $8, $1}}' {input.agp} > {output.bed}
        bedtools getfasta -nameOnly -fi {input.fasta} -bed {output.bed} > {output.fasta}
        """

rule copy_liftover_fasta:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            'juicebox',
            "{sample}_{ec_type}_converted.fasta",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            "{sample}_{ec_type}_converted.fasta",
        ),
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "cp {input} {output}"