rule all:
    input:
        expand(
            f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}_converted.fasta",
            ec_type=["ec", "noec"],
        ),
        expand(
            f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}.hic",
            ec_type=["ec", "noec"],
        ),


rule index_hifi_assembly:
    input:
        f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa",
    output:
        f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa.fai",
        f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa.bwt.2bit.64",
    singularity:
        "docker://aewebb/bwa-mem2:v2.2.1"
    resources:
        mem_mb=32000,
    threads: 1
    shell:
        """
        bwa-mem2 index {input}
        samtools faidx {input}
        """


rule align_hic_reads_bwa:
    input:
        read_fastq=f"HiC/FASTQ/{config['species']}_HiC_{{read}}.fq.gz",
        assembly_fasta=f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa",
        index_fasta=f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa.fai",
    output:
        temp(f"HiC/BAM/Aligned/{config['species']}_{{read}}.aligned.bam"),
    singularity:
        "docker://aewebb/bwa-mem2:v2.2.1"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "bwa-mem2 mem -t {threads} {input.assembly_fasta} {input.read_fastq} | samtools view --threads {threads} -bh -o {output}"


rule filter_hic_reads:
    input:
        f"HiC/BAM/Aligned/{config['species']}_{{read}}.aligned.bam",
    output:
        temp(f"HiC/BAM/Aligned/{config['species']}_{{read}}.filtered.bam"),
    singularity:
        "docker://aewebb/arima_mapping:05222024"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "samtools view -h {input} | filter_five_end.pl | samtools view --threads {threads} -bh -o {output}"


rule combine_hic_reads:
    input:
        r1_bam=f"HiC/BAM/Aligned/{config['species']}_R1.filtered.bam",
        r2_bam=f"HiC/BAM/Aligned/{config['species']}_R2.filtered.bam",
        index_fasta=f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa",
    output:
        temp(f"HiC/BAM/Sorted/{config['species']}.sorted.bam"),
    params:
        mapq_filter="10",
    singularity:
        "docker://aewebb/arima_mapping:05222024"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "two_read_bam_combiner.pl {input.r1_bam} {input.r2_bam} samtools {params.mapq_filter} | samtools view -bh -t {input.index_fasta} | samtools sort -@ {threads} -o {output}"


rule add_read_groups:
    input:
        f"HiC/BAM/Sorted/{config['species']}.sorted.bam",
    output:
        temp(f"HiC/BAM/Sorted/{config['species']}.sorted_groups.bam"),
    params:
        species=config["species"],
    singularity:
        "docker://aewebb/gatk4:v4.6.1.0"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gatk AddOrReplaceReadGroups INPUT={input} OUTPUT={output} ID={params.species} LB={params.species} SM={params.species} PL=ILLUMINA PU=none"


rule mark_duplicates:
    input:
        f"HiC/BAM/Sorted/{config['species']}.sorted_groups.bam",
    output:
        bam=f"HiC/BAM/Sorted/{config['species']}.dedup.bam",
        metrics=f"HiC/BAM/Sorted/{config['species']}.dedup.metrics.txt",
    params:
        tmp_dir=".tmp",
    singularity:
        "docker://aewebb/gatk4:v4.6.1.0"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        """
        gatk MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} TMP_DIR={params.tmp_dir} ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
        samtools index {output.bam}
        """


rule dedup_bam_stats:
    input:
        f"HiC/BAM/Sorted/{config['species']}.dedup.bam",
    output:
        f"HiC/BAM/Sorted/{config['species']}.dedup.stats.txt",
    singularity:
        "docker://aewebb/arima_mapping:05222024"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "get_stats.pl {input} > {output}"


rule yahs_ec:
    input:
        dedup_bam=f"HiC/BAM/Sorted/{config['species']}.dedup.bam",
        index_fasta=f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa",
    output:
        temp(
            f"Assembly/YaHS/{config['species']}_{config['assembly_version']}_yahsout_ec.bin"
        ),
        temp(
            f"Assembly/YaHS/{config['species']}_{config['assembly_version']}_yahsout_ec_scaffolds_final.agp"
        ),
    params:
        out_prefix=f"Assembly/YaHS/{config['species']}_{config['assembly_version']}_yahsout_ec",
    singularity:
        "docker://aewebb/yahs:v1.2.2"
    resources:
        mem_mb=32000,
    threads: 1
    shell:
        "yahs {input.index_fasta} {input.dedup_bam} -o {params.out_prefix}"


rule yahs_noec:
    input:
        dedup_bam=f"HiC/BAM/Sorted/{config['species']}.dedup.bam",
        index_fasta=f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa",
    output:
        temp(
            f"Assembly/YaHS/{config['species']}_{config['assembly_version']}_yahsout_noec.bin"
        ),
        temp(
            f"Assembly/YaHS/{config['species']}_{config['assembly_version']}_yahsout_noec_scaffolds_final.agp"
        ),
    params:
        out_prefix=f"Assembly/YaHS/{config['species']}_{config['assembly_version']}_yahsout_noec",
    singularity:
        "docker://aewebb/yahs:v1.2.2"
    resources:
        mem_mb=32000,
    threads: 1
    shell:
        "yahs --no-contig-ec {input.index_fasta} {input.dedup_bam} -o {params.out_prefix}"


rule yahs_juicer_pre:
    input:
        reads_bin=f"Assembly/YaHS/{config['species']}_{config['assembly_version']}_yahsout_{{ec_type}}.bin",
        assembly_agp=f"Assembly/YaHS/{config['species']}_{config['assembly_version']}_yahsout_{{ec_type}}_scaffolds_final.agp",
        assembly_index=f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa.fai",
    output:
        txt=temp(
            f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}.txt"
        ),
        agp=temp(
            f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}.liftover.agp"
        ),
        log=f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}_juicer_pre.log",
    params:
        out_prefix=f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}",
    singularity:
        "docker://aewebb/yahs:v1.2.2"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "juicer pre -a -o {params.out_prefix} {input.reads_bin} {input.assembly_agp} {input.assembly_index} 2> {output.log}"


rule juicer_tools_pre:
    input:
        txt=f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}.txt",
        log=f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}_juicer_pre.log",
    output:
        hic=f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}.hic",
    singularity:
        "docker://aewebb/juicer_tools:v1.19.02"
    resources:
        mem_mb=32000,
    threads: 1
    shell:
        """
        assembly_size=$(grep 'PRE_C_SIZE' {input.log} | awk '{{print $3}}')
        java -jar /opt/juicer_tools.jar pre {input.txt} {output.hic} <(echo "assembly ${{assembly_size}}")
        """


rule create_converted_fasta:
    input:
        agp=f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}.liftover.agp",
        fasta=f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa",
    output:
        bed=temp(
            f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}.bed"
        ),
        fasta=f"Assembly/Juicebox/{config['species']}_{config['assembly_version']}_{{ec_type}}_converted.fasta",
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
