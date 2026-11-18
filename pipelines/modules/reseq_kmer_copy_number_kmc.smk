ruleorder: kmc_sample_paired_end > kmc_sample_single_end


rule all:
    input:
        f"KmerCopyNumber/{config['species']}_{config['assembly_version']}.kmer_depth.tsv",


rule index_assembly_samtools:
    input:
        config["assembly_input"],
    output:
        f"{config['assembly_input']}.fai",
    log:
        "logs/samtools/kmer_depth.index_assembly.log",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "samtools faidx {input} &> {log}"


rule target_region_samtools:
    input:
        assembly=config["assembly_input"],
        index=f"{config['assembly_input']}.fai",
    output:
        "KmerCopyNumber/Target/target.fa",
    log:
        "logs/samtools/kmer_depth.target_region.log",
    params:
        region=f"{config['chrom']}:{config['start']}-{config['end']}",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "samtools faidx {input.assembly} {params.region} > {output} 2> {log}"


rule kmc_assembly:
    input:
        config["assembly_input"],
    output:
        multiext("KmerCopyNumber/Assembly/assembly", ".kmc_pre", ".kmc_suf"),
    log:
        "logs/kmc/kmer_depth.assembly.log",
    params:
        k=config["k"],
        db_prefix=subpath(output[0], strip_suffix=".kmc_pre"),
        tmp_dir="KmerCopyNumber/tmp/assembly",
        ram_gb=lambda wildcards, resources: resources.mem_mb // 1000,
    singularity:
        "docker://aewebb/kmc:v3.2.4"
    resources:
        mem_mb=16000,
    threads: 8
    shell:
        """
        mkdir -p {params.tmp_dir}
        kmc -k{params.k} -ci1 -fm -t{threads} -m{params.ram_gb} {input} {params.db_prefix} {params.tmp_dir} &> {log}
        rm -rf {params.tmp_dir}
        """


rule kmc_target:
    input:
        "KmerCopyNumber/Target/target.fa",
    output:
        multiext("KmerCopyNumber/Target/target", ".kmc_pre", ".kmc_suf"),
    log:
        "logs/kmc/kmer_depth.target.log",
    params:
        k=config["k"],
        db_prefix=subpath(output[0], strip_suffix=".kmc_pre"),
        tmp_dir="KmerCopyNumber/tmp/target",
        ram_gb=lambda wildcards, resources: resources.mem_mb // 1000,
    singularity:
        "docker://aewebb/kmc:v3.2.4"
    resources:
        mem_mb=4000,
    threads: 2
    shell:
        """
        mkdir -p {params.tmp_dir}
        kmc -k{params.k} -ci1 -fm -t{threads} -m{params.ram_gb} {input} {params.db_prefix} {params.tmp_dir} &> {log}
        rm -rf {params.tmp_dir}
        """


rule kmc_target_assembly_counts:
    input:
        target=multiext("KmerCopyNumber/Target/target", ".kmc_pre", ".kmc_suf"),
        assembly=multiext("KmerCopyNumber/Assembly/assembly", ".kmc_pre", ".kmc_suf"),
    output:
        temp("KmerCopyNumber/Target/target.assembly_counts.txt"),
    log:
        "logs/kmc/kmer_depth.target.assembly_counts.log",
    params:
        target_db=lambda wildcards, input: input.target[0].removesuffix(".kmc_pre"),
        assembly_db=lambda wildcards, input: input.assembly[0].removesuffix(".kmc_pre"),
        tmp_db="KmerCopyNumber/tmp/target.assembly_counts",
    singularity:
        "docker://aewebb/kmc:v3.2.4"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        """
        kmc_tools simple {params.target_db} {params.assembly_db} intersect {params.tmp_db} -ocright &> {log}
        kmc_tools transform {params.tmp_db} dump {output} &>> {log}
        rm -f {params.tmp_db}.kmc_pre {params.tmp_db}.kmc_suf
        """


rule filter_target_kmers_awk:
    input:
        "KmerCopyNumber/Target/target.assembly_counts.txt",
    output:
        "KmerCopyNumber/Query/target.fa",
    params:
        max_count=config["max_assembly_count"],
    localrule: True
    shell:
        "awk -v max={params.max_count} '$2 <= max {{print \">\" NR; print $1}}' {input} > {output}"


rule kmc_query:
    input:
        "KmerCopyNumber/Query/target.fa",
    output:
        multiext("KmerCopyNumber/Query/target", ".kmc_pre", ".kmc_suf"),
    log:
        "logs/kmc/kmer_depth.query.log",
    params:
        k=config["k"],
        db_prefix=subpath(output[0], strip_suffix=".kmc_pre"),
        tmp_dir="KmerCopyNumber/tmp/query",
        ram_gb=lambda wildcards, resources: resources.mem_mb // 1000,
    singularity:
        "docker://aewebb/kmc:v3.2.4"
    resources:
        mem_mb=4000,
    threads: 2
    shell:
        """
        mkdir -p {params.tmp_dir}
        kmc -k{params.k} -ci1 -fm -t{threads} -m{params.ram_gb} {input} {params.db_prefix} {params.tmp_dir} &> {log}
        rm -rf {params.tmp_dir}
        """


rule kmc_sample_single_end:
    input:
        r1_reads="FASTQ/Trimmed/{sample}_R1.fq.gz",
    output:
        multiext("KmerCopyNumber/Samples/{sample}", ".kmc_pre", ".kmc_suf"),
    log:
        "logs/kmc/{sample}.log",
    params:
        k=config["k"],
        min_count=config["min_count"],
        db_prefix=subpath(output[0], strip_suffix=".kmc_pre"),
        tmp_dir="KmerCopyNumber/tmp/{sample}",
        reads_list="KmerCopyNumber/tmp/{sample}.list",
        ram_gb=lambda wildcards, resources: resources.mem_mb // 1000,
    singularity:
        "docker://aewebb/kmc:v3.2.4"
    resources:
        mem_mb=16000,
    threads: 8
    shell:
        """
        mkdir -p {params.tmp_dir}
        printf '%s\\n' {input.r1_reads} > {params.reads_list}
        kmc -k{params.k} -ci{params.min_count} -cs100000 -t{threads} -m{params.ram_gb} -fq @{params.reads_list} {params.db_prefix} {params.tmp_dir} &> {log}
        rm -rf {params.tmp_dir} {params.reads_list}
        """


rule kmc_sample_paired_end:
    input:
        r1_reads="FASTQ/Trimmed/{sample}_R1.fq.gz",
        r2_reads="FASTQ/Trimmed/{sample}_R2.fq.gz",
    output:
        multiext("KmerCopyNumber/Samples/{sample}", ".kmc_pre", ".kmc_suf"),
    log:
        "logs/kmc/{sample}.log",
    params:
        k=config["k"],
        min_count=config["min_count"],
        db_prefix=subpath(output[0], strip_suffix=".kmc_pre"),
        tmp_dir="KmerCopyNumber/tmp/{sample}",
        reads_list="KmerCopyNumber/tmp/{sample}.list",
        ram_gb=lambda wildcards, resources: resources.mem_mb // 1000,
    singularity:
        "docker://aewebb/kmc:v3.2.4"
    resources:
        mem_mb=16000,
    threads: 8
    shell:
        """
        mkdir -p {params.tmp_dir}
        printf '%s\\n' {input.r1_reads} {input.r2_reads} > {params.reads_list}
        kmc -k{params.k} -ci{params.min_count} -cs100000 -t{threads} -m{params.ram_gb} -fq @{params.reads_list} {params.db_prefix} {params.tmp_dir} &> {log}
        rm -rf {params.tmp_dir} {params.reads_list}
        """


rule kmc_sample_histogram:
    input:
        multiext("KmerCopyNumber/Samples/{sample}", ".kmc_pre", ".kmc_suf"),
    output:
        "KmerCopyNumber/Histograms/{sample}.hist.txt",
    log:
        "logs/kmc/{sample}.histogram.log",
    params:
        db_prefix=lambda wildcards, input: input[0].removesuffix(".kmc_pre"),
    singularity:
        "docker://aewebb/kmc:v3.2.4"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "kmc_tools transform {params.db_prefix} histogram {output} -cx100000 &> {log}"


rule kmc_sample_target_counts:
    input:
        sample=multiext("KmerCopyNumber/Samples/{sample}", ".kmc_pre", ".kmc_suf"),
        query=multiext("KmerCopyNumber/Query/target", ".kmc_pre", ".kmc_suf"),
    output:
        "KmerCopyNumber/Counts/{sample}.txt",
    log:
        "logs/kmc/{sample}.target_counts.log",
    params:
        sample_db=lambda wildcards, input: input.sample[0].removesuffix(".kmc_pre"),
        query_db=lambda wildcards, input: input.query[0].removesuffix(".kmc_pre"),
        tmp_db="KmerCopyNumber/tmp/{sample}.target_counts",
    singularity:
        "docker://aewebb/kmc:v3.2.4"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        """
        kmc_tools simple {params.sample_db} {params.query_db} intersect {params.tmp_db} -ocleft &> {log}
        kmc_tools transform {params.tmp_db} dump {output} &>> {log}
        rm -f {params.tmp_db}.kmc_pre {params.tmp_db}.kmc_suf
        """


rule summarize_kmer_depth:
    input:
        counts=expand("KmerCopyNumber/Counts/{sample}.txt", sample=config["samples"]),
        query="KmerCopyNumber/Query/target.fa",
        histograms=expand(
            "KmerCopyNumber/Histograms/{sample}.hist.txt", sample=config["samples"]
        ),
    output:
        f"KmerCopyNumber/{config['species']}_{config['assembly_version']}.kmer_depth.tsv",
    localrule: True
    run:
        import statistics

        with open(input.query) as fh:
            n_query = sum(1 for line in fh if line.startswith(">"))

        def target_median(sample):
            with open(f"KmerCopyNumber/Counts/{sample}.txt") as fh:
                found = [int(line.split()[1]) for line in fh]
            # Query k-mers absent from the reads (below min_count) are counted as 0
            return statistics.median(found + [0] * (n_query - len(found))), len(found)

        def histogram_peak(sample):
            # Single-copy (genome-wide) k-mer coverage: first maximum after the error tail
            with open(f"KmerCopyNumber/Histograms/{sample}.hist.txt") as fh:
                hist = {int(l.split()[0]): int(l.split()[1]) for l in fh if l.strip()}
            ys = [hist.get(x, 0) for x in range(1, max(hist) + 1)] if hist else []
            if not ys:
                return 0
            i = 0
            while i + 1 < len(ys) and ys[i + 1] <= ys[i]:
                i += 1
            return max(range(i, len(ys)), key=lambda j: ys[j]) + 1

        rows = []
        for sample in config["samples"]:
            median, found = target_median(sample)
            peak = histogram_peak(sample)
            rows.append([sample, found, median, peak, median / peak if peak else None])

        ratios = [r[4] for r in rows if r[4] is not None]
        cohort_median = statistics.median(ratios) if ratios else None

        header = [
            "sample",
            "n_target_kmers",
            "n_target_kmers_in_reads",
            "target_median",
            "genome_peak",
            "depth_ratio",
            "rel_to_cohort_median",
        ]
        with open(output[0], "w") as out:
            out.write("\t".join(header) + "\n")
            for sample, found, median, peak, ratio in rows:
                ratio_str = round(ratio, 3) if ratio is not None else "NA"
                rel_str = (
                    round(ratio / cohort_median, 3)
                    if ratio is not None and cohort_median
                    else "NA"
                )
                out.write(
                    "\t".join(
                        map(
                            str,
                            [sample, n_query, found, median, peak, ratio_str, rel_str],
                        )
                    )
                    + "\n"
                )
