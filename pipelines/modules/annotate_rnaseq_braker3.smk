rule all:
    input:
        "Annotations/BRAKER3/braker.gff3",
        "Annotations/BRAKER3/braker.codingseq",
        "Annotations/BRAKER3/braker.aa",


rule annotate_braker3:
    input:
        masked_assembly=f"Assembly/{config['species']}_{config['assembly_version']}.fa.masked",
        merged_bam=f"RNAseq/BAM/{config['species']}_{config['assembly_version']}.bam",
        protein_hints="Homology/ProteinHints.fa",
        augustus_check="Downloads/augustus/.config.chk",
    output:
        "Annotations/BRAKER3/braker.gff3",
        "Annotations/BRAKER3/braker.codingseq",
        "Annotations/BRAKER3/braker.aa",
    params:
        annotations_dir="Annotations/BRAKER3/",
        augustus_config=os.path.abspath("Downloads/augustus/config"),
    singularity:
        "docker://teambraker/braker3:v3.0.7.6"
    resources:
        mem_mb=32000,
    threads: 20
    shell:
        "braker.pl --genome {input.masked_assembly} --prot_seq {input.protein_hints} --bam {input.merged_bam} -gff3 --softmasking --threads {threads} --workingdir {params.annotations_dir} --AUGUSTUS_CONFIG_PATH {params.augustus_config}"
