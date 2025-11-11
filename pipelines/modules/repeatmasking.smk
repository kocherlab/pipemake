rule all:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa.masked",


rule repeat_modeler:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        f"Assembly/RepeatModeler/Families/{config['species']}_{config['assembly_version']}-families.fa",
    params:
        repeatmodeler_db=f"Assembly/RepeatModeler/DB/{config['species']}_{config['assembly_version']}_DB",
        db_dir="Assembly/RepeatModeler/DB",
        wd_dir="Assembly/RepeatModeler/WorkingDirectory",
    resources:
        mem_mb=40000,
    threads: 20
    singularity:
        "docker://dfam/tetools:1.94"
    shell:
        r"""
        mkdir -p {params.db_dir}
        mkdir -p {params.wd_dir}
        BuildDatabase -name {params.repeatmodeler_db} {input}
        RepeatModeler -database {params.repeatmodeler_db} -threads {threads} -LTRStruct
        rm_working_dir=$(grep 'Using output directory' {params.repeatmodeler_db}-rmod.log | grep -o '[^/]*$')
        tar -czf {params.wd_dir}/RepeatModeler_WD.tar.gz $rm_working_dir
        rm -rf $rm_working_dir
        mv {params.repeatmodeler_db}-families.fa {output}
        mv {params.repeatmodeler_db}-families.stk {params.wd_dir}
        mv {params.repeatmodeler_db}-rmod.log {params.wd_dir}
        """


rule repeat_masker:
    input:
        assembly=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
        families=f"Assembly/RepeatModeler/Families/{config['species']}_{config['assembly_version']}-families.fa",
    output:
        os.path.join(
            f"Assembly/RepeatModeler/MaskedAssembly/{config['species']}_{config['assembly_version']}.fa.masked",
        ),
    params:
        mask_dir="Assembly/RepeatModeler/MaskedAssembly",
    resources:
        mem_mb=24000,
    threads: 12
    singularity:
        "docker://dfam/tetools:1.94"
    shell:
        "RepeatMasker -par {threads} -dir {params.mask_dir} -lib {input.families} {input.assembly}"


rule softmask:
    input:
        assembly=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
        masked_assembly=f"Assembly/RepeatModeler/MaskedAssembly/{config['species']}_{config['assembly_version']}.fa.masked",
    output:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa.masked",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "softmask --input-fasta {input.assembly} --hard-masked-fasta {input.masked_assembly} --output-fasta {output}"
