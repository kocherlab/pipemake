rule all:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa.masked",


ruleorder: joint_softmask_assemblies_pipemake > store_repeatmasker_assembly > softmask_windowmasker_assembly


rule joint_softmask_assemblies_pipemake:
    input:
        repeatmodeler_masked=f"Assembly/RepeatModeler/MaskedAssembly/{config['species']}_{config['assembly_version']}.fa.masked",
        windowmasker_masked=f"Assembly/WindowMasker/{config['species']}_{config['assembly_version']}.fa.masked",
    output:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa.masked",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "joint-softmask --masked-fastas {input.repeatmodeler_masked} {input.windowmasker_masked} --output-fasta {output}"
