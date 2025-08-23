rule all:
    input:
        "Indices/sortmerna/.idx.chk",


rule sortmerna_index:
    output:
        index_chk=temp("Indices/sortmerna/.idx.chk"),
        work_dir=temp(directory("Indices/.sortmerna_work_dir")),
    params:
        index_dir="Indices/sortmerna",
        sortmerna_db=config["sortmerna_db"],
    singularity:
        "docker://aewebb/sortmerna:v4.3.7"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "sortmerna -index 1 --ref /opt/DBs/{params.sortmerna_db} --idx-dir {params.index_dir} --workdir {output.work_dir} && touch {output.index_chk}"
