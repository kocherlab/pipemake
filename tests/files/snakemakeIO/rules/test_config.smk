rule script_call:
    input:
        r1_reads=os.path.join(
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    output:
        r1_reads=os.path.join(
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    params:
        test_param_static=config["test1"],
        test_param_path=config["params"]["test2"],
        test_param_optional_1=f"1" if "test3" in config else "0",
        test_param_optional_2=f'--test {config["test4"]}' if "test4" in config else "",
        test_param_optional_3=(
            f"--test {config['test5']}" if "test5" in config["params"] else ""
        ),
        test_param_optional_3=f"1" if "test6" in config["params"] else "0",
        config_test=config["test7"],
        config_test2=config["params"]["test8"],
    resources:
        mem_mb=16000,
    threads: 4
    script:
        "../scripts/test.py"
