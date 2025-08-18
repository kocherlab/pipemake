import os
import pytest
import tempfile
import unittest.mock


from pipemake.pipemake import main


@pytest.mark.parametrize(
    "wildcard_str",
    [
        "tests/files/wildcardIO/{samples}_{reads}.fq.gz",
    ],
)
def test_pipemake_main_wildcard_wo_error(wildcard_str):
    test_dir = tempfile.mkdtemp()
    os.environ["PM_SNAKEMAKE_DIR"] = "tests/files/snakemakeIO"

    # Assign the workflor prefix
    workflow_prefix = os.path.join(test_dir, "test_workflow")

    # Create the test path
    test_prefix = os.path.join(test_dir, "test_dir")
    os.makedirs(test_prefix, exist_ok=True)

    # Assign the command line arguments
    test_cmd = [
        "",
        "fastq-filter",
        "--fastq-wildcard",
        wildcard_str,
        "--workflow-dir",
        workflow_prefix,
        "--resource-yml",
        "--path-test",
        test_prefix,
    ]

    # Run pipemake with the test command
    with unittest.mock.patch("sys.argv", test_cmd):
        main()

    # Check if the workflow files were created
    assert os.path.isfile(f"{workflow_prefix}/Snakefile")
    assert os.path.isfile(f"{workflow_prefix}/config.yml")
    assert os.path.isdir(f"{workflow_prefix}")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/pipeline.log")
    assert os.path.isdir(f"{workflow_prefix}/pipemake")
    assert os.path.isdir(f"{workflow_prefix}/pipemake/backups")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/backups/Snakefile.bkp")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/backups/config.yml.bkp")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/backups/resources.yml.bkp")
    assert os.path.isdir(f"{workflow_prefix}/pipemake/modules")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/modules/fastq_filter_fastp.smk")
    assert not os.path.isfile(f"{workflow_prefix}/pipemake/modules/test_script.smk")
    assert os.path.isdir(f"{workflow_prefix}/pipemake/modules")
    assert os.path.isdir(f"{workflow_prefix}/FASTQ/Unfiltered")


@pytest.mark.parametrize(
    "table_str",
    [
        "tests/files/fileIO/test_table.tsv",
    ],
)
def test_pipemake_main_table_wo_error(table_str):
    test_dir = tempfile.mkdtemp()
    os.environ["PM_SNAKEMAKE_DIR"] = "tests/files/snakemakeIO"

    # Assign the workflow prefix
    workflow_prefix = os.path.join(test_dir, "test_workflow")

    # Create the test path
    test_prefix = os.path.join(test_dir, "test_dir")
    os.makedirs(test_prefix, exist_ok=True)

    # Assign the command line arguments
    test_cmd = [
        "",
        "fastq-filter",
        "--fastq-table",
        table_str,
        "--workflow-dir",
        workflow_prefix,
        "--resource-yml",
        "--dir-test",
        "tests/files/fileIO",
        "--path-test",
        test_prefix,
    ]

    # Run pipemake with the test command
    with unittest.mock.patch("sys.argv", test_cmd):
        main()

    # Check if the workflow files were created
    assert os.path.isfile(f"{workflow_prefix}/Snakefile")
    assert os.path.isfile(f"{workflow_prefix}/config.yml")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/pipeline.log")
    assert os.path.isdir(f"{workflow_prefix}/pipemake")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/backups/Snakefile.bkp")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/backups/config.yml.bkp")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/backups/resources.yml.bkp")
    assert os.path.isdir(f"{workflow_prefix}/pipemake/backups")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/modules/fastq_filter_fastp.smk")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/modules/test_script.smk")
    assert os.path.isdir(f"{workflow_prefix}/pipemake/modules")
    assert os.path.isdir(f"{workflow_prefix}/FASTQ/Unfiltered")

    with open(f"{workflow_prefix}/config.yml") as f:
        for line in f:
            print(line)


@pytest.mark.parametrize(
    "table_str",
    [
        "tests/files/fileIO/test_table2.tsv",
    ],
)
def test_pipemake_main_table_w_error(table_str):
    test_dir = tempfile.mkdtemp()
    os.environ["PM_SNAKEMAKE_DIR"] = "tests/files/snakemakeIO"

    # Assign the workflor prefix
    workflow_prefix = os.path.join(test_dir, "test_workflow")

    # Assign the command line arguments
    test_cmd = [
        "",
        "fastq-filter",
        "--fastq-table",
        table_str,
        "--workflow-prefix",
        workflow_prefix,
        "--resource-yml",
    ]

    # Run pipemake with the test command
    with unittest.mock.patch("sys.argv", test_cmd):
        with pytest.raises(Exception):
            main()


def test_pipemake_main_help_wo_error():
    os.environ["PM_SNAKEMAKE_DIR"] = "tests/files/snakemakeIO"

    # Assign the command line arguments
    test_cmd1 = [""]
    test_cmd2 = ["", "fastq-filter"]

    # Run pipemake with the test command
    with unittest.mock.patch("sys.argv", test_cmd1):
        with pytest.raises(SystemExit) as excinfo:
            main()
            excinfo.value.code == 2

    # Run pipemake with the test command
    with unittest.mock.patch("sys.argv", test_cmd2):
        with pytest.raises(SystemExit) as excinfo:
            main()
            excinfo.value.code == 2


@pytest.mark.parametrize(
    "wildcard_str",
    [
        "tests/files/wildcardIO/{samples}_{reads}.fq.gz",
    ],
)
def test_pipemake_main_links_wo_error(wildcard_str):
    test_dir = tempfile.mkdtemp()
    os.environ["PM_SNAKEMAKE_DIR"] = "tests/files/snakemakeIO_links"

    # Assign the workflor prefix
    workflow_prefix = os.path.join(test_dir, "test_workflow")

    # Create the test path
    test_prefix = os.path.join(test_dir, "test_dir")
    os.makedirs(test_prefix, exist_ok=True)

    # Assign the command line arguments
    test_cmd = [
        "",
        "trim-fastqs",
        "--fastq-wildcard",
        wildcard_str,
        "--workflow-dir",
        workflow_prefix,
    ]

    # Run pipemake with the test command
    with unittest.mock.patch("sys.argv", test_cmd):
        main()

    # Check if the workflow files were created
    assert os.path.isfile(f"{workflow_prefix}/Snakefile")
    assert os.path.isfile(f"{workflow_prefix}/config.yml")
    assert os.path.isdir(f"{workflow_prefix}")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/pipeline.log")
    assert os.path.isdir(f"{workflow_prefix}/pipemake")
    assert os.path.isdir(f"{workflow_prefix}/pipemake/backups")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/backups/Snakefile.bkp")
    assert os.path.isfile(f"{workflow_prefix}/pipemake/backups/config.yml.bkp")
    assert os.path.isdir(f"{workflow_prefix}/pipemake/modules")
    assert os.path.isfile(
        f"{workflow_prefix}/pipemake/modules/fastq_process_wildcards.smk"
    )
    assert os.path.isfile(f"{workflow_prefix}/pipemake/modules/fastq_trim_fastp.smk")
    assert not os.path.isfile(
        f"{workflow_prefix}/pipemake/modules/fastq_sra_paired_end_download.smk"
    )
    assert not os.path.isfile(
        f"{workflow_prefix}/pipemake/modules/fastq_sra_single_end_download.smk"
    )
    assert os.path.isdir(f"{workflow_prefix}/pipemake/modules")
    assert os.path.isdir(f"{workflow_prefix}/FASTQ/Input")
