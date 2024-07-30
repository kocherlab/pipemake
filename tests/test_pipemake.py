import os
import pytest
import tempfile
import unittest

from unittest.mock import patch

from pipemake.pipemake import *

@pytest.mark.parametrize(
    "wildcard_str",
    [
        'tests/files/wildcardIO/{sample}_{read}.fq.gz',
    ],
)
def test_pipemake_main_wildcard_wo_error (wildcard_str):
    test_dir = tempfile.mkdtemp()
    os.environ['PM_SNAKEMAKE_DIR'] = 'tests/files/snakemakeIO'

    # Assign the workflor prefix
    workflow_prefix = os.path.join(test_dir, 'test_workflow')

    # Assign the command line arguments
    test_cmd = ['', 'fastq-filter', '--fastq-wildcard', wildcard_str, '--workflow-prefix', workflow_prefix, '--resource-yml']

    # Run pipemake with the test command
    with unittest.mock.patch('sys.argv', test_cmd):
        main()

    # Check if the workflow files were created
    assert os.path.isfile(f'{workflow_prefix}.smk')
    assert os.path.isfile(f'{workflow_prefix}.yml')
    assert os.path.isdir(f'{workflow_prefix}')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/pipeline.log')
    assert os.path.isdir(f'{workflow_prefix}/pipemake')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/backups/test_workflow.smk.bkp')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/backups/test_workflow.yml.bkp')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/backups/test_workflow.resources.yml.bkp')
    assert os.path.isdir(f'{workflow_prefix}/pipemake/backups')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/modules/fastq_filter_fastp.smk')
    assert os.path.isdir(f'{workflow_prefix}/pipemake/modules')
    assert os.path.isdir(f'{workflow_prefix}/FASTQ/Unfiltered')

@pytest.mark.parametrize(
    "table_str",
    [
        'tests/files/seqIO/test_table.tsv',
    ],
)
def test_pipemake_main_table_wo_error (table_str):
    test_dir = tempfile.mkdtemp()
    os.environ['PM_SNAKEMAKE_DIR'] = 'tests/files/snakemakeIO'

    # Assign the workflor prefix
    workflow_prefix = os.path.join(test_dir, 'test_workflow')

    # Assign the command line arguments
    test_cmd = ['', 'fastq-filter', '--fastq-table', table_str, '--workflow-prefix', workflow_prefix, '--resource-yml']

    # Run pipemake with the test command
    with unittest.mock.patch('sys.argv', test_cmd):
        main()

    # Check if the workflow files were created
    assert os.path.isfile(f'{workflow_prefix}.smk')
    assert os.path.isfile(f'{workflow_prefix}.yml')
    assert os.path.isdir(f'{workflow_prefix}')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/pipeline.log')
    assert os.path.isdir(f'{workflow_prefix}/pipemake')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/backups/test_workflow.smk.bkp')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/backups/test_workflow.yml.bkp')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/backups/test_workflow.resources.yml.bkp')
    assert os.path.isdir(f'{workflow_prefix}/pipemake/backups')
    assert os.path.isfile(f'{workflow_prefix}/pipemake/modules/fastq_filter_fastp.smk')
    assert os.path.isdir(f'{workflow_prefix}/pipemake/modules')
    assert os.path.isdir(f'{workflow_prefix}/FASTQ/Unfiltered')

@pytest.mark.parametrize(
    "table_str",
    [
        'tests/files/seqIO/test_table2.tsv',
    ],
)
def test_pipemake_main_table_w_error (table_str):
    test_dir = tempfile.mkdtemp()
    os.environ['PM_SNAKEMAKE_DIR'] = 'tests/files/snakemakeIO'

    # Assign the workflor prefix
    workflow_prefix = os.path.join(test_dir, 'test_workflow')

    # Assign the command line arguments
    test_cmd = ['', 'fastq-filter', '--fastq-table', table_str, '--workflow-prefix', workflow_prefix, '--resource-yml']

    # Run pipemake with the test command
    with unittest.mock.patch('sys.argv', test_cmd):
        with pytest.raises(Exception):
            main()


def test_pipemake_main_table_w_error ():
    
    os.environ['PM_SNAKEMAKE_DIR'] = 'tests/files/snakemakeIO'

    # Assign the command line arguments
    test_cmd1 = ['']
    test_cmd2 = ['', 'fastq-filter']

    # Run pipemake with the test command
    with unittest.mock.patch('sys.argv', test_cmd1):
        with pytest.raises(SystemExit) as excinfo:
            main()
            excinfo.value.code == 2
        
    # Run pipemake with the test command
    with unittest.mock.patch('sys.argv', test_cmd2):
        with pytest.raises(SystemExit) as excinfo:
            main()
            excinfo.value.code == 2