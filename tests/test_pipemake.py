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
def test_pipemake_main (wildcard_str):
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