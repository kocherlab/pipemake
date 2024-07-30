import pytest
import tempfile
import yaml
import filecmp

from pipemake.snakemakeIO import *

@pytest.mark.parametrize(
    'job_prefix, pipeline_storage_dir, resource_yml, scale_threads, scale_mem, indent_style, overwrite',
    [
        ('test', 'tests/files/snakemakeIO2', None, 0, 0, '\t', True),
        ('test', 'tests/files/snakemakeIO', 'None', 0, 0, '\t', True),
        ('test', 'tests/files/snakemakeIO', None, '0', 0, '\t', True),
        ('test', 'tests/files/snakemakeIO', None, 0, '0', '\t', True),
        ('test', 'tests/files/snakemakeIO', None, 0, '0', '\t', True),
        ('test', 'tests/files/snakemakeIO', None, 0, 0, ',', True),
        ('test', 'tests/files/snakemakeIO', None, 0, 0, '\t', None),
    ]
)
def test_SnakePipelineIO_w_error (job_prefix, pipeline_storage_dir, resource_yml, scale_threads, scale_mem, indent_style, overwrite):
    test_dir = tempfile.mkdtemp()
    singularity_dir = tempfile.mkdtemp()

    snakemake_job_prefix = os.path.join(test_dir, job_prefix)

    # Test if the function raises an error
    with pytest.raises(Exception):
        SnakePipelineIO.open(snakemake_job_prefix = snakemake_job_prefix, pipeline_storage_dir  = pipeline_storage_dir, pipeline_job_dir = test_dir, resource_yml = resource_yml, scale_threads = scale_threads, scale_mem = scale_mem, singularity_dir = singularity_dir, indent_style = indent_style, overwrite = overwrite)

@pytest.mark.parametrize(
    'job_prefix, pipeline_storage_dir, resource_yml, scale_threads, scale_mem, indent_style, overwrite',
    [
        ('test', 'tests/files/snakemakeIO', None, 1.0, 1.0, '\t', True)
    ]
)
def test_SnakePipelineIO_wo_error (job_prefix, pipeline_storage_dir, resource_yml, scale_threads, scale_mem, indent_style, overwrite):
    test_dir = tempfile.mkdtemp()
    singularity_dir = tempfile.mkdtemp()

    snakemake_job_prefix = os.path.join(test_dir, job_prefix)
    
    # Create a test pipeline
    test_pipeline = SnakePipelineIO.open(snakemake_job_prefix = snakemake_job_prefix, workflow_prefix = snakemake_job_prefix, pipeline_storage_dir  = pipeline_storage_dir, pipeline_job_dir = test_dir, resource_yml = resource_yml, scale_threads = scale_threads, scale_mem = scale_mem, singularity_dir = singularity_dir, indent_style = indent_style, overwrite = overwrite)

    # Add a module to the pipeline
    test_pipeline.addModule('fastq_filter_fastp.smk')

    # Create the output directories
    unfiltered_fastq_dir = os.path.join(test_dir, 'fastq')
    filtered_fastq_dir = os.path.join(test_dir, 'fastq_filtered')

    # Write the pipeline config
    test_pipeline.writeConfig({'min_length': 100, 'samples': ['sample1', 'sample2'], 'threads': 4, 'mem': 16, 'unfiltered_fastq_dir' : unfiltered_fastq_dir, 'filtered_fastq_dir' : filtered_fastq_dir, 'workflow_prefix': snakemake_job_prefix})

    # Write the pipeline
    test_pipeline.writePipeline()

    # Build the singularity containers
    test_pipeline.buildSingularityContainers()

    test_pipeline.close()

    # Check if the pipeline was written
    assert os.path.exists(os.path.join(test_dir, 'test.smk'))

    test_module_path = os.path.join(test_pipeline._module_job_dir, 'fastq_filter_fastp.smk')

    # Check if a string is within the file
    with open(os.path.join(test_dir, 'test.smk'), 'r') as f:
        test_pipeline_content = f.read()
        assert "expand(os.path.join(config['paths']['workflow_prefix'], config['paths']['filtered_fastq_dir'], '{sample}.json'), sample=config['samples'])" in test_pipeline_content
        assert f'include: "{test_module_path}"' in test_pipeline_content

    # Check if the config file was written
    assert os.path.exists(os.path.join(test_dir, 'test.yml'))

    # Check if the singularity containers were built
    assert os.path.isfile(os.path.join(singularity_dir, 'pipemake_utils.v0.1.27.sif'))

    # Read the config yaml
    with open(os.path.join(test_dir, 'test.yml'), 'r') as yaml_file:
        config_dict = yaml.safe_load(yaml_file)
        assert 'samples' in config_dict
        assert ['sample1', 'sample2'] == config_dict['samples']
        assert 'min_length' in config_dict
        assert 100 == config_dict['min_length']
        assert 'fastp_single_end' in config_dict['resources']
        assert 'mem_mb' in config_dict['resources']['fastp_single_end']
        assert 16000 == config_dict['resources']['fastp_single_end']['mem_mb']
        assert 'threads' in config_dict['resources']['fastp_single_end']
        assert 4 == config_dict['resources']['fastp_single_end']['threads']
        assert 'fastp_pair_end' in config_dict['resources']
        assert 'mem_mb' in config_dict['resources']['fastp_pair_end']
        assert 16000 == config_dict['resources']['fastp_pair_end']['mem_mb']
        assert 'threads' in config_dict['resources']['fastp_pair_end']
        assert 4 == config_dict['resources']['fastp_pair_end']['threads']
        assert 'unfiltered_fastq_dir' in config_dict['paths']
        assert unfiltered_fastq_dir == config_dict['paths']['unfiltered_fastq_dir']
        assert 'filtered_fastq_dir' in config_dict['paths']
        assert filtered_fastq_dir == config_dict['paths']['filtered_fastq_dir']

@pytest.mark.parametrize(
    'smk_filename',
    [
        ('tests/files/snakemakeIO/modules/fastq_filter_fastp2.smk')
    ]
)
def test_SnakeFileIO_w_error (smk_filename):
    # Test if the function raises an error
    with pytest.raises(Exception):
        SnakeFileIO.open(smk_filename = smk_filename)

@pytest.mark.parametrize(
    'smk_filename',
    [
        ('tests/files/snakemakeIO/modules/fastq_filter_fastp.smk')
    ]
)
def test_SnakeFileIO_wo_error (smk_filename):

    # Create a test directory
    test_dir = tempfile.mkdtemp()

    # Create a test module
    test_module = SnakeFileIO.open(smk_filename = smk_filename, singularity_dir = test_dir)

    # Create a out filename
    out_filename = os.path.join(test_dir, 'test_fastq_filter_fastp.smk')
    cmp_filename = os.path.join(test_dir, 'cmp_fastq_filter_fastp.smk')
    
    # Generate the module file and confirm it exists
    test_module.write(out_filename)
    assert os.path.exists(out_filename)

    with open(cmp_filename, 'w') as cmp_file:
        with open('tests/files/snakemakeIO/test_fastq_filter_fastp.smk', 'r') as test_file:
            for test_line in test_file:
                if 'URL' in test_line: 
                    cmp_file.write(test_line.format(URL = test_dir))
                else: 
                    cmp_file.write(test_line)
                    
    # Check if the module was created correctly
    assert filecmp.cmp(out_filename, cmp_filename)

@pytest.mark.parametrize(
    'rule_filename',
    [
        ('tests/files/snakemakeIO/rules/fastp_pair_end.smk')
    ]
)
def test_SnakeRuleIO_wo_error (rule_filename):

    # Create a test directory
    test_dir = tempfile.mkdtemp()

    # Create a test rule
    rule_str = ''

    # Read the rule and store it as a string
    with open(rule_filename, 'r') as rule_file:
        for rule_line in rule_file: rule_str += rule_line

    # Test if the function raises an error
    test_rule = SnakeRuleIO.read(rule_str = rule_str, singularity_dir = test_dir, indent_style = '\t')

    # Create a test rule
    test_str = ''

    # Read the rule and store it as a string
    with open('tests/files/snakemakeIO/test_fastp_pair_end.smk', 'r') as test_file:
        for test_line in test_file:
            if 'URL' in test_line: test_str += test_line.format(URL = test_dir)
            else: test_str += test_line
        test_str += '\n'

    print(test_rule._rule_text)

    print(test_str)

    # Check if the rule was created correctly
    assert test_rule._rule_text == test_str
    assert test_rule.rule_name == 'fastp_pair_end'
    assert test_rule._rule_resource_params == {'threads': 4, 'mem_mb': 16000}
    assert test_rule._rule_config_params == {('paths', 'workflow_prefix'), ('paths', 'unfiltered_fastq_dir'), ('paths', 'filtered_fastq_dir'), ('min_length',)}

@pytest.mark.parametrize(
    'rule_name, attribute_type, attribute_text',
    [
        ('fastp_single_end', 'params', "min_length=config['min_length']"),
        ('fastp_pair_end', 'resources', '\t\tmem_mb=16000')
    ]
)
def test_SnakeAttributeIO (rule_name, attribute_type, attribute_text):

    # Assign the attribute
    test_attribute = SnakeAttributeIO.process(rule_name = rule_name, attribute_type = attribute_type, attribute_text = attribute_text, indent_style = '\t')
    
    # Assign the attribute argument
    attribute_text_arg = attribute_text.split('=')[0].strip()

    # Check if the attributes were created correctly
    if attribute_type == 'params':
        with pytest.raises(Exception): test_attribute.updateSnakeResource()
    elif attribute_type == 'resources':
        assert test_attribute.updateSnakeResource() == f"\t\tmem_mb=config['{attribute_type}']['{rule_name}']['{attribute_text_arg}']"