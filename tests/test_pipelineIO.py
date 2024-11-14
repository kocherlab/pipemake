import pytest

from pipemake.pipelineIO import ConfigPipelinesIO, ConfigPipelineIO


@pytest.mark.parametrize(
    "config_path",
    ["tests/files/snakemakeIO/"],
)
def test_ConfigPipelinesIO_fromDirectory_wo_error(config_path):
    # Read in the pipeline configs
    configs = ConfigPipelinesIO.fromDirectory(config_path)

    # Check the pipeline configs
    assert isinstance(configs, ConfigPipelinesIO)
    assert len(configs) == 1
    assert "fastq-filter" in configs
    assert isinstance(configs["fastq-filter"], ConfigPipelineIO)


@pytest.mark.parametrize(
    "config_path",
    ["tests/files/snakemakeIO/configs"],
)
def test_ConfigPipelinesIO_fromDirectory_w_error(config_path):
    with pytest.raises(Exception):
        ConfigPipelinesIO.fromDirectory(config_path)


@pytest.mark.parametrize(
    "config_yaml",
    ["tests/files/snakemakeIO/configs/fastq_filter.yml"],
)
def test_ConfigPipelineIO_fromYAML_wo_error(config_yaml):
    # Read in the pipeline configs
    config = ConfigPipelineIO.fromYAML(config_yaml)

    # Check the pipeline configs
    assert isinstance(config, ConfigPipelineIO)
    assert config.name == "fastq-filter"
    assert config.version == 1.0


@pytest.mark.parametrize(
    "config_yaml",
    ["tests/files/snakemakeIO/configs/fastq_filter.yaml"],
)
def test_ConfigPipelineIO_fromYAML_w_error(config_yaml):
    with pytest.raises(Exception):
        ConfigPipelineIO.fromYAML(config_yaml)


@pytest.mark.parametrize(
    "config_yaml",
    ["tests/files/snakemakeIO/configs/fastq_filter.yml"],
)
def test_ConfigPipelineIO_fromYAML_parser_dict(config_yaml):
    # Read in the pipeline configs
    config = ConfigPipelineIO.fromYAML(config_yaml)

    assert "help" in config.parser_dict
    assert "arg-groups" in config.parser_dict
    assert "basic" in config.parser_dict["arg-groups"]
    assert "paths" in config.parser_dict["arg-groups"]
    assert "args" in config.parser_dict["arg-groups"]["basic"]
    assert "args" in config.parser_dict["arg-groups"]["paths"]
    assert "mutually-exclusive-args" in config.parser_dict["arg-groups"]["basic"]
    assert (
        "input-parser"
        in config.parser_dict["arg-groups"]["basic"]["mutually-exclusive-args"]
    )
    assert "fastq-wildcard" in config.parser_dict["arg-groups"]["basic"]["args"]
    assert "type" in config.parser_dict["arg-groups"]["basic"]["args"]["fastq-wildcard"]
    assert "help" in config.parser_dict["arg-groups"]["basic"]["args"]["fastq-wildcard"]
    assert (
        "mutually-exclusive"
        in config.parser_dict["arg-groups"]["basic"]["args"]["fastq-wildcard"]
    )


@pytest.mark.parametrize(
    "config_yaml",
    ["tests/files/snakemakeIO/configs/fastq_filter.yml"],
)
def test_ConfigPipelineIO_fromYAML_setup_dict(config_yaml):
    # Read in the pipeline configs
    config = ConfigPipelineIO.fromYAML(config_yaml)

    assert "fastq_input" in config._setup_dict
    assert "wildcard-method" in config._setup_dict["fastq_input"]
    assert "input" in config._setup_dict["fastq_input"]["wildcard-method"]
    assert "args" in config._setup_dict["fastq_input"]["wildcard-method"]["input"]
    assert set(
        [
            "workflow_prefix",
            "fastq-wildcard",
            "fastq-standardized-wildcard",
            "unfiltered-fastq-dir",
        ]
    ) == set(config._setup_dict["fastq_input"]["wildcard-method"]["input"]["args"])
    assert "standardize" in config._setup_dict["fastq_input"]["wildcard-method"]
    assert (
        "method" in config._setup_dict["fastq_input"]["wildcard-method"]["standardize"]
    )
    assert "args" in config._setup_dict["fastq_input"]["wildcard-method"]["standardize"]
    assert (
        "out_dir"
        in config._setup_dict["fastq_input"]["wildcard-method"]["standardize"]["args"]
    )
    assert (
        "gzipped"
        in config._setup_dict["fastq_input"]["wildcard-method"]["standardize"]["args"]
    )


@pytest.mark.parametrize(
    "config_yaml",
    ["tests/files/snakemakeIO/configs/fastq_filter.yml"],
)
def test_ConfigPipelineIO_fromYAML_snakefiles(config_yaml):
    # Read in the pipeline configs
    config = ConfigPipelineIO.fromYAML(config_yaml)

    assert "fastq_filter_fastp" in config._snakefiles
