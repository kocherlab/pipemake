import os
import pytest
import tempfile
import itertools

from pipemake.processIO import ProcessIO, processInput


@pytest.mark.parametrize(
    "wildcard_str",
    [
        "tests/files/wildcardIO/{samples}_{reads}.fq",
    ],
)
def test_processIO_fromWildcardStr_w_error(wildcard_str):
    with pytest.raises(Exception):
        ProcessIO.fromWildcardStr(wildcard_str)


@pytest.mark.parametrize(
    "wildcard_str",
    ["tests/files/wildcardIO/{samples}_{reads}.fq.gz"],
)
def test_processIO_fromWildcardStr_wo_error(wildcard_str):
    ProcessIO.fromWildcardStr(wildcard_str)


@pytest.mark.parametrize(
    "wildcard_str, standardized_wildcard",
    [
        (
            "tests/files/wildcardIO/{samples}_{reads}.fq.gz",
            "{samples}_{reads}.test.fq.gz",
        )
    ],
)
def test_processIO_fromWildcardStr_standardizeInput(
    wildcard_str, standardized_wildcard
):
    test_dir = tempfile.mkdtemp()
    test_input = processInput(
        method="wildcard_str",
        args={
            "wildcard_str": wildcard_str,
            "standardized_filename": standardized_wildcard,
            "out_dir": test_dir,
        },
    )
    test_input.standardize()

    for sample, read in list(itertools.product(["test1", "test2"], ["R1", "R2"])):
        standardized_filename = f"{sample}_{read}.test.fq.gz"
        assert os.path.isfile(os.path.join(test_dir, standardized_filename))


@pytest.mark.parametrize(
    "wildcard_str",
    [("tests/files/wildcardIO/{samples}_{reads}.fq.gz")],
)
def test_processIO_fromWildcardStr(wildcard_str):
    test_input = processInput(
        method="wildcard_str",
        args={"wildcard_str": wildcard_str, "sample_keywords": ["samples"]},
    )

    assert test_input.returnPaths() == [os.path.abspath(os.path.dirname(wildcard_str))]

    samples_dict = test_input.returnSamples()
    assert "samples" in samples_dict
    assert set(samples_dict["samples"]) == set(["test1", "test2"])


@pytest.mark.parametrize(
    "file_str",
    ["tests/files/fileIO/test2_R1.fq", "tests/files/fileIO/test2_R2.fq"],
)
def test_processIO_fromFileStr_w_error(file_str):
    with pytest.raises(IOError):
        ProcessIO.fromFileStr(file_str)


@pytest.mark.parametrize(
    "file_str",
    ["tests/files/fileIO/test1_R1.fq", "tests/files/fileIO/test1_R2.fq"],
)
def test_processIO_fromFileStr_wo_error(file_str):
    ProcessIO.fromFileStr(file_str)


@pytest.mark.parametrize(
    "file_str, standardized_filename",
    [
        ("tests/files/fileIO/test1_R1.fq", "test1_R1.test.fq.gz"),
        ("tests/files/fileIO/test1_R2.fq", "test1_R2.test.fq.gz"),
    ],
)
def test_processIO_fromFileStr_standardizeInput_wo_error(
    file_str, standardized_filename
):
    test_dir = tempfile.mkdtemp()
    test_input = processInput(
        method="file_str",
        args={
            "file_str": file_str,
            "standardized_filename": standardized_filename,
            "out_dir": test_dir,
        },
    )
    test_input.standardize()
    assert os.path.isfile(os.path.join(test_dir, standardized_filename))


@pytest.mark.parametrize(
    "file_str, standardized_filename",
    [
        ("tests/files/fileIO/test1_R1.fq", "test1_R1.test.fq.gz"),
        ("tests/files/fileIO/test1_R2.fq", "test1_R2.test.fq.gz"),
    ],
)
def test_processIO_fromFileStr_standardizeInput_w_error(
    file_str, standardized_filename
):
    test_input = processInput(
        method="file_str",
        args={
            "file_str": file_str,
            "standardized_filename": standardized_filename,
        },
    )
    assert test_input.returnPaths() == [os.path.abspath(os.path.dirname(file_str))]
    with pytest.raises(Exception):
        test_input.returnSamples()


@pytest.mark.parametrize(
    "table_file",
    ["tests/files/fileIO/test_table2.tsv"],
)
def test_processIO_fromTableFile_w_error(table_file):
    with pytest.raises(IOError):
        ProcessIO.fromTableFile(table_file)


@pytest.mark.parametrize(
    "table_file",
    ["tests/files/fileIO/test_table.tsv"],
)
def test_processIO_fromTableFile_wo_error(table_file):
    ProcessIO.fromTableFile(table_file, sample_keywords=["samples"])


@pytest.mark.parametrize(
    "table_file, standardized_wildcard",
    [
        ("tests/files/fileIO/test_table.tsv", "{samples}_{reads}.test.fq.gz"),
    ],
)
def test_processIO_fromTableFile_standardizeInput(table_file, standardized_wildcard):
    test_dir = tempfile.mkdtemp()
    test_input = processInput(
        method="table_file",
        args={
            "table_file": table_file,
            "standardized_filename": standardized_wildcard,
            "out_dir": test_dir,
            "sample_keywords": ["samples"],
        },
    )

    test_input.standardize()

    for sample, read in list(itertools.product(["test1", "test2"], ["R1", "R2"])):
        standardized_filename = f"{sample}_{read}.test.fq.gz"
        assert os.path.isfile(os.path.join(test_dir, standardized_filename))


@pytest.mark.parametrize(
    "table_file",
    [
        ("tests/files/fileIO/test_table.tsv"),
    ],
)
def test_processIO_fromTableFile(table_file):
    test_input = processInput(
        method="table_file",
        args={"table_file": table_file, "sample_keywords": ["samples"]},
    )

    assert test_input.returnPaths() == [os.path.abspath(os.path.dirname(table_file))]

    samples_dict = test_input.returnSamples()
    assert "samples" in samples_dict
    assert set(samples_dict["samples"]) == set(["test1", "test2"])
    assert "reads" in samples_dict
    assert set(samples_dict["reads"]) == set(["R1", "R2"])
