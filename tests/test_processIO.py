import os
import pytest
import tempfile
import itertools

from pipemake.processIO import *

@pytest.mark.parametrize(
    "wildcard_str",
    [
        'tests/files/wildcardIO/{sample}_{read}.fq',
    ],
)
def test_processIO_fromWildcardStr_w_error (wildcard_str):
    with pytest.raises(Exception):
        ProcessIO.fromWildcardStr(wildcard_str)

@pytest.mark.parametrize(
    "wildcard_str",
    [
        'tests/files/wildcardIO/{sample}_{read}.fq.gz'
    ],
)
def test_processIO_fromWildcardStr_wo_error (wildcard_str):
    ProcessIO.fromWildcardStr(wildcard_str)

@pytest.mark.parametrize(
    "wildcard_str, standardized_wildcard",
    [
        ('tests/files/wildcardIO/{sample}_{read}.fq.gz', '{sample}_{read}.test.fq.gz')
    ],
)
def test_processIO_fromWildcardStr_standardizeInput (wildcard_str, standardized_wildcard):
    test_dir = tempfile.mkdtemp()
    standardizeInput(method = 'wildcard-str', args = {'wildcard_str': wildcard_str, 'standardized_filename': standardized_wildcard, 'out_dir': test_dir})

    for sample, read in list(itertools.product(['test1', 'test2'], ['R1', 'R2'])):
        standardized_filename = f'{sample}_{read}.test.fq.gz'
        assert os.path.isfile(os.path.join(test_dir, standardized_filename))

@pytest.mark.parametrize(
    "wildcard_str",
    [
        ('tests/files/wildcardIO/{sample}_{read}.fq.gz')
    ],
)
def test_processIO_fromWildcardStr (wildcard_str):
    assert returnPaths(method = 'wildcard-str', args = {'wildcard_str': wildcard_str}) == [os.path.abspath(os.path.dirname(wildcard_str))]
    assert set(returnSamples(method = 'wildcard-str', args = {'wildcard_str': wildcard_str, 'sample_wildcard': 'sample'})) == set(['test1', 'test2'])

@pytest.mark.parametrize(
    "file_str",
    [
        "tests/files/seqIO/test2_R1.fq",
        "tests/files/seqIO/test2_R2.fq"
    ],
)
def test_processIO_fromFileStr_w_error (file_str):
    with pytest.raises(IOError):
        ProcessIO.fromFileStr(file_str)

@pytest.mark.parametrize(
    "file_str",
    [
        "tests/files/seqIO/test1_R1.fq",
        "tests/files/seqIO/test1_R2.fq"
    ],
)
def test_processIO_fromFileStr_wo_error (file_str):
    ProcessIO.fromFileStr(file_str)

@pytest.mark.parametrize(
    "file_str, standardized_filename",
    [
        ("tests/files/seqIO/test1_R1.fq", "test1_R1.test.fq.gz"),
        ("tests/files/seqIO/test1_R2.fq", "test1_R2.test.fq.gz")
    ],
)
def test_processIO_fromFileStr_standardizeInput (file_str, standardized_filename):
    test_dir = tempfile.mkdtemp()
    standardizeInput(method = 'file-str', args = {'input_filename': file_str, 'standardized_filename': standardized_filename, 'out_dir': test_dir})
    assert os.path.isfile(os.path.join(test_dir, standardized_filename))

@pytest.mark.parametrize(
    "file_str, standardized_filename",
    [
        ("tests/files/seqIO/test1_R1.fq", "test1_R1.test.fq.gz"),
        ("tests/files/seqIO/test1_R2.fq", "test1_R2.test.fq.gz")
    ],
)
def test_processIO_fromFileStr_standardizeInput (file_str, standardized_filename):
    assert returnPaths(method = 'file-str', args = {'input_filename': file_str, 'standardized_filename': standardized_filename}) == [os.path.abspath(os.path.dirname(file_str))]
    with pytest.raises(Exception):
        returnSamples(method = 'file-str', args = {'input_filename': file_str})

@pytest.mark.parametrize(
    "table_file",
    [
        "tests/files/seqIO/test_table2.tsv"
    ],
)
def test_processIO_fromTableFile_w_error (table_file):
    with pytest.raises(IOError):
        ProcessIO.fromTableFile(table_file)

@pytest.mark.parametrize(
    "table_file",
    [
        "tests/files/seqIO/test_table.tsv"
    ],
)
def test_processIO_fromTableFile_wo_error (table_file):
    ProcessIO.fromTableFile(table_file)

@pytest.mark.parametrize(
    "table_file, standardized_wildcard",
    [
        ("tests/files/seqIO/test_table.tsv", "{sample}_{read}.test.fq.gz"),
    ],
)
def test_processIO_fromTableFile_standardizeInput (table_file, standardized_wildcard):
    test_dir = tempfile.mkdtemp()
    standardizeInput(method = 'table-file', args = {'table_filename': table_file, 'standardized_filename': standardized_wildcard, 'out_dir': test_dir})

    for sample, read in list(itertools.product(['test1', 'test2'], ['R1', 'R2'])):
        standardized_filename = f'{sample}_{read}.test.fq.gz'
        assert os.path.isfile(os.path.join(test_dir, standardized_filename))

@pytest.mark.parametrize(
    "table_file",
    [
        ("tests/files/seqIO/test_table.tsv"),
    ],
)
def test_processIO_fromTableFile (table_file):
    assert returnPaths(method = 'table-file', args = {'table_filename': table_file}) == [os.path.dirname(table_file)]
    assert set(returnSamples(method = 'table-file', args = {'table_filename': table_file, 'sample_wildcard': 'sample'})) == set(['test1', 'test2'])
