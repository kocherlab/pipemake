import pytest
import tempfile

from pipemake.wildcardIO import *

@pytest.mark.parametrize(
    "wildcard_str",
    [
        'tests/files/wildcardIO/{sample}_{read}.fq',
    ],
)
def test_WildcardIO_fromStr_w_error (wildcard_str):
    with pytest.raises(Exception):
        WildcardIO.fromStr(wildcard_str)

@pytest.mark.parametrize(
    "wildcard_str",
    [
        'tests/files/wildcardIO/{sample}_{read}.fq.gz'
    ],
)
def test_WildcardIO_fromStr (wildcard_str):
    test_wildcard = WildcardIO.fromStr(wildcard_str)
    
    assert test_wildcard.wildcard_str == wildcard_str
    assert 'sample' in test_wildcard.wildcard_dict
    assert 'read' in test_wildcard.wildcard_dict

    sample_set = set(test_wildcard.wildcard_dict['sample'])
    assert sample_set == set(['test1', 'test2']) and len(sample_set) == 2

    read_set = set(test_wildcard.wildcard_dict['read'])
    assert read_set == set(['R1', 'R2']) and len(read_set) == 2

@pytest.mark.parametrize(
    "wildcard_str",
    [
        'tests/files/wildcardIO/{sample}_{read}.fq.gz'
    ],
)
def test_WildcardIO_samples_w_error (wildcard_str):
    with pytest.raises(Exception):
        test_wildcard = WildcardIO.fromStr(wildcard_str)
        assert test_wildcard.samples == ['test1', 'test2']

@pytest.mark.parametrize(
    "wildcard_str, sample_wildcard",
    [
        ('tests/files/wildcardIO/{sample}_{read}.fq.gz', 'sample')
    ],
)
def test_WildcardIO_samples_wo_error (wildcard_str, sample_wildcard):
    test_wildcard = WildcardIO.fromStr(wildcard_str, sample_wildcard = sample_wildcard)
    sample_set = set(test_wildcard.samples)
    assert sample_set == set(['test1', 'test2']) and len(sample_set) == 2

@pytest.mark.parametrize(
    "wildcard_str, sample_wildcard",
    [
        ('tests/files/wildcardIO/{sample}_{read}.fq.gz', 'sample')
    ],
)
def test_WildcardIO_samples_wo_error (wildcard_str, sample_wildcard):
    test_wildcard = WildcardIO.fromStr(wildcard_str, sample_wildcard = sample_wildcard)
    sample_set = set(test_wildcard.samples)
    assert sample_set == set(['test1', 'test2']) and len(sample_set) == 2

@pytest.mark.parametrize(
    "wildcard_str, sample_wildcard, copy_method",
    [
        ('tests/files/wildcardIO/{sample}_{read}.fq.gz', 'sample', 'symbolic_link'),
        ('tests/files/wildcardIO/{sample}_{read}.fq.gz', 'sample', 'copy')
    ],
)
def test_WildcardIO_standardizedFiles (wildcard_str, sample_wildcard, copy_method):
    test_dir = tempfile.mkdtemp()
    test_wildcard = WildcardIO.fromStr(wildcard_str, sample_wildcard = sample_wildcard)
    test_wildcard.standardizedFiles('{sample}_{read}.test.fq.gz', out_dir = test_dir, copy_method = copy_method)

    for sample, read in list(itertools.product(['test1', 'test2'], ['R1', 'R2'])):
        standardized_filename = f'{sample}_{read}.test.fq.gz'
        if copy_method == 'symbolic_link': assert os.path.islink(os.path.join(test_dir, standardized_filename))
        elif copy_method == 'copy': assert os.path.isfile(os.path.join(test_dir, standardized_filename))
        else: raise Exception(f'Unsupported copy method: {copy_method}')

@pytest.mark.parametrize(
    "wildcard_str, sample_wildcard, copy_method",
    [
        ('tests/files/wildcardIO/{sample}_{read}.fq.gz', 'sample', 'symbolic_link'),
        ('tests/files/wildcardIO/{sample}_{read}.fq.gz', 'sample', 'copy')
    ],
)
def test_WildcardIO_returnPaths (wildcard_str, sample_wildcard, copy_method):
    test_wildcard = WildcardIO.fromStr(wildcard_str, sample_wildcard = sample_wildcard)
    if copy_method == 'copy': assert test_wildcard.returnPaths(copy_method) == []
    else: assert test_wildcard.returnPaths(copy_method) == [os.path.abspath(os.path.dirname(wildcard_str))]