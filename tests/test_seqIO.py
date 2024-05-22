import pytest
import tempfile

from pipemake.seqIO import *

@pytest.mark.parametrize(
    "filename, expected",
    [
        ("tests/files/seqIO/bgzip.gz", True),
        ("tests/files/seqIO/gzip.gz", True),
        ("tests/files/seqIO/text.txt", False)
    ],
)
def test_checkIfGzipped (filename, expected):
    assert checkIfGzipped(filename) == expected

@pytest.mark.parametrize(
    "filename",
    [
        "tests/files/seqIO/test1_R1.fq",
        "tests/files/seqIO/test1_R2.fq"
    ],
)
def test_SeqFileIO_exists_wo_error (filename):
    try:
        SeqFileIO.create(filename)
    except IOError:
        pytest.fail("SeqFileIO.create() raised an error when opening a valid file.")

@pytest.mark.parametrize(
    "filename",
    [
        "tests/files/seqIO/test2_R1.fq",
        "tests/files/seqIO/test2_R2.fq"
    ],
)
def test_SeqFileIO_exists_w_error (filename):
    with pytest.raises(IOError):
        SeqFileIO.create(filename)


@pytest.mark.parametrize(
    "filename, is_gzipped, standardized_filename, gzipped, copy_method",
    [
        ("tests/files/seqIO/test1_R1.fq", False, "test1_R1.fq", True, "symbolic_link"),
        ("tests/files/seqIO/test1_R1.fq.gz", True, "test1_R1.fq", True, "symbolic_link"),
        ("tests/files/seqIO/test1_R1.fq", False,"test1_R1.fq", False, "symbolic_link"),
        ("tests/files/seqIO/test1_R1.fq.gz", True, "test1_R1.fq", False, "symbolic_link"),
        ("tests/files/seqIO/test1_R1.fq", False,"test1_R1.fq", None, "symbolic_link"),
        ("tests/files/seqIO/test1_R1.fq.gz", True, "test1_R1.fq", None, "symbolic_link"),
        ("tests/files/seqIO/test1_R1.fq.gz", True, "test1_R1.fq", True, "copy"),
        ("tests/files/seqIO/test1_R1.fq", False, "test1_R1.fq", True, "copy"),
        ("tests/files/seqIO/test1_R1.fq", False, "test1_R1.fq", False, "copy")
    ],
)
def test_SeqFileIO_standardize (filename, is_gzipped, standardized_filename, gzipped, copy_method):
    test_dir = tempfile.mkdtemp()
    try:
        test_seqfile = SeqFileIO.create(filename)
        test_seqfile.standardize(standardized_filename, out_dir = test_dir, gzipped = gzipped, copy_method = copy_method)

        # Test the copy method
        if copy_method == 'copy': 
            assert os.path.isfile(os.path.join(test_dir, standardized_filename))
            return

        # Test the symbolic link methods
        if is_gzipped == True and gzipped == True:
            assert os.path.islink(os.path.join(test_dir, standardized_filename))
        
        elif is_gzipped == False and gzipped == True:
            assert not os.path.islink(os.path.join(test_dir, standardized_filename))
            assert os.path.isfile(os.path.join(test_dir, standardized_filename))

        elif is_gzipped == False and gzipped == False:
            assert os.path.islink(os.path.join(test_dir, standardized_filename))
        
        elif is_gzipped == True and gzipped == False:
            assert not os.path.islink(os.path.join(test_dir, standardized_filename))
            assert os.path.isfile(os.path.join(test_dir, standardized_filename))

        elif is_gzipped == True and gzipped == None:
            assert os.path.islink(os.path.join(test_dir, standardized_filename))

        elif is_gzipped == False and gzipped == None:
            assert os.path.islink(os.path.join(test_dir, standardized_filename))

    except IOError:
        pytest.fail("SeqFileIO.standardize() raised an error.")


@pytest.mark.parametrize(
    "filename, copy_method",
    [
        ("tests/files/seqIO/test1_R1.fq", "symbolic_link"),
        ("tests/files/seqIO/test1_R1.fq", "copy")
    ],
)
def test_SeqFileIO_returnPaths (filename, copy_method):
    try:
        test_seqfile = SeqFileIO.create(filename)
        if copy_method == 'copy': assert test_seqfile.returnPaths(copy_method) == []
        else: assert test_seqfile.returnPaths(copy_method) == [os.path.dirname(test_seqfile.filename)]
        
    except IOError:
        pytest.fail("SeqFileIO.returnPaths() raised an error.")


@pytest.mark.parametrize(
    "filename, standardized_filename, gzipped, copy_method",
    [
        ("tests/files/seqIO/test1_R1.fq.gz", "test1_R1.fq", True, "symbolic_link"),
        ("tests/files/seqIO/test1_R1.fq.gz", "test1_R1.fq", True, "copy")
    ],
)
def test_SeqFileIO_args (filename, standardized_filename, gzipped, copy_method):
    test_dir = tempfile.mkdtemp()
    try:
        test_seqfile = SeqFileIO.create(filename)
        test_seqfile.standardize(standardized_filename, out_dir = test_dir, gzipped = gzipped, copy_method = copy_method)
        if copy_method == 'copy': assert test_seqfile.args == {}
        else: assert test_seqfile.args == {'bind':os.path.dirname(test_seqfile.filename)}
        
    except IOError:
        pytest.fail("SeqFileIO.returnPaths() raised an error.")


@pytest.mark.parametrize(
    "filename",
    [
        "tests/files/seqIO/test_table.tsv"
    ],
)
def test_SeqFileIO_args (filename, standardized_filename, gzipped, copy_method):
    test_dir = tempfile.mkdtemp()
    try:
        test_seqfile = SeqFileIO.create(filename)
        test_seqfile.standardize(standardized_filename, out_dir = test_dir, gzipped = gzipped, copy_method = copy_method)
        if copy_method == 'copy': assert test_seqfile.args == {}
        else: assert test_seqfile.args == {'bind':os.path.dirname(test_seqfile.filename)}
        
    except IOError:
        pytest.fail("SeqFileIO.returnPaths() raised an error.")


        fromFilenameStr





