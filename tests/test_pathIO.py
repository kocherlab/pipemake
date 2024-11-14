import pytest
import os
import tempfile

from pipemake.pathIO import DirIO


@pytest.mark.parametrize(
    "path_name",
    ["tests/files/pathIO/test_dir"],
)
def test_DirIO_exists_wo_error(path_name):
    DirIO.create(path_name)


@pytest.mark.parametrize(
    "path_name",
    ["tests/files/pathIO/test_dir2"],
)
def test_DirIO_exists_w_error(path_name):
    with pytest.raises(IOError):
        DirIO.create(path_name)


@pytest.mark.parametrize(
    "path_name",
    ["tests/files/pathIO/test_dir"],
)
def test_DirIO_standardize_wo_error(path_name):
    test_dir = tempfile.mkdtemp()
    test_Dir = DirIO.create(path_name)
    test_Dir.standardize("new_dir", workflow_prefix=test_dir)
    assert os.path.isdir(os.path.join(test_dir, "new_dir"))
    assert os.path.isfile(os.path.join(test_dir, "new_dir", "fasta_table.tsv"))

    with open(os.path.join(test_dir, "new_dir", "fasta_table.tsv")) as table_file:
        assert table_file.readline() == "genomes\tfilename\n"
        assert (
            table_file.readline() == "genome1\ttests/files/pathIO/test_dir/genome1.fa\n"
        )
        assert (
            table_file.readline() == "genome2\ttests/files/pathIO/test_dir/genome2.fa"
        )


@pytest.mark.parametrize(
    "path_name",
    ["tests/files/pathIO/test_dir"],
)
def test_DirIO_standardize_copy_wo_error(path_name):
    test_dir = tempfile.mkdtemp()
    test_Dir = DirIO.create(path_name)
    test_Dir.standardize("new_dir", workflow_prefix=test_dir, copy_method="copy")
    assert os.path.isdir(os.path.join(test_dir, "new_dir"))
    assert os.path.isfile(os.path.join(test_dir, "new_dir", "fasta_table.tsv"))

    with open(os.path.join(test_dir, "new_dir", "fasta_table.tsv")) as table_file:
        assert table_file.readline() == "genomes\tfilename\n"
        assert (
            table_file.readline() == "genome1\ttests/files/pathIO/test_dir/genome1.fa\n"
        )
        assert (
            table_file.readline() == "genome2\ttests/files/pathIO/test_dir/genome2.fa"
        )


@pytest.mark.parametrize(
    "path_name",
    ["tests/files/pathIO/test_dir"],
)
def test_DirIO_returnPaths_link_wo_error(path_name):
    test_Dir = DirIO.create(path_name)
    assert test_Dir.returnPaths() == [os.path.abspath(path_name)]


@pytest.mark.parametrize(
    "path_name",
    ["tests/files/pathIO/test_dir"],
)
def test_DirIO_returnPaths_copy_wo_error(path_name):
    test_Dir = DirIO.create(path_name)
    assert test_Dir.returnPaths(copy_method="copy") == []
