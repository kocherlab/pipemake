import pytest
import os
import tempfile
import itertools

from pipemake.fileIO import checkIfGzipped, FileIO, TableIO


@pytest.mark.parametrize(
    "filename, expected",
    [
        ("tests/files/fileIO/bgzip.gz", True),
        ("tests/files/fileIO/gzip.gz", True),
        ("tests/files/fileIO/text.txt", False),
    ],
)
def test_checkIfGzipped(filename, expected):
    assert checkIfGzipped(filename) == expected


@pytest.mark.parametrize(
    "filename",
    ["tests/files/fileIO/test1_R1.fq", "tests/files/fileIO/test1_R2.fq"],
)
def test_FileIO_exists_wo_error(filename):
    try:
        FileIO.create(filename)
    except IOError:
        pytest.fail("FileIO.create() raised an error when opening a valid file.")


@pytest.mark.parametrize(
    "filename",
    ["tests/files/fileIO/test2_R1.fq", "tests/files/fileIO/test2_R2.fq"],
)
def test_FileIO_exists_w_error(filename):
    with pytest.raises(IOError):
        FileIO.create(filename)


@pytest.mark.parametrize(
    "filename, is_gzipped, standardized_filename, gzipped, copy_method",
    [
        ("tests/files/fileIO/test1_R1.fq", False, "test1_R1.fq", True, "symbolic_link"),
        (
            "tests/files/fileIO/test1_R1.fq.gz",
            True,
            "test1_R1.fq",
            True,
            "symbolic_link",
        ),
        (
            "tests/files/fileIO/test1_R1.fq",
            False,
            "test1_R1.fq",
            False,
            "symbolic_link",
        ),
        (
            "tests/files/fileIO/test1_R1.fq.gz",
            True,
            "test1_R1.fq",
            False,
            "symbolic_link",
        ),
        ("tests/files/fileIO/test1_R1.fq", False, "test1_R1.fq", None, "symbolic_link"),
        (
            "tests/files/fileIO/test1_R1.fq.gz",
            True,
            "test1_R1.fq",
            None,
            "symbolic_link",
        ),
        ("tests/files/fileIO/test1_R1.fq.gz", True, "test1_R1.fq", True, "copy"),
        ("tests/files/fileIO/test1_R1.fq", False, "test1_R1.fq", True, "copy"),
        ("tests/files/fileIO/test1_R1.fq", False, "test1_R1.fq", False, "copy"),
    ],
)
def test_FileIO_standardize(
    filename, is_gzipped, standardized_filename, gzipped, copy_method
):
    test_dir = tempfile.mkdtemp()
    test_File = FileIO.create(filename)
    test_File.standardize(
        standardized_filename,
        out_dir=test_dir,
        gzipped=gzipped,
        copy_method=copy_method,
    )

    # Test the copy method
    if copy_method == "copy":
        assert os.path.isfile(os.path.join(test_dir, standardized_filename))
        return

    # Test the symbolic link methods
    if is_gzipped is True and gzipped is True:
        assert os.path.islink(os.path.join(test_dir, standardized_filename))

    elif is_gzipped is False and gzipped is True:
        assert not os.path.islink(os.path.join(test_dir, standardized_filename))
        assert os.path.isfile(os.path.join(test_dir, standardized_filename))

    elif is_gzipped is False and gzipped is False:
        assert os.path.islink(os.path.join(test_dir, standardized_filename))

    elif is_gzipped is True and gzipped is False:
        assert not os.path.islink(os.path.join(test_dir, standardized_filename))
        assert os.path.isfile(os.path.join(test_dir, standardized_filename))

    elif is_gzipped is True and gzipped is None:
        assert os.path.islink(os.path.join(test_dir, standardized_filename))

    elif is_gzipped is False and gzipped is None:
        assert os.path.islink(os.path.join(test_dir, standardized_filename))


@pytest.mark.parametrize(
    "filename, copy_method",
    [
        ("tests/files/fileIO/test1_R1.fq", "symbolic_link"),
        ("tests/files/fileIO/test1_R1.fq", "copy"),
    ],
)
def test_FileIO_returnPaths(filename, copy_method):
    test_File = FileIO.create(filename)
    if copy_method == "copy":
        assert test_File.returnPaths(copy_method) == []
    else:
        assert test_File.returnPaths(copy_method) == [
            os.path.dirname(test_File.filename)
        ]


@pytest.mark.parametrize(
    "filename",
    ["tests/files/fileIO/test_table2.tsv"],
)
def test_TableIO_no_table_error(filename):
    with pytest.raises(IOError):
        TableIO.fromFilenameStr(filename)


@pytest.mark.parametrize(
    "filename",
    ["tests/files/fileIO/test_table.tsv"],
)
def test_TableIO_fromFilenameStr(filename):
    test_seqtable = TableIO.fromFilenameStr(filename, sample_column="samples")
    assert "samples" in test_seqtable.samples
    assert set(test_seqtable.samples["samples"]) == set(["test1", "test2"])
    assert "reads" in test_seqtable.samples
    assert set(test_seqtable.samples["reads"]) == set(["R1", "R2"])
    assert test_seqtable._sample_column == "samples"
    assert test_seqtable._file_columns == {"reads"}
    assert test_seqtable._table_columns == {"samples", "reads"}


@pytest.mark.parametrize(
    "filename",
    ["tests/files/fileIO/fasta_table.tsv"],
)
def test_TableIO_unique_column(filename):
    test_dir = tempfile.mkdtemp()
    test_seqtable = TableIO.fromFilenameStr(filename, sample_column="genomes")
    assert "genomes" in test_seqtable.samples
    assert set(test_seqtable.samples["genomes"]) == set(["genome1", "genome2"])
    assert test_seqtable._sample_column == "genomes"
    assert test_seqtable._file_columns == {"filename"}
    assert test_seqtable._table_columns == {"genomes"}

    test_seqtable.standardizedFiles(
        "{genomes}.fa", out_dir=test_dir, copy_method="symbolic_link"
    )
    assert os.path.islink(os.path.join(test_dir, "genome1.fa"))
    assert os.path.islink(os.path.join(test_dir, "genome2.fa"))

    for genome in ["genome1", "genome2"]:
        with open(os.path.join(test_dir, f"{genome}.fa")) as genome_file:
            assert (
                genome_file.read()
                == f">{genome.capitalize()}.1\nATCG\n>{genome.capitalize()}.2\nATCG"
            )


@pytest.mark.parametrize(
    "filename",
    ["tests/files/fileIO/id_table.tsv"],
)
def test_TableIO_no_file_column(filename):
    test_seqtable = TableIO.fromFilenameStr(filename, sample_column="SRA")
    assert "SRA" in test_seqtable.samples
    assert set(test_seqtable.samples["SRA"]) == set(
        ["SRR000001", "SRR000002", "SRR000003", "SRR000004"]
    )
    assert test_seqtable._sample_column == "SRA"
    assert test_seqtable._file_columns == set()
    assert test_seqtable._table_columns == {"SRA"}


@pytest.mark.parametrize(
    "filename, copy_method",
    [
        ("tests/files/fileIO/test_table.tsv", "symbolic_link"),
        ("tests/files/fileIO/test_table.tsv", "copy"),
    ],
)
def test_TableIO_standardizedFiles(filename, copy_method):
    test_dir = tempfile.mkdtemp()
    test_seqtable = TableIO.fromFilenameStr(filename, sample_column="samples")
    test_seqtable.standardizedFiles(
        "{samples}_{reads}.test.fq.gz", out_dir=test_dir, copy_method=copy_method
    )

    for samples, reads in list(itertools.product(["test1", "test2"], ["R1", "R2"])):
        standardized_filename = f"{samples}_{reads}.test.fq.gz"
        if copy_method == "symbolic_link":
            assert os.path.islink(os.path.join(test_dir, standardized_filename))
        elif copy_method == "copy":
            assert os.path.isfile(os.path.join(test_dir, standardized_filename))
        else:
            raise Exception(f"Unsupported copy method: {copy_method}")


@pytest.mark.parametrize(
    "filename, copy_method",
    [
        ("tests/files/fileIO/test_table.tsv", "symbolic_link"),
        ("tests/files/fileIO/test_table.tsv", "copy"),
    ],
)
def test_TableIO_returnPaths(filename, copy_method):
    test_seqtable = TableIO.fromFilenameStr(filename, sample_column="samples")
    if copy_method == "copy":
        assert test_seqtable.returnPaths(copy_method) == []
    else:
        assert test_seqtable.returnPaths(copy_method) == [os.path.dirname(filename)]


@pytest.mark.parametrize(
    "filename",
    ["tests/files/fileIO/test_table.tsv"],
)
def test_TableIO_returnSamples(filename):
    test_seqtable = TableIO.fromFilenameStr(filename, sample_column="samples")
    assert test_seqtable.returnSamples()["samples"] == ["test1", "test2"]
