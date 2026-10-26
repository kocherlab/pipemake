import os
import pytest

from pipemake.workflowIO import confirmSafeToOverwrite, collectExistingPaths


def makeWorkflowDir(path):
    """Create a directory that looks like a previous pipemake workflow"""
    (path / "pipemake").mkdir(parents=True)
    (path / "pipemake" / "pipeline.log").write_text("log")
    return path


def test_confirmSafeToOverwrite_previous_workflow(tmp_path):
    workflow = makeWorkflowDir(tmp_path / "workflow")
    confirmSafeToOverwrite(str(workflow))

    # Paths are normalized before being checked
    confirmSafeToOverwrite(f"{tmp_path}/other/../workflow/")


def test_confirmSafeToOverwrite_input_outside_workflow(tmp_path):
    workflow = makeWorkflowDir(tmp_path / "workflow")
    input_file = tmp_path / "input.vcf"
    input_file.write_text("vcf")
    confirmSafeToOverwrite(str(workflow), input_paths=[str(input_file)])


def test_confirmSafeToOverwrite_not_pipemake_workflow(tmp_path):
    other_dir = tmp_path / "results"
    other_dir.mkdir()
    (other_dir / "data.vcf").write_text("vcf")
    with pytest.raises(Exception, match="previous pipemake workflow"):
        confirmSafeToOverwrite(str(other_dir))


def test_confirmSafeToOverwrite_missing_directory(tmp_path):
    with pytest.raises(Exception, match="not a directory"):
        confirmSafeToOverwrite(str(tmp_path / "missing"))


def test_confirmSafeToOverwrite_file(tmp_path):
    a_file = tmp_path / "file.txt"
    a_file.write_text("text")
    with pytest.raises(Exception, match="not a directory"):
        confirmSafeToOverwrite(str(a_file))


@pytest.mark.parametrize("workflow_arg", [".", "./", ".."])
def test_confirmSafeToOverwrite_current_directory(tmp_path, monkeypatch, workflow_arg):
    # Even a real pipemake workflow is protected when it is the current directory,
    # or one of its parents
    workflow = makeWorkflowDir(tmp_path / "workflow")
    (workflow / "sub").mkdir()
    monkeypatch.chdir(workflow if workflow_arg != ".." else workflow / "sub")
    with pytest.raises(Exception, match="current directory"):
        confirmSafeToOverwrite(workflow_arg)


def test_confirmSafeToOverwrite_home_directory(tmp_path, monkeypatch):
    workflow = makeWorkflowDir(tmp_path / "workflow")
    monkeypatch.setenv("HOME", str(workflow))
    with pytest.raises(Exception, match="home directory"):
        confirmSafeToOverwrite(str(workflow))

    # A parent of the home directory is also protected
    monkeypatch.setenv("HOME", str(workflow / "user"))
    (workflow / "user").mkdir()
    with pytest.raises(Exception, match="home directory"):
        confirmSafeToOverwrite(str(workflow))


def test_confirmSafeToOverwrite_filesystem_root():
    with pytest.raises(Exception, match="filesystem root"):
        confirmSafeToOverwrite(os.sep)


def test_confirmSafeToOverwrite_contains_input(tmp_path):
    workflow = makeWorkflowDir(tmp_path / "workflow")
    input_file = workflow / "input.vcf"
    input_file.write_text("vcf")
    with pytest.raises(Exception, match="contains an input"):
        confirmSafeToOverwrite(str(workflow), input_paths=[str(input_file)])


def test_confirmSafeToOverwrite_symlink(tmp_path):
    workflow = makeWorkflowDir(tmp_path / "workflow")
    link = tmp_path / "link"
    link.symlink_to(workflow)
    with pytest.raises(Exception, match="symbolic link"):
        confirmSafeToOverwrite(str(link))


def test_collectExistingPaths(tmp_path):
    (tmp_path / "reads").mkdir()
    existing_file = tmp_path / "input.vcf"
    existing_file.write_text("vcf")

    arg_dict = {
        "vcf": str(existing_file),
        "wildcard": f"{tmp_path}/reads/{{samples}}_{{reads}}.fq.gz",
        "glob": f"{tmp_path}/reads/*.fq.gz",
        "list": [str(existing_file), "not_a_path"],
        "missing": str(tmp_path / "missing.txt"),
        "number": 1.0,
        "none": None,
        "workflow_dir": str(tmp_path),
    }

    paths = collectExistingPaths(arg_dict, omit=["workflow_dir"])
    assert paths == sorted([str(existing_file), str(tmp_path / "reads")])
