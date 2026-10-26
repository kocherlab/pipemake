import os
import re


def collectExistingPaths(arg_dict, omit=()):
    """Return the existing paths referenced by the string values of an argument dict.

    Wildcard values (e.g. "dir/{samples}_R1.fq.gz") are reduced to the directory
    before the first wildcard character, so wildcard inputs are also detected.
    """
    paths = set()
    for arg_name, arg_value in arg_dict.items():
        if arg_name in omit:
            continue

        values = arg_value if isinstance(arg_value, (list, tuple, set)) else [arg_value]
        for value in values:
            if not isinstance(value, str) or not value:
                continue

            # Reduce wildcard values to the literal directory prefix
            wildcard = re.search(r"[*?\[{]", value)
            if wildcard:
                value = os.path.dirname(value[: wildcard.start()])

            if value and os.path.exists(value):
                paths.add(value)

    return sorted(paths)


def confirmSafeToOverwrite(workflow_dir, input_paths=(), job_dir="pipemake"):
    """Raise an exception unless workflow_dir is safe to delete.

    A directory is only deleted if it looks like a previous pipemake workflow
    (it contains <job_dir>/pipeline.log), it is not a symlink, it is not the
    filesystem root, the home directory (or one of its parents), or the current
    directory (or one of its parents), and it does not contain any input path.
    """

    def isWithin(child, parent):
        return child == parent or child.startswith(parent.rstrip(os.sep) + os.sep)

    target = os.path.realpath(workflow_dir)
    cwd = os.path.realpath(os.getcwd())
    home = os.path.realpath(os.path.expanduser("~"))

    if os.path.islink(workflow_dir):
        raise Exception(f"Refusing to overwrite a symbolic link: {workflow_dir}")

    if not os.path.isdir(target):
        raise Exception(f"Refusing to overwrite, not a directory: {workflow_dir}")

    if target == os.path.realpath(os.sep):
        raise Exception("Refusing to overwrite the filesystem root")

    if isWithin(home, target):
        raise Exception(
            f"Refusing to overwrite {target}: it is, or contains, your home directory"
        )

    if isWithin(cwd, target):
        raise Exception(
            f"Refusing to overwrite {target}: it is, or contains, the current directory"
        )

    for input_path in input_paths:
        if isWithin(os.path.realpath(input_path), target):
            raise Exception(
                f"Refusing to overwrite {target}: it contains an input ({input_path})"
            )

    if not os.path.isfile(os.path.join(target, job_dir, "pipeline.log")):
        raise Exception(
            f"Refusing to overwrite {target}: it does not look like a previous pipemake workflow "
            f"(missing {os.path.join(job_dir, 'pipeline.log')})"
        )
