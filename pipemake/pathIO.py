import os
import shutil
import logging


class DirIO:
    def __init__(self, path_name, **kwargs):
        # Confirm the directory exists
        if not os.path.isdir(path_name):
            raise IOError(
                f"Unable to locate: {path_name}. Please confirm the directory exists."
            )

        # Convert the directory to an absolute path
        abs_path_name = os.path.abspath(path_name)

        # Confirm the directory using the absolute path
        if not os.path.isdir(abs_path_name):
            raise IOError(
                f"Unable to locate: {abs_path_name}. Please confirm the directory exists."
            )

        # Assign the basic arguments
        self.path_name = abs_path_name

        # Stores path if using a link
        self.link_path = ""

    @classmethod
    def create(cls, *args, **kwargs):
        return cls(*args, **kwargs)

    def standardize(
        self,
        standardized_directory,
        out_dir="",
        workflow_prefix="",
        work_dir="",
        copy_method="symbolic_link",
        **kwargs,
    ):
        # Assign the destination directory
        dest_directory = standardized_directory

        # Create path as needed
        if out_dir:
            dest_directory = os.path.join(out_dir, dest_directory)
        if workflow_prefix:
            dest_directory = os.path.join(workflow_prefix, dest_directory)
        if work_dir:
            dest_directory = os.path.join(work_dir, dest_directory)

        # Create the output directory, if needed
        if not os.path.exists(os.path.dirname(dest_directory)):
            os.makedirs(os.path.dirname(dest_directory))

        # Copy the directory
        if copy_method == "symbolic_link":
            os.symlink(self.path_name, dest_directory)
        elif copy_method == "copy":
            shutil.copytree(self.path_name, dest_directory)
        else:
            raise Exception(f"Unknown copy method: {copy_method}")

        logging.info(
            f"Standardized: {self.path_name} to {dest_directory}. Method: {copy_method}"
        )

    def returnPaths(self, copy_method="symbolic_link", **kwargs):
        if copy_method == "copy":
            return []
        elif copy_method == "symbolic_link":
            return [self.path_name]
        else:
            raise Exception(f"Unsupported copy method: {copy_method}")
