from pipemake.fileIO import FileIO, TableIO
from pipemake.wildcardIO import WildcardIO
from pipemake.pathIO import DirIO


class ProcessIO:
    def __init__(
        self, processIO, standardize_func=None, standardize_type=None, **kwargs
    ):
        self.processIO = processIO
        self.standardize_call = getattr(self.processIO, standardize_func)
        self.standardize_type = standardize_type
        self._kwargs = kwargs

    @classmethod
    def fromWildcardStr(cls, wildcard_str="", **kwargs):
        return cls(
            WildcardIO.fromStr(wildcard_str, **kwargs),
            standardize_func="standardizedFiles",
            standardize_type="file",
            **kwargs,
        )

    @classmethod
    def fromTableFile(cls, table_file="", **kwargs):
        return cls(
            TableIO.fromFilenameStr(table_file, **kwargs),
            standardize_func="standardizedFiles",
            standardize_type="file",
            **kwargs,
        )

    @classmethod
    def fromFileStr(cls, file_str="", **kwargs):
        return cls(
            FileIO.create(file_str, **kwargs),
            standardize_func="standardize",
            standardize_type="file",
            **kwargs,
        )

    @classmethod
    def fromDirStr(cls, dir_str="", **kwargs):
        return cls(
            DirIO.create(dir_str, **kwargs),
            standardize_func="standardize",
            standardize_type="dir",
            **kwargs,
        )

    def standardize(self):
        if self.standardize_type == "file":
            self._standardize_file(**self._kwargs)
        elif self.standardize_type == "dir":
            self._standardize_dir(**self._kwargs)
        else:
            raise Exception(f"Standardize type not recognized: {self.standardize_type}")

    def _standardize_file(self, standardized_filename="", **kwargs):
        self.standardize_call(standardized_filename, **kwargs)

    def _standardize_dir(self, standardized_directory="", **kwargs):
        self.standardize_call(standardized_directory, **kwargs)

    def returnSamples(self):
        return self.processIO.returnSamples()

    def returnPaths(self):
        return self.processIO.returnPaths(**self._kwargs)


def processInput(method="", args={}):
    # Convert method to lowercase and replace hyphens with underscores
    method = method.replace("-", "_").lower()

    # Create the standardization call
    if method == "wildcard_str":
        processed_input = ProcessIO.fromWildcardStr(**args)
    elif method == "table_file":
        processed_input = ProcessIO.fromTableFile(**args)
    elif method == "file_str":
        processed_input = ProcessIO.fromFileStr(**args)
    elif method == "dir_str":
        processed_input = ProcessIO.fromDirStr(**args)
    else:
        raise Exception(f"No standardization method given for: {method}")

    return processed_input
