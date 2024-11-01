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

    @classmethod
    def fromWildcardStr(cls, wildcard_str="", **kwargs):
        return cls(
            WildcardIO.fromStr(wildcard_str, **kwargs),
            standardize_func="standardizedFiles",
            standardize_type="file",
        )

    @classmethod
    def fromTableFile(cls, table_filename="", **kwargs):
        return cls(
            TableIO.fromFilenameStr(table_filename, **kwargs),
            standardize_func="standardizedFiles",
            standardize_type="file",
        )

    @classmethod
    def fromFileStr(cls, input_filename="", **kwargs):
        return cls(
            FileIO.create(input_filename, **kwargs),
            standardize_func="standardize",
            standardize_type="file",
        )

    @classmethod
    def fromDirStr(cls, path_name="", **kwargs):
        return cls(
            DirIO.create(path_name, **kwargs),
            standardize_func="standardize",
            standardize_type="dir",
        )

    def standardize(self, **kwargs):
        if self.standardize_type == "file":
            self._standardize_file(**kwargs)
        elif self.standardize_type == "dir":
            self._standardize_dir(**kwargs)
        else:
            raise Exception(f"Standardize type not recognized: {self.standardize_type}")

    def _standardize_file(self, standardized_filename="", **kwargs):
        self.standardize_call(standardized_filename, **kwargs)

    def _standardize_dir(self, standardized_directory="", **kwargs):
        self.standardize_call(standardized_directory, **kwargs)

    def returnSamples(self, **kwargs):
        return self.processIO.returnSamples()

    def returnPaths(self, **kwargs):
        return self.processIO.returnPaths(**kwargs)


def standardizeInput(method="", args={}):
    # Create the standardization call
    if method == "wildcard-str":
        standardize_input_call = ProcessIO.fromWildcardStr(**args)
    elif method == "table-file":
        standardize_input_call = ProcessIO.fromTableFile(**args)
    elif method == "file-str":
        standardize_input_call = ProcessIO.fromFileStr(**args)
    elif method == "dir-str":
        standardize_input_call = ProcessIO.fromDirStr(**args)
    else:
        raise Exception(f"No standardization method given for: {method}")

    # Standardize the input
    standardize_input_call.standardize(**args)


def returnPaths(method="", args={}):
    # Create the standardization call
    if method == "wildcard-str":
        return_path_call = ProcessIO.fromWildcardStr(**args)
    elif method == "table-file":
        return_path_call = ProcessIO.fromTableFile(**args)
    elif method == "file-str":
        return_path_call = ProcessIO.fromFileStr(**args)
    elif method == "dir-str":
        return_path_call = ProcessIO.fromDirStr(**args)
    else:
        raise Exception(f"No standardization method given for: {method}")

    return return_path_call.returnPaths(**args)


def returnSamples(method="", args={}):
    # Create the return samples call
    if method == "wildcard-str":
        return_samples_call = ProcessIO.fromWildcardStr(**args)
    elif method == "table-file":
        return_samples_call = ProcessIO.fromTableFile(**args)
    else:
        raise Exception(f"No sample method given for: {method}")

    return return_samples_call.returnSamples()
