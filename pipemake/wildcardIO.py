import os
import logging
import itertools

from pipemake.fileIO import FileIO

from snakemake.io import glob_wildcards


class WildcardIO:
    def __init__(self, wildcard_str, wildcard_dict, sample_wildcards=[], **kwargs):
        # Confirm wildcards were assigned
        if not list(wildcard_dict):
            raise Exception(f"Unable to find wildcards within: {wildcard_str}")

        # Assign the basic arguments
        self.wildcard_str = wildcard_str
        self.wildcard_dict = wildcard_dict

        # Check if the wildcard dict is empty
        missing_wildcards = [
            _w_name for _w_name, _w_value in self.wildcard_dict.items() if not _w_value
        ]
        if missing_wildcards:
            # Report that all wildcards are missing, likely due to incorrect filenames
            if len(missing_wildcards) == len(self.wildcard_dict):
                raise Exception(
                    "Unable to assign by wildcards. Please confirm the filenames are correct."
                )

            # Report the missing wildcards, not sure if this is possible
            else:
                raise Exception(
                    f"Unable to find wildcard values for: {', '.join(missing_wildcards)}"
                )

        if not sample_wildcards:
            self.samples = {}
        else:
            self.samples = {_col: wildcard_dict[_col] for _col in sample_wildcards}

    @classmethod
    def fromStr(cls, wildcard_str, **kwargs):
        wildcard_dict = {}

        # Use glob to retrieve wildcard entries
        wildcard_data = glob_wildcards(wildcard_str)

        # Loop the wildcard names - i.e. fields
        for wildcard_name in wildcard_data._fields:
            # Report the IDs in order without duplicates
            wildcard_dict[wildcard_name] = list(
                dict.fromkeys(getattr(wildcard_data, wildcard_name))
            )

        return cls(wildcard_str, wildcard_dict, **kwargs)

    def standardizedFiles(self, standardized_wildcard, **kwargs):
        # Create list of the wildcard name and values
        wildcard_names, wildcard_values = zip(*self.wildcard_dict.items())

        # Format the wildcard str for the sample and standardized file
        for sample_wildcard_dict in [
            dict(zip(wildcard_names, v)) for v in itertools.product(*wildcard_values)
        ]:
            sample_filename = self.wildcard_str.format(**sample_wildcard_dict)

            # Confirm the sample file exists
            if not os.path.isfile(sample_filename):
                logging.warning(
                    f"Unable to find file: {sample_filename}. This may occur if wildcards are not always combined."
                )
                continue

            # Check if any standardized wildcards are missing
            missing_wildcards = [
                _w
                for _w in glob_wildcards(standardized_wildcard)._fields
                if _w not in sample_wildcard_dict
            ]

            # Confirm the wildcards provided are in the sample wildcard dict
            if missing_wildcards:
                raise Exception(
                    f"Unable to assign [{', '.join(missing_wildcards)}] in {standardized_wildcard}. Please confirm the same wildcards are used in the sample and standardized filenames."
                )

            standardized_filename = standardized_wildcard.format(**sample_wildcard_dict)

            # Standardize the file
            sample_file = FileIO.create(sample_filename)
            sample_file.standardize(standardized_filename, **kwargs)

    def returnPaths(self, copy_method="symbolic_link", **kwargs):
        if copy_method == "copy":
            return []
        elif copy_method == "symbolic_link":
            path_name = os.path.dirname(self.wildcard_str)
            if not os.path.isdir(path_name):
                raise Exception(f"Unable to assign path: {path_name}")
            return [os.path.abspath(path_name)]
        else:
            raise Exception(f"Unsupported copy method: {copy_method}")

    def returnSamples(self):
        return self.samples
