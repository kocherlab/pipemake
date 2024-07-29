import os
import itertools

from pipemake.seqIO import SeqFileIO

from snakemake.io import glob_wildcards

class WildcardIO ():
	def __init__ (self, wildcard_str, wildcard_dict, sample_wildcard = '', **kwargs):

		# Confirm wildcards were assigned
		if not list(wildcard_dict):
			raise Exception(f'Unable to find wildcards within: {wildcard_str}')

		# Assign the basic arguments
		self.wildcard_str = wildcard_str
		self.wildcard_dict = wildcard_dict

		# Check if the wildcard dict is empty
		for wildcard_name, wildcard_values in self.wildcard_dict.items():
			if not wildcard_values: raise Exception(f'Unable to find wildcard values for: {wildcard_name}')

		if not sample_wildcard: self.samples = []
		else: self.samples = self.wildcard_dict[sample_wildcard]

	@classmethod
	def fromStr (cls, wildcard_str, **kwargs):

		wildcard_dict = {}

		# Use glob to retrieve wildcard entries
		wildcard_data = glob_wildcards(wildcard_str)

		# Loop the wildcard names - i.e. fields
		for wildcard_name in wildcard_data._fields:

			# Report the IDs in order without duplicates
			wildcard_dict[wildcard_name] = list(dict.fromkeys(getattr(wildcard_data, wildcard_name)))

		return cls(wildcard_str, wildcard_dict, **kwargs)

	def standardizedFiles (self, standardized_wildcard, **kwargs):

		print(f'wildcard_dict: {self.wildcard_dict}')

		print(kwargs)

		# Create list of the wildcard name and values
		wildcard_names, wildcard_values = zip(*self.wildcard_dict.items())
		
		# Format the wildcard str for the sample and standardized file
		for sample_wildcard_dict in [dict(zip(wildcard_names, v)) for v in itertools.product(*wildcard_values)]:
			sample_filename = self.wildcard_str.format(**sample_wildcard_dict)
			standardized_filename = standardized_wildcard.format(**sample_wildcard_dict)
			
			# Standardize the file
			sample_file = SeqFileIO.create(sample_filename)
			sample_file.standardize(standardized_filename, **kwargs)

	def returnPaths (self, copy_method = 'symbolic_link', **kwargs):
		if copy_method == 'copy': return []
		elif copy_method == 'symbolic_link': 
			path_name = os.path.dirname(self.wildcard_str)
			if not os.path.isdir(path_name): raise Exception(f'Unable to assign path: {path_name}')
			return [os.path.abspath(path_name)]
		else: raise Exception(f'Unsupported copy method: {copy_method}')

	def returnSamples (self):
		return self.samples