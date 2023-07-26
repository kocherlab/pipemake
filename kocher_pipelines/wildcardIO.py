import os
import sys
import copy
import itertools

from kocher_pipelines.seqIO import SeqFileIO

from snakemake.io import expand, glob_wildcards

class WildcardIO ():
	def __init__ (self, wildcard_str, wildcard_dict, sample_wildcard = '', **kwargs):

		# Confirm wildcards were assigned
		if not list(wildcard_dict):
			raise Exception(f'Unable to find wildcards within: {wildcard_str}')

		# Assign the basic arguments
		self.wildcard_str = wildcard_str
		self.wildcard_dict = wildcard_dict
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

'''
class WildcardIO2 ():
	def __init__ (self, wildcard_str, wildcard_dict, seq_format):

		# Assign the wildcard format class
		self.wildcard_class = getattr(sys.modules[__name__], seq_format)()

		# Check for unsupported wildcare names
		unsupported_wildcare_names = list(set(wildcard_dict) - self.wildcard_class.reserved_wildcards)
	
		# Confirm the keys are supported
		if unsupported_wildcare_names:
			raise Exception (f'Found wildcard entries that are not supported: {unsupported_wildcare_names}')

		# Assign the basic arguments
		self.wildcard_str = wildcard_str
		self.wildcard_dict = wildcard_dict

	@classmethod
	def fromStr (cls, wildcard_str, seq_format, **kwargs):

		wildcard_dict = {}

		# Use glob to retrieve wildcard entries
		wildcard_data = glob_wildcards(wildcard_str)

		# Loop the wildcard names - i.e. fields
		for wildcard_name in wildcard_data._fields:

			# Report the IDs in order without duplicates
			wildcard_dict[wildcard_name] = list(dict.fromkeys(getattr(wildcard_data, wildcard_name)))

		return cls(wildcard_str, wildcard_dict, seq_format, **kwargs)

	def createStandardizedFiles (self, dest, gzipped = True, **kwargs):

		# Create the destination directory
		if not os.path.exists(dest): os.makedirs(dest)

		# Create the the standardized files at the destination
		for seq_file in self.wildcard_class.assignStandardizedFiles(self.wildcard_str, self.wildcard_dict):
			seq_file.standardize(dest, gzipped = gzipped, **kwargs)

class fastq ():
	def __init__ (self):

		self.sample_wildcard = 'sample'
		self.read_wildcard = 'read'
		self.allow_duplcate_samples = True

		self.required_wildcards = [self.sample_wildcard]
		self.optional_wildcards = [self.read_wildcard]

		self.standard_templates = {(self.sample_wildcard,): f'{{{self.sample_wildcard}}}_R1.fq.gz',
								   (self.sample_wildcard, self.read_wildcard): f'{{{self.sample_wildcard}}}_{{{self.read_wildcard}}}.fq.gz'}
	@property
	def reserved_wildcards (self):
		return set(self.required_wildcards + self.optional_wildcards)

	def assignStandardizedFiles (self, wildcard_str, wildcard_dict):

		# Assign the expected number of files
		if self.read_wildcard not in wildcard_dict: file_count = 1
		else: file_count = len(wildcard_dict[self.read_wildcard])

		# Loop the samples
		for sample in wildcard_dict[self.sample_wildcard]:

			# Reduce the dict to only return the current sample
			sample_wildcard_dict = copy.deepcopy(wildcard_dict)
			sample_wildcard_dict[self.sample_wildcard] = [sample]

			# Check the file count
			sample_filenames = expand(wildcard_str, **sample_wildcard_dict)
			if len(sample_filenames) != file_count: raise Exception(f'Wildcard assignment error for {sample}. Unable to match expected number of files ({file_count}) using {sample_filenames}')

			# Assign the standardized filenames
			standardized_template = self.standard_templates[tuple(wildcard_dict)]
			standardized_filenames = expand(standardized_template, **sample_wildcard_dict)

			# Loop the sample/read
			for _sf, _stdf in zip(sample_filenames, standardized_filenames):
				yield SeqFileIO.create(_sf, 'fastq', standardized_filename = _stdf)
'''