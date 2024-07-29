import os
import gzip
import shutil
import logging

import pandas as pd

from snakemake.io import get_wildcard_names

def checkIfGzipped (filename):

	with open (filename, 'rb') as check_file:
		check_line = check_file.readline()
		bgzip_header = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00'
		gzip_header = b'\x1f\x8b'
		if check_line[:len(bgzip_header)] == bgzip_header: return True
		if check_line[:len(gzip_header)] == gzip_header: return True
	return False

class SeqFileIO ():
	def __init__ (self, seq_filename, **kwargs):

		# Confirm the file exists
		if not os.path.isfile(seq_filename):
			raise IOError (f'Unable to open: {seq_filename}. Please confirm the file exists.')

		# Convert the filename to an absolute path
		abs_seq_filename = os.path.abspath(seq_filename)

		# Confirm the file using the absolute path
		if not os.path.isfile(abs_seq_filename):
			raise IOError (f'Unable to open: {abs_seq_filename}. Please confirm the file exists.')

		# Assign the basic arguments
		self.filename = abs_seq_filename
		
		# Assign if the file is compressed
		self.gzipped = checkIfGzipped(self.filename)

		# Stores path if using a link
		self.link_path = ''
	
	@ property
	def args (self):
		
		# Create dict to store optional args
		arg_dict = {}
		
		# Assign potential args
		if self.link_path: arg_dict['bind'] = self.link_path
		
		return arg_dict

	@classmethod
	def create (cls, *args, **kwargs):
		return cls(*args, **kwargs)

	def standardize (self, standardized_filename, out_dir = '', workflow_prefix = '', work_dir = '', gzipped = None, copy_method = 'symbolic_link', **kwargs):

		print(standardized_filename, out_dir, workflow_prefix, work_dir, gzipped, copy_method)

		# Assign the destination filename
		dest_filename = standardized_filename

		# Create path as needed
		if out_dir: dest_filename = os.path.join(out_dir, dest_filename)
		if workflow_prefix: dest_filename = os.path.join(workflow_prefix, dest_filename)
		if work_dir: dest_filename = os.path.join(work_dir, dest_filename)

		# Create the output directory, if needed
		if not os.path.exists(os.path.dirname(dest_filename)):
			os.makedirs(os.path.dirname(dest_filename))

		# Copy files with the same gzipped status
		if gzipped == None or self.gzipped == gzipped:
			if copy_method == 'symbolic_link':
				self.link_path = os.path.dirname(self.filename)
				os.symlink(self.filename, dest_filename)
			elif copy_method == 'copy': shutil.copy(self.filename, dest_filename)
			else: raise Exception(f'Unknown copy method: {copy_method}')
			
		# Copy to a decompressed status
		elif self.gzipped and not gzipped:
			with gzip.open(self.filename, 'r') as f_in, open(dest_filename, 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)

		# Copy to a compressed status
		elif not self.gzipped and gzipped:
			with open(self.filename, 'r') as f_in, gzip.open(dest_filename, 'wt') as f_out:
				shutil.copyfileobj(f_in, f_out)

		# Raise exception (not possible?)
		else: raise Exception (f'Unable to copy file: {self.filename}')

		logging.info(f"Standardized: {self.filename} to {dest_filename}. Method: {copy_method}")

	def returnPaths (self, copy_method = 'symbolic_link', **kwargs):
		if copy_method == 'copy': return []
		elif copy_method == 'symbolic_link': return [os.path.dirname(self.filename)]
		else: raise Exception(f'Unsupported copy method: {copy_method}')

class SeqTableIO ():
	def __init__ (self, table_dataframe, sample_column = 'sample', **kwargs):

		# Confirm wildcards were assigned
		if not isinstance(table_dataframe, pd.DataFrame): 
			raise Exception(f'Table must be stored as a pandas DataFrame: {table_dataframe}')
		
		# Assign the basic arguments
		self._table_dataframe = table_dataframe
		self._sample_column = sample_column
		self._file_columns = set()
		self._table_columns = set([sample_column])

		if sample_column not in self._table_dataframe.columns:
			raise Exception(f'Unable to assign sample column ({self._sample_column}) in DataFrame columns: {self._table_dataframe.columns}')

		# Assign arguments from the columns of the dataframe
		for col in self._table_dataframe.columns:
			if col == self._sample_column: self.samples = list(self._table_dataframe[self._sample_column])
			elif ':' in col:
				wildcard_col, _ = col.split(':')
				self._file_columns.add(wildcard_col)
			else:
				self._table_columns.add(col)

		# Check the samples are unique
		if len(self.samples) != len(list(set(self.samples))):
			raise Exception(f'Sample names are not unique')

		if len(self._file_columns) > 1: 
			raise Exception(f'Files may exist within a single column')

		# Create a set of all the columns
		self._table_columns = self._table_columns.union(self._file_columns)
	
	@property
	def _file_column (self):
		if len(self._file_columns) > 1:
			raise Exception(f'Files may exist within a single column')
		return list(self._file_columns)[0]
	
	@classmethod
	def fromFilenameStr (cls, table_filename, sep = '\t', **kwargs):
		return cls(table_dataframe = pd.read_csv(table_filename, sep = sep, index_col = False, dtype = str).fillna(''), **kwargs)

	def standardizedFiles (self, standardized_wildcard, **kwargs):

		# Assign the standardized wildcard names
		standardized_wildcard_names = get_wildcard_names(standardized_wildcard)

		# Confirm the dataframe column are being standardized
		for col in self._table_columns:
			if col not in standardized_wildcard_names:
				raise Exception(f'Unable to assign column ({col}) for standardization: {standardized_wildcard}')
			standardized_wildcard_names.remove(col)
		
		# Confirm all wildcards were used
		if len(standardized_wildcard_names) > 0:
			raise Exception(f"Not all wildcards found in table: {', '.join(list(standardized_wildcard_names))}")

		# Loop the dataframe by row
		for _, row_dict in self._table_dataframe.iterrows():

			# Assign the non-filename wildcards for the row
			str_wildcards = {_c:_v for _c, _v in row_dict.items() if self._file_column not in _c}

			# Loop the row columns (i.e. indices)
			for file_col, sample_filename in row_dict.items():

				# Confirm the row column (i.e. index) is a filename
				if self._file_column not in file_col: continue

				# Skip if sample filename is blank (e.g. No R2 filename)
				if not sample_filename: continue

				# Confirm the file exists
				if not os.path.isfile(sample_filename):
					raise IOError(f"Unable to locate file. Column: {file_col}. Filename: {sample_filename}")

				# Create the wilcard dict for the sample
				file_wildcard_dict = dict([file_col.split(':')])
				file_wildcard_dict.update(str_wildcards)

				# Create the standardized filename
				standardized_filename = standardized_wildcard.format(**file_wildcard_dict)

				# Standardize the file
				sample_file = SeqFileIO.create(sample_filename)
				sample_file.standardize(standardized_filename, **kwargs)

	def returnPaths (self, copy_method = 'symbolic_link', **kwargs):
		if copy_method == 'copy': return []
		elif copy_method == 'symbolic_link':

			# Create list to store paths
			path_set = set()

			# Loop the dataframe by row
			for _, row_dict in self._table_dataframe.iterrows():
				for file_col, sample_filename in row_dict.items():

					# Confirm the row column (i.e. index) is a filename
					if self._file_column not in file_col: continue

					# Skip if sample filename is blank (e.g. No R2 filename)
					if not sample_filename: continue

					# Confirm the file exists
					if not os.path.isfile(sample_filename):
						raise IOError(f"Unable to locate file. Column: {file_col}. Filename: {sample_filename}")

					# Assign and store the directory
					path_name = os.path.dirname(sample_filename)
					if not os.path.isdir(path_name): raise Exception(f'Unable to assign path: {path_name}')
					path_set.add(path_name)

			return list(path_set)

		else: raise Exception(f'Unsupported copy method: {copy_method}')

	def returnSamples (self):
		return self.samples