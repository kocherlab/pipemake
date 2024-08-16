from pipemake.seqIO import SeqFileIO, SeqTableIO
from pipemake.wildcardIO import WildcardIO

class ProcessIO ():
	def __init__ (self, processIO, standardize_func = None, **kwargs):
		self.processIO = processIO
		self.standardize_call = getattr(self.processIO, standardize_func)
		#self.samples = []

	@classmethod
	def fromWildcardStr (cls, wildcard_str = '', **kwargs):
		return cls(WildcardIO.fromStr(wildcard_str, **kwargs), standardize_func = 'standardizedFiles')

	@classmethod
	def fromTableFile (cls, table_filename = '', **kwargs):
		return cls(SeqTableIO.fromFilenameStr(table_filename, **kwargs), standardize_func = 'standardizedFiles')

	@classmethod
	def fromFileStr (cls, input_filename = '', **kwargs):
		return cls(SeqFileIO.create(input_filename, **kwargs), standardize_func = 'standardize')

	def standardize (self, standardized_filename = '', **kwargs):
		self.standardize_call(standardized_filename, **kwargs)

	def returnSamples (self, **kwargs):
		return self.processIO.returnSamples()

	def returnPaths (self, **kwargs):
		return self.processIO.returnPaths(**kwargs)

def standardizeInput (method = '', args = {}):
	
	# Create the standardization call
	if method == 'wildcard-str': standardize_input_call = ProcessIO.fromWildcardStr(**args)
	elif method == 'table-file': standardize_input_call = ProcessIO.fromTableFile(**args)
	elif method == 'file-str': standardize_input_call = ProcessIO.fromFileStr(**args)
	else: raise Exception(f'No standardization method given for: {method}')
	
	# Standardize the input
	standardize_input_call.standardize(**args)

def returnPaths (method = '', args = {}):

	# Create the standardization call
	if method == 'wildcard-str': return_path_call = ProcessIO.fromWildcardStr(**args)
	elif method == 'table-file': return_path_call = ProcessIO.fromTableFile(**args)
	elif method == 'file-str': return_path_call = ProcessIO.fromFileStr(**args)
	else: raise Exception(f'No standardization method given for: {method}')
	
	return return_path_call.returnPaths(**args)

def returnSamples (method = '', args = {}):

	# Create the return samples call
	if method == 'wildcard-str': return_samples_call = ProcessIO.fromWildcardStr(**args)
	elif method == 'table-file': return_samples_call = ProcessIO.fromTableFile(**args)
	elif method == 'file-str': raise Exception('Not implemented')
	else: raise Exception(f'No standardization method given for: {method}')
	
	return return_samples_call.returnSamples()


	


