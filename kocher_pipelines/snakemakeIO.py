import os
import logging

from collections import defaultdict

class SnakePipelineIO (): 
	def __init__ (self, smkp_prefix, indent_style = '\t', overwrite = True, **kwargs):

		# Assign the basic arguments
		self._smkp_prefix = smkp_prefix
		self.smkp_filename = f'{smkp_prefix}.smk'

		# Confirm the file exists
		if os.path.isfile(self.smkp_filename) and not overwrite:
			raise IOError (f'SnakePipeline file already exists: {self.smkp_filename}. Please rename the SnakePipeline file or overwrite')
	
		# Assign the basic arguments
		self._pipe_file = open(self.smkp_filename, 'w')
		self._indent_style = indent_style

		# Create arguments to populate
		self._module_filenames = []
		self._all_input = []
		self._config_params = defaultdict(list)

	def __enter__ (self):
		return self

	def __exit__ (self, exc_type, exc_value, traceback):
		self.close()

	@classmethod
	def open (cls, *args, **kwargs):
		return cls(*args, **kwargs)
	
	def addSnakeModule (self, smkm_filename, **kwargs):

		# Open the snakemake module file
		smkm_file = SnakeFileIO.open(smkm_filename)

		# Create the rule dict
		smkm_rule_all_dict = smkm_file.returnRuleDict('all')

		# Check that the rule all dict is valid
		if 'input' not in smkm_rule_all_dict: raise Exception(f'Unable to parse rule all in: {smkm_filename}')

		# Create the config dict
		smkm_config_dict = smkm_file.returnRuleDict('config')

		### CAN THIS BE OPTIONAL?
		# Check that the config dict is valid 
		if 'params' not in smkm_config_dict: raise Exception(f'Unable to parse config in: {smkm_filename}')
		
		# Update the pipeline with the module
		self._module_filenames.append(smkm_filename)
		self._all_input.extend(smkm_rule_all_dict['input'])
		for _a, _v in smkm_config_dict.items(): self._config_params[_a].extend(_v)

		# Copy the pipeline module
		smkm_file.copy(out_filename = os.path.basename(smkm_filename), **kwargs)

		logging.info(f"Module added to pipeline: {smkm_filename}")

	def returnConfigParams (self):
		
		# Return the config parameters
		return self._config_params

	def writePipeline (self):

		# Assign the config file and working directory
		self._pipe_file.write(f"configfile: '{self._smkp_prefix}.yml'\n\n")
		self._pipe_file.write("workdir: config['workdir']\n\n")

		# Create the rule all block
		self._pipe_file.write('rule all:\n')
		self._pipe_file.write(f'{self._indent_style}input:\n')
		
		# Add the rule all input
		indented_input = [f'{self._indent_style}{self._indent_style}{_input}' for _input in self._all_input]
		self._pipe_file.write(',\n'.join(indented_input) + '\n\n')
		
		# Create the input block of modules
		for module_filename in self._module_filenames :
			self._pipe_file.write(f'include: "{module_filename}"\n')

		logging.info(f"Pipeline written successfully")
	
	def close(self):
		self._pipe_file.close()
	
class SnakeFileIO ():
	def __init__ (self, smk_filename, **kwargs):

		# Confirm the file exists
		if not os.path.isfile(smk_filename):
			raise IOError (f'Unable to open: {smk_filename}. Please confirm the file exists.')

		# Assign the basic arguments
		self.filename = smk_filename
		self._indent_style = None
		self._exclude_rules = ['config', 'refs', 'all']

		# Assign the indent style
		self._assignIndent()

	def _assignIndent (self):
		
		# Bool to store if within a rule
		in_rule_block = False

		with open(self.filename) as smk_file:
			for smk_line in smk_file:

				# Skip blank lines
				if not smk_line.strip(): continue

				# Confirm indent style once assigned
				if self._indent_style:

					# Skip rules as they have no indent
					if smk_line.startswith('rule'): continue
					
					if not smk_line.startswith(self._indent_style):
						raise Exception(f'Inconsistent indent style "{smk_line[0]}"')

					# Get the count of the indent style
					indent_count = len(smk_line) - len(smk_line.lstrip())

					# Calc the indent remainder
					indent_remainder = indent_count % len(self._indent_style)

					# Check indent remainder
					if indent_remainder:
						raise Exception(f'Indent style remainder: {indent_remainder}')
									
				# Assign the indent style
				elif in_rule_block:
					
					# Get the count of the leading whiespace
					leading_whitespace_count = (len(smk_line) - len(smk_line.lstrip()))

					# Check for table indent style
					if smk_line[0] == '\t':
						self._indent_style = '\t' * leading_whitespace_count
						logging.info(f"Indent style: Tab(s). Count: {leading_whitespace_count}")
					
					# Check for space indent style
					elif smk_line[0] == ' ': 
						self._indent_style = ' ' * leading_whitespace_count
						logging.info(f"Indent style: Spaces(s). Count: {leading_whitespace_count}")
					
					# Report error, if other
					else:
						raise Exception('Unable to assign indent style')

				# Check if within first rule block 
				if smk_line.startswith('rule'): in_rule_block = True
			
	@classmethod
	def open (cls, *args, **kwargs):
		return cls(*args, **kwargs)
	
	def copy (self, out_prefix = '', out_filename = '', out_dir = '', exclude_rules = [], **kwargs):
	
		# Populate the rules to exclude
		exclude_rules.extend(self._exclude_rules)

		# Craete bool to assign if within a exclusion block
		within_exclusion_block = False

		if not out_prefix and not out_filename:
			raise Exception(f'No output method given for snakefile')

		# Create the output file
		if out_prefix: out_path = f'{out_prefix}.smk'
		else: out_path = out_filename

		# Update the path with the output directory, if given
		if out_dir: out_path = os.path.join(out_dir, out_path)

		# Create the output file
		snakemake_output_file = open(out_path, 'w')

		# Parse the snakefile
		with open(self.filename) as smk_file:
			for smk_line in smk_file:

				# Split the line by the indent style (rule, attribute, param)
				smk_list = smk_line.rstrip().split(self._indent_style)

				# Include non-rule functions
				if smk_list[0] and smk_list[0].split(' ')[0] != 'rule': within_exclusion_block = False

				# Check if within a rule block
				elif smk_list[0]:

					# Confirm the rule block
					if smk_list[0].split(' ') == 'rule' and not smk_list[0].endswith(':'):
						raise Exception (f'Error copying rule: {smk_list[0]}')

					# Assign the rule name
					rule_name = smk_list[0].split(' ')[1].strip(':')
					
					# Assign the exclusion bool
					if rule_name in exclude_rules: within_exclusion_block = True
					else: within_exclusion_block = False	

				# Write output if not being excluded
				if not within_exclusion_block: snakemake_output_file.write(smk_line)

		logging.info(f"Copied snakemake module: {self.filename} to {out_path}")

	def returnRuleDict (self, rule_name):

		# Create the param dict
		param_dict = defaultdict(list)

		# Create str and int for param assignment check
		param_attribute_name = ''
		param_attribute_level = ''

		# Create bool to check if in the config block
		within_rule_block = False

		# Parse the snakefile
		with open(self.filename) as smk_file:
			for smk_line in smk_file:
				
				# Skip blank lines
				if not smk_line.strip(): continue

				# Split the line by the indent style (rule, attribute, param)
				smk_list = smk_line.rstrip().split(self._indent_style)

				# Check if within a rule
				if smk_list[0]:
					if rule_name in smk_list[0]: within_rule_block = True
					elif within_rule_block: break
					else: within_rule_block = False
					continue
				
				# Check if within the config block
				if within_rule_block:
					
					# Loop the snakefile line by level
					for smk_level, smk_item in enumerate(smk_list):
						if not smk_item: continue
						
						# Assign the param attribute
						if smk_item.endswith(':'):
							param_attribute_name = smk_item[:-1]
							param_attribute_level = smk_level

						# Assign the param
						else:
							if (smk_level - param_attribute_level) != 1:
								raise Exception (f'Parameter assignment error')
							param_dict[param_attribute_name].append(smk_item)

		logging.info(f"Created param dict for rule: {rule_name} from {self.filename}")

		return param_dict

