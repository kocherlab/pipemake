import os
import yaml
import logging

from collections import defaultdict

class SnakePipelineIO (): 
	def __init__ (self, snakemake_job_prefix = '', pipeline_storage_dir  = '', pipeline_config_dir = '', work_dir = '', resource_yml = None, scale_threads = 0, scale_mem = 0, indent_style = '\t', overwrite = True, **kwargs):

		# Assign the basic arguments
		self._snakemake_job_prefix = snakemake_job_prefix

		# Assign the snakemake pipeline filename
		if snakemake_job_prefix.endswith('.smk'): self.smkp_filename = snakemake_job_prefix
		else: self.smkp_filename = f'{snakemake_job_prefix}.smk'

		# Confirm the file exists
		if os.path.isfile(self.smkp_filename) and not overwrite:
			raise IOError (f'SnakePipeline file already exists: {self.smkp_filename}. Please rename the SnakePipeline file or overwrite')
	
		# Assign the basic arguments
		self._pipe_file = open(self.smkp_filename, 'w')
		self._resource_yml = resource_yml
		self._scale_threads = scale_threads
		self._scale_mem = scale_mem
		self._mem_types = ['mem_b', 'mem_kb', 'mem_mb', 'mem_gb', 'mem_tb', 'mem_pb', 'mem_kib', 'mem_mib', 'mem_gib', 'mem_tib', 'mem_pib']
		self._indent_style = indent_style

		# Assign the module storage directory
		self._module_storage_dir = os.path.join(pipeline_storage_dir, 'modules')
		if not os.path.exists(self._module_storage_dir): os.makedirs(self._module_storage_dir)

		# Assign the module config directory
		self._module_config_dir = os.path.join(pipeline_config_dir, 'modules')
		if not os.path.exists(self._module_config_dir): os.makedirs(self._module_config_dir)

		# Assign the work directory
		self._work_dir = work_dir
		
		# Create arguments to populate
		self._module_filenames = []
		self._all_input = []
		self._config_params = defaultdict(list)
		self._resource_params = defaultdict(lambda: defaultdict(str))

	def __enter__ (self):
		return self

	def __exit__ (self, exc_type, exc_value, traceback):
		self.close()

	@classmethod
	def open (cls, *args, **kwargs):
		return cls(*args, **kwargs)

	def addSnakeModule (self, smkm_filename, **kwargs):

		# Assign the snakemake module filename
		if not smkm_filename.endswith('.smk'): smkm_filename = f'{smkm_filename}.smk'

		# Create the snakemake module filepath
		smkm_filepath = os.path.join(self._module_storage_dir, smkm_filename)

		if not os.path.isfile(smkm_filepath):
			raise IOError (f'Unable to open: {smkm_filepath}. Please confirm the file exists.')

		# Open the snakemake module file
		smkm_file = SnakeFileIO.open(smkm_filepath)

		# Create the rule dict
		smkm_rule_all_dict = smkm_file.returnBlockDict('all')

		# Check that the rule all dict is valid
		if 'input' not in smkm_rule_all_dict: raise Exception(f'Unable to parse rule all in: {smkm_filename}')

		# Create the config dict
		smkm_config_dict = smkm_file.returnBlockDict('config', pipeline_module_block = True)

		# Create the resource dict
		smkm_resource_dict = smkm_file.returnBlockDict('resources', pipeline_module_block = True, dict_w_defaults = True)

		### CAN THIS BE OPTIONAL?
		# Check that the config dict is valid 
		if 'params' not in smkm_config_dict: raise Exception(f'Unable to parse config in: {smkm_filename}')
		
		# Update the pipeline with the module
		self._module_filenames.append(smkm_filename)
		self._all_input.extend(smkm_rule_all_dict['input'])
		
		# Create the config dict
		for _a, _v in smkm_config_dict.items(): self._config_params[_a].extend(_v)
		
		# Create and scale the resource dict
		for _r, _adict in smkm_resource_dict.items():
			for _a, _v in _adict.items():
				if _a == 'threads' and self._scale_threads: _v = self._scaler(_v, self._scale_threads)
				elif 'mem_' in _a and _a.lower() in self._mem_types and self._scale_mem: _v = self._scaler(_v, self._scale_mem)
				self._resource_params[_r][_a] = _v

		# Copy the pipeline module
		smkm_file.copy(out_filename = os.path.basename(smkm_filename), out_dir = self._module_config_dir, **kwargs)

		logging.info(f"Module added to pipeline: {smkm_filename}")

	def writeConfig (self, pipeline_args):
		
		# Create yml dicts
		yml_config_dict = {'workdir': self._work_dir}
		yml_resource_dict = {'resources': {}}

		# Populate the config yml dict with the pipeline args
		for group, config_list in self._config_params.items():
			if group == 'params':
				for config in config_list:
					yml_config_dict[config] = pipeline_args[config.replace('-', '_')]
			else:
				group_dict = {}
				for config in config_list:
					group_dict[config] = pipeline_args[config.replace('-', '_')]
				yml_config_dict[group] = group_dict

		# Populate the resource yml dict
		for rule_name, rule_resource_dict in self._resource_params.items():
			if not rule_resource_dict: continue

			# Create the rule resource dict
			yml_resource_dict['resources'][rule_name] = {}

			# Populate the rule resource dict
			for resource_name, resource_value in rule_resource_dict.items():
				yml_resource_dict['resources'][rule_name][resource_name] = resource_value

		# Check if the resource yml should be a separate file
		if not self._resource_yml: yml_config_dict.update(yml_resource_dict)
		else: self._createYml(yml_resource_dict, f"{self._snakemake_job_prefix}.resources.yml")

		self._createYml(yml_config_dict, f"{self._snakemake_job_prefix}.yml")

	def writePipeline (self):

		# Assign the config file(s) and working directory
		self._pipe_file.write(f"configfile: '{self._snakemake_job_prefix}.yml'\n")
		if self._resource_yml: self._pipe_file.write(f"configfile: '{self._snakemake_job_prefix}.resources.yml'\n")
		self._pipe_file.write("\nworkdir: config['workdir']\n\n")

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
	
	@staticmethod
	def _scaler (value, scaler, value_type = float, output_type = int, **kwargs):

		# Confirm a scaler was given
		if not scaler: raise Exception(f'No scaler given for: {value}')

		# Convert the value to the correct type
		try: value = value_type(value)
		except: raise Exception(f'Unable to convert value to type: {value}')

		# Scale and return the input
		return output_type(value * scaler)

	@staticmethod
	def _createYml (yml_dict, yml_filename, **kwargs):

		class MyDumper(yaml.Dumper):
			def increase_indent(self, flow=False, indentless=False):
				return super(MyDumper, self).increase_indent(flow, False)

		yml_file = open(yml_filename, 'w')
		yml_file.write(yaml.dump(yml_dict, Dumper = MyDumper, sort_keys = False))
		yml_file.close()
		
class SnakeFileIO ():
	def __init__ (self, smk_filename, **kwargs):

		# Confirm the file exists
		if not os.path.isfile(smk_filename):
			raise IOError (f'Unable to open: {smk_filename}. Please confirm the file exists.')

		# Assign the basic arguments
		self.filename = smk_filename
		self._indent_style = None
		self._exclude_rules = ['all']
		
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
				if smk_line.split(' ')[0] == 'rule': in_rule_block = True
			
	@classmethod
	def open (cls, *args, **kwargs):
		return cls(*args, **kwargs)
	
	def copy (self, out_prefix = '', out_filename = '', out_dir = '', exclude_rules = [], **kwargs):
	
		# Populate the rules to exclude
		exclude_rules.extend(self._exclude_rules)

		# Craete bool to assign if within a exclusion block
		within_exclusion_block = False

		if not out_prefix and not out_filename:
			raise Exception(f'No output method given for snakefile: {self.filename}')

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

				# Assign the class name
				class_name = smk_list[0].split(' ')[0]

				# Include non-rule functions
				if smk_list[0] and class_name not in ['rule', 'module']: within_exclusion_block = False
	
				# Check if within a rule block
				elif smk_list[0]:

					# Confirm the module block
					if class_name != 'rule' and not smk_list[0].endswith(':'):
						raise Exception (f'Error copying rule: {smk_list[0]}')

					# Assign the rule name
					block_name = smk_list[0].split(' ')[1].strip(':')
					
					# Assign the exclusion bool
					if block_name in exclude_rules: within_exclusion_block = True
					elif class_name == 'module': within_exclusion_block = True
					else: within_exclusion_block = False

					if class_name == 'module' and not within_exclusion_block:
						raise Exception (f'Module copy error: {smk_list[0]}')

				# Write output if not being excluded
				if not within_exclusion_block: 
					snakemake_output_file.write(smk_line)

		logging.info(f"Copied snakemake module: {self.filename} to {out_path}")

	def returnBlockDict (self, block_name, pipeline_module_block = False, dict_w_defaults = False):

		# Create the param dict
		if not dict_w_defaults: param_dict = defaultdict(list)
		else: param_dict = defaultdict(lambda: defaultdict(str))

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
					if not pipeline_module_block and block_name in smk_list[0]: within_rule_block = True
					elif pipeline_module_block and smk_list[0].startswith('module') and block_name in smk_list[0]: within_rule_block = True
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
							param_attribute_name = smk_item[:-1].strip()
							param_attribute_level = smk_level

						# Assign the param
						else:
							if (smk_level - param_attribute_level) != 1:
								raise Exception (f'Parameter assignment error')

							# Assign the default value, if possible
							if ':' not in smk_item: smk_value = None
							else: 
								smk_item, smk_value = smk_item.split(':')
								smk_item = smk_item.strip()
								smk_value = smk_value.strip()

							# Check if default values are expected and assigned
							if dict_w_defaults and not smk_value:
								raise Exception (f'No default value assigned for: {param_attribute_name}: {smk_item}')
							
							# Assign the param to the dict
							if dict_w_defaults: param_dict[param_attribute_name][smk_item] = smk_value
							else: param_dict[param_attribute_name].append(smk_item)

		logging.info(f"Created param dict for rule: {block_name} from {self.filename}")

		return param_dict
	