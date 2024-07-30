import os
import re
import copy
import yaml
import shutil
import logging

from collections import defaultdict

from pipemake.singularityIO import Singularity

class SnakePipelineIO (): 
	def __init__ (self, workflow_prefix = '', pipeline_storage_dir  = '', pipeline_job_dir = '', work_dir = '', singularity_dir = '', resource_yml = None, scale_threads = 0, scale_mem = 0, indent_style = '\t', overwrite = True, **kwargs):

		# Assign the basic arguments
		self._workflow_prefix = workflow_prefix
				
		# Assign the snakemake pipeline filename
		if workflow_prefix.endswith('.smk'): self.smkp_filename = workflow_prefix
		else: self.smkp_filename = f'{workflow_prefix}.smk'

		# Check if the overwrite is a bool or None
		if not isinstance(overwrite, bool): raise Exception (f'Invalid overwrite type: {overwrite}')

		# Confirm the file exists
		if os.path.isfile(self.smkp_filename) and not overwrite:
			raise IOError (f'SnakePipeline file already exists: {self.smkp_filename}. Please rename the SnakePipeline file or overwrite')

		# Confirm the resource yml is bool or None
		if not isinstance(resource_yml, (bool, type(None))): raise Exception (f'Invalid resource yml type: {resource_yml}')

		# Confirm the scale threads and mem are integers
		if not isinstance(scale_threads, float): raise Exception (f'Invalid scale threads type: {scale_threads}')
		if not isinstance(scale_mem, float): raise Exception (f'Invalid scale mem type: {scale_mem}')

		# Confirm the indent style
		if indent_style not in [' ', '\t']: raise Exception (f'Specified indent style not supported: {indent_style}')
	
		# Assign the basic arguments
		self._pipe_file = open(self.smkp_filename, 'w')
		self._resource_yml = resource_yml
		self._scale_threads = scale_threads
		self._scale_mem = scale_mem
		self._mem_types = ['mem_b', 'mem_kb', 'mem_mb', 'mem_gb', 'mem_tb', 'mem_pb', 'mem_kib', 'mem_mib', 'mem_gib', 'mem_tib', 'mem_pib']
		self._indent_style = indent_style

		# Assign the module storage directory, and confirm it exists
		self._module_storage_dir = os.path.join(pipeline_storage_dir, 'modules')
		if not os.path.exists(self._module_storage_dir): raise IOError (f'Unable to open: {self._module_storage_dir}. Please confirm the directory exists.')

		# Assign the module config directory
		self._module_job_dir = os.path.join(pipeline_job_dir, 'modules')
		if not os.path.exists(self._module_job_dir): os.makedirs(self._module_job_dir)

		# Assign the backup directory
		self._backup_dir = os.path.join(pipeline_job_dir, 'backups')
		if not os.path.exists(self._backup_dir): os.makedirs(self._backup_dir)
		
		# Assign the work directory
		self._work_dir = work_dir
		
		# Create arguments to populate
		self._module_filenames = []
		self._output_list = []
		self._config_params = set()
		self._resource_params = defaultdict(lambda: defaultdict(int))
		self._singularity_dir = singularity_dir
		self._pipeline_singularity_dict = {}

	def __enter__ (self):
		return self

	def __exit__ (self, exc_type, exc_value, traceback):
		self.close()

	@classmethod
	def open (cls, *args, **kwargs):
		return cls(*args, **kwargs)

	def addModule (self, module_filename, **kwargs):

		# Update the module filename, if necessary
		if not module_filename.endswith('.smk'): module_filename = f'{module_filename}.smk'

		# Create the filepath to the stored module
		storage_module_path = os.path.join(self._module_storage_dir, module_filename)

		# Confirm the stored module exists
		if not os.path.isfile(storage_module_path):
			raise IOError (f'Unable to open: {storage_module_path}. Please confirm the file exists.')		

		# Open and process the stored module file
		smk_module = SnakeFileIO.open(storage_module_path, singularity_dir = self._singularity_dir)

		# Loops the file rules in test
		for rule in smk_module._file_rules:

			# Save the config params and output list from the rule
			self._output_list.extend(rule._rule_output_list)
			self._config_params.update(rule._rule_config_params)

			# Scale and save the resource params
			for resource_type, resource_value in rule._rule_resource_params.items():
				if resource_type == 'threads' and self._scale_threads: resource_value = self._scaler(resource_value, self._scale_threads)
				elif 'mem_' in resource_type and resource_type.lower() in self._mem_types and self._scale_mem: resource_value = self._scaler(resource_value, self._scale_mem)
				self._resource_params[rule.rule_name][resource_type] = resource_value

		# Create the job module filepath
		storage_module_path = os.path.join(self._module_job_dir, module_filename)

		# Create the job module file
		smk_module.write(storage_module_path)

		# Add the module filename to the list
		self._module_filenames.append(module_filename)

		# Store the singularity containers from the module
		for singularity_path, singularity_container in smk_module._file_singularity_dict.items():
			if singularity_path not in self._pipeline_singularity_dict:
				self._pipeline_singularity_dict[singularity_path] = singularity_container

		logging.info(f"Module added to pipeline: {module_filename}")

	def buildSingularityContainers (self):
		
		# Loop the singularity containers
		for singularity_container in self._pipeline_singularity_dict.values():

			# Download the singularity container
			singularity_container.download()

	def writeConfig (self, pipeline_args):

		# Create a list of pipeline args to make sure they are all used
		pipeline_args_unused = set(pipeline_args.keys())
		
		# Create yml dicts
		yml_config_dict = {}
		yml_resource_dict = {'resources': {}}

		# Check if the workdir is in the pipeline args
		if self._work_dir: yml_config_dict['workdir'] = self._work_dir

		# Create a list of config params to sort by length
		config_params_list = list(self._config_params)
		config_params_list.sort(key = lambda t: len(t))

		print(pipeline_args)
		
		# Populate the config yml dict with the pipeline args
		for config_param in config_params_list:

			# Check if the config param is within a group
			if len(config_param) > 1: group, config_arg = config_param
			
			# Assign the config arg from the config param
			else:
				group = None
				config_arg = config_param[0]

			# Create the equivalent pipeline arg for the config arg
			pipeline_arg = config_arg.replace('-', '_')

			# Raise an error if the pipeline arg is not in the pipeline args
			if pipeline_arg not in pipeline_args: raise Exception (f'Pipeline arg not found: {pipeline_arg}')

			# Remove the pipeline arg from the unused list
			pipeline_args_unused.remove(pipeline_arg)
				
			# Check if the group is not in the yml dict
			if group and group not in yml_config_dict: yml_config_dict[group] = {}

			if not group: yml_config_dict[config_arg] = pipeline_args[pipeline_arg]
			else: yml_config_dict[group][config_arg] = pipeline_args[pipeline_arg]

		# Check if there are unused pipeline args
		#if pipeline_args_unused:

			# Report the unused pipeline args
			#for unused_arg in pipeline_args_unused: logging.warning(f"Configuration argument not found in Snakemake modules: {unused_arg}")

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
		else: self._createYml(yml_resource_dict, f"{self._workflow_prefix}.resources.yml")

		self._createYml(yml_config_dict, f"{self._workflow_prefix}.yml")


	def writePipeline (self):

		# Assign the config file(s) and working directory
		self._pipe_file.write(f"configfile: '{self._workflow_prefix}.yml'\n")
		if self._resource_yml: self._pipe_file.write(f"configfile: '{self._workflow_prefix}.resources.yml'\n")
		if self._work_dir: self._pipe_file.write("\nworkdir: config['workdir']\n\n")

		# Create the rule all block
		self._pipe_file.write('rule all:\n')
		self._pipe_file.write(f'{self._indent_style}input:\n')
		
		# Add the rule all input
		indented_input = [f'{self._indent_style}{self._indent_style}{_output}' for _output in self._output_list]
		self._pipe_file.write(',\n'.join(indented_input) + '\n\n')
		
		# Create the input block of modules
		for module_filename in self._module_filenames:
			self._pipe_file.write(f'include: "{os.path.join(self._module_job_dir, module_filename)}"\n')

		logging.info(f"Pipeline written successfully")
	
	def close(self):
		self._pipe_file.close()

		# Assign the basename of the workflow prefix
		workflow_basename = os.path.basename(self._workflow_prefix)

		# Create backups of the pipeline files
		shutil.copy(f"{self._workflow_prefix}.smk", os.path.join(self._backup_dir, f"{workflow_basename}.smk.bkp"))
		shutil.copy(f"{self._workflow_prefix}.yml", os.path.join(self._backup_dir, f"{workflow_basename}.yml.bkp"))
		if os.path.isfile(f"{self._workflow_prefix}.resources.yml"): 
			shutil.copy(f"{self._workflow_prefix}.resources.yml", os.path.join(self._backup_dir, f"{workflow_basename}.resources.yml.bkp"))

		logging.info(f"Pipeline backups created successfully")
	
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
	def __init__ (self, smk_filename, singularity_dir = '', **kwargs):

		# Confirm the file exists
		if not os.path.isfile(smk_filename):
			raise IOError (f'Unable to open: {smk_filename}. Please confirm the file exists.')

		# Assign the basic arguments
		self.filename = smk_filename
		self._singularity_dir = singularity_dir
		self._indent_style = None
		self._output_rule = 'all'
		self._exclude_rules = [self._output_rule]
		self._file_rules = []
		self._non_rule_text = ''
		
		# Assign the indent style
		self._assignIndent()

		# Parse the snakefile
		self._parseFile()
	
	@property
	def _file_singularity_dict (self):
		file_singularity_dict = {}
		for rule in self._file_rules:
			for singularity_path, singularity_container in rule._rule_singularity_dict.items():
				if singularity_path not in file_singularity_dict:
					file_singularity_dict[singularity_path] = singularity_container
		return file_singularity_dict

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
					if smk_line.startswith('rule') or smk_line.startswith('checkpoint') or smk_line.startswith('def'): continue
					
					if not smk_line.startswith(self._indent_style):
						raise Exception(f'Inconsistent indent style "{smk_line[0]}" in: {self.filename}')

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
			
	def _parseFile (self):

		# Create a string to store the rule block
		rule_block = ''

		# Create bool to check if in the config block
		within_rule_block = False

		# Parse the snakefile
		with open(self.filename) as smk_file:
			for smk_line in smk_file:

				# Skip blank lines
				if not smk_line.strip(): continue

				# Split the line by the indent style (rule, attribute, param)
				smk_list = smk_line.rstrip().split(self._indent_style)

				# Assign if line is a rule or checkpoint
				is_rule_line = ((smk_list[0].startswith('rule') or smk_list[0].startswith('checkpoint')) and smk_line.rstrip().endswith(':'))

				# Check if the line is the end of a rule
				if smk_list[0] and not is_rule_line:
					within_rule_block = False
					
				# Check if the line is the start of a rule
				if is_rule_line:

					within_rule_block = True

					# Check if a previous rule block was found
					if rule_block:
						self._file_rules.append(SnakeRuleIO.read(rule_block, singularity_dir = self._singularity_dir, indent_style = self._indent_style, output_rule = self._output_rule))
						rule_block = ''
				
				# Add the line to the rule block
				if within_rule_block and smk_line.strip(): rule_block += smk_line

				# Save text that is not within a rule
				if not within_rule_block and smk_line.strip(): self._non_rule_text += smk_line

			# Parse the final rule block
			if rule_block: 
				self._file_rules.append(SnakeRuleIO.read(rule_block, singularity_dir = self._singularity_dir, indent_style = self._indent_style, output_rule = self._output_rule))

	def write (self, filename):

		# Open the file
		with open(filename, 'w') as smk_file:

			# Check if there is any non-rule text
			if self._non_rule_text: smk_file.write(self._non_rule_text)

			# Loop the file rules
			for rule in self._file_rules:

				# Write the rule, if included in the output
				if rule.in_output: 
					smk_file.write(f'\n{rule}')

	@classmethod
	def open (cls, *args, **kwargs):
		return cls(*args, **kwargs)

class SnakeRuleIO ():
	def __init__ (self, rule_list, singularity_dir = '', indent_style = None, output_rule = '', **kwargs):

		# Assign the rule argument to populate
		self.rule_name = rule_list[0].split()[1][: -1]
		self.in_output = False if self.rule_name == output_rule else True
		self._singularity_dir = singularity_dir
		self._indent_style = indent_style
		self._rule_config_params = set()
		self._rule_singularity_dict = {}
		self._rule_resource_params = defaultdict(int)
		self._rule_text = rule_list[0] + '\n'
		self._rule_output_list = []
		
		# Assign the fixed arguments
		self._resource_attributes = ['threads', 'resources']
		
		# Parse the rule
		self._parseRule(rule_list[1:])

	def __str__ (self):
		return self._rule_text
	
	def __repr__ (self):
		return self._rule_text

	def _parseRule (self, rule_list):

		def processAttribute (rule_attribute, rule_line):

			# Assign the rule attribute
			processed_attribute = SnakeAttributeIO.process(self.rule_name, rule_attribute, rule_line, singularity_dir = self._singularity_dir, indent_style = self._indent_style)

			# Update the SnakeFile, if a resource or container
			if processed_attribute.is_resource: rule_line = processed_attribute.updateSnakeResource()
			elif processed_attribute.is_container: rule_line = processed_attribute.updateSnakeContainer()

			# Update the rule params
			if processed_attribute.is_resource: self._rule_resource_params.update(processed_attribute.updateParams())
			elif processed_attribute.is_container: self._rule_singularity_dict.update(processed_attribute.updateParams())
			else: self._rule_config_params.update(processed_attribute.updateParams())

			# Check if the rule is output, confirm the attribute is input, and update the rule output
			if not self.in_output and rule_attribute == 'input':
				output_line = rule_line.strip()[:-1] if rule_line.strip().endswith(',') else rule_line.strip()
				self._rule_output_list.append(output_line)
	
			return rule_line

		# Assign the current rule attribute
		rule_attribute = ''

		# Loop the rule block and parse the attributes
		for rule_line in rule_list:
			if not rule_line.strip(): continue

			# Split the line by the indent style
			rule_line_list = rule_line.rstrip().split(self._indent_style)

			# Assign the attribute level
			attribute_level = _attributeLevel(rule_line_list)

			# Check if the line is a rule attribute (without an assigned value)
			if attribute_level == 1 and rule_line.rstrip().endswith(':'):
				rule_attribute = rule_line_list[1][:-1]
			
			# Check if the line is a rule attribute (with an assigned value)
			elif attribute_level == 1 and not rule_line.rstrip().endswith(':'):
				
				# Split the attribute, confirm a value is assigned
				split_attribute = [_a.strip() for _a in rule_line_list[1].split(':', 1)]
				if len(split_attribute) != 2: raise Exception (f'Rule attribute error. Unable to assign attribute: {rule_line_list[1]}')

				# Process the attribute
				rule_line = processAttribute(split_attribute[0], rule_line)

			# Check if the line is a atrribute value line and if so, process the attribute
			elif attribute_level == 2: 
				rule_line = processAttribute(rule_attribute, rule_line)

			# Raise an error if the attribute level is not 1 or 2
			else: raise Exception (f'Attribute level error: {rule_line}')

			# Update the rule text
			self._rule_text += rule_line + '\n'

	@classmethod
	def read (cls, rule_str, *args, **kwargs):
		return cls(rule_str.splitlines(), *args, **kwargs)
	
class SnakeAttributeIO ():
	def __init__ (self, rule_name, attribute_type, attribute_text, singularity_dir = '', indent_style = None, resource_attributes = ['threads', 'resources'], **kwargs):

		# Assign the required attributes
		self._rule_name = rule_name
		self._type = attribute_type
		self._singularity_dir = singularity_dir
		self._indent_style = indent_style
		self._original_text = attribute_text
		self._multiline_text = True if attribute_text.endswith(',') else False
		self._resource_attributes = resource_attributes
		self.is_resource = True if attribute_type in self._resource_attributes else False
		self.is_container = True if attribute_type == 'singularity' else False
		self._resource_assignment_type = ':'
		self._container = None

		# Assign the config attributes
		self._resource_assignment_dict = defaultdict(int)
		self._resource_replacment_dict = defaultdict(str)
		self._config_assignment_set = set()
		
		# Process the attribute
		if self.is_container: self._parseContainer()
		elif self.is_resource: self._parseResource()
		else: self._parseConfig()

	@classmethod
	def process (cls, *args, **kwargs):
		return cls(*args, **kwargs)

	def _parseResource (self):

		# Create the working text
		working_text = copy.deepcopy(self._original_text) if not self._multiline_text else self._original_text[:-1]

		# Create a string to store the resource name
		resource_name = ''

		# Assign the resource as threads, for a single argument assignment
		if self._type == 'threads':

			# Confirm the threads is threads assignment
			if _attributeLevel(working_text.split(self._indent_style)) != 1:
				raise Exception (f'Config threads error. Incorrect number of indents: {self._original_text}')

			# Confirm a single colon sign is present
			if working_text.count(':') != 1:
				raise Exception (f'Config threads error. Expected colon: {self._original_text}')
			
			# Update the resource name and working text, by splitting the text by the assignment type
			resource_name, working_text = [_v.strip() for _v in working_text.strip().split(self._resource_assignment_type, 1)]

		# Check if the resource is a resources, for single or multiple argument assignment
		elif self._type == 'resources':

			# Confirm the resource is resource assignment
			if _attributeLevel(working_text.split(self._indent_style)) != 2:
				raise Exception (f'Config resource error. Incorrect number of indents: {self._original_text}')

			# Confirm a single equal sign is present
			if working_text.count('=') != 1:
				raise Exception (f'Config resource error. Expected equals siqn: {self._original_text}')

			# Update the resource assignment type
			self._resource_assignment_type = '='

			# Update the resource name and working text, by splitting the text by the assignment type
			resource_name, working_text = [_v.strip() for _v in working_text.strip().split(self._resource_assignment_type, 1)]

		# Check for a simple resource (integer) assignment
		if working_text.isdigit():

			# Confirm the resource has a name
			if not resource_name: raise Exception (f'Config resource error. No resource name: {self._original_text}')

			# Add the resource assignment to the dict
			self._resource_assignment_dict[resource_name] = int(working_text)
			self._resource_replacment_dict[working_text] = f"config['resources']['{self._rule_name}']['{resource_name}']"
		
		# Check for a complex resource assignment
		else:

			# Find all occurences of the word resources followed by one or more sets of square brackets ending with a colon
			for resource_match in re.finditer(r'resources(\[(.*?)\])+\:', working_text):

				# Get the resource str
				resource_str = resource_match.group(0)

				# Split the param str by either square brackets using regex
				resource_list = [_p for _p in re.split(r'\[|\]', resource_str[:-1].replace('"','').replace("'",'')) if _p][1:]

				# Check that the resource list is a single argument
				if len(resource_list) != 1: raise Exception (f'Config resource error. Too many arguments in complex assignment: {self._original_text}')

				# Assign the resource argument
				resource_name = resource_list[0]

				# Get text after resource_str in the working_text
				resource_suffix = working_text[working_text.find(resource_str) + len(resource_str):]

				# Get the first match between curly brackets without returning the brackets
				resource_suffix_match = re.search(r' ?{(.*?)}', resource_suffix)

				# Assign the match str and value
				resource_value_str = resource_suffix_match.group(0)
				resource_value = resource_suffix_match.group(1)

				# Add the resource assignment to the dict
				self._resource_assignment_dict[resource_name] = int(resource_value)
				self._resource_replacment_dict[resource_str + resource_value_str] = f"config['resources']['{self._rule_name}']['{resource_name}']"

	def _parseConfig (self):

		# Find all occurences of the word config follpwed by one or more sets of square brackets
		for config_match in re.finditer(r'config(\[(.*?)\])+', self._original_text):

			# Get the config match str
			config_str = config_match.group(0)

			# Split the param str by either square brackets using regex
			config_list = [_p for _p in re.split(r'\[|\]', config_str.replace('"','').replace("'",'')) if _p][1:]

			# Check if the match is associated with a value
			config_has_value = False if self._original_text.find(f'{config_str}:') == -1 else True

			# Check if the param does not have a value, i.e. config argument
			if config_has_value: raise Exception (f'Config param error. Has value: {self._original_text}')

			# Confirm the config list isn't more than two parameters
			if len(config_list) > 2: raise Exception (f'Config param error. Too many param levels: {self._original_text}')

			# Add the config assignment to the set
			self._config_assignment_set.add(tuple(config_list))

	def _parseContainer (self):

		# Check if the attribute is not a container
		if not self.is_container: raise Exception (f'Attribute assignment error. Not a container: {self.is_container} {self._original_text}')

		# Check if the attribute is a multiline text
		if self._multiline_text: raise Exception (f'Container assignment error. Cannot be multiline text: {self._original_text}')

		# Assign the container text
		container_url = self._original_text.strip().split('"')[1]
		
		# Store the container
		self._container = Singularity.fromURL(container_url, self._singularity_dir)
	
	def updateSnakeResource (self):

		# Check if the attribute is not a resource
		if not self.is_resource: raise Exception (f'Attribute assignment error. Not a resource: {self.is_resource} {self._original_text}')
		
		# Loop the resource assignments
		for resource_original_text, resource_assignment in self._resource_replacment_dict.items():

			# Replace the resource assignment with the string
			self._original_text = self._original_text.replace(resource_original_text, resource_assignment)

		return self._original_text
	
	def updateSnakeContainer (self):

		# Skip updating if the singularity dir is not provided
		if not self._singularity_dir: return self._original_text

		# Check if the attribute is not a container
		if not self.is_container: raise Exception (f'Attribute assignment error. Not a container: {self.is_container} {self._original_text}')

		# Return the container path
		return self._container.updateContainer(self._original_text)

	def updateParams (self):

		# Check if the attribute is a container, return the container dict
		if self.is_container: 
			if not self._singularity_dir: return {}
			else: return {self._container.returnPath(): self._container}

		# Check if the attribute is a resource, return the resource assignment dict
		elif self.is_resource: return self._resource_assignment_dict

		# If not a resource or container, return the config assignment set
		else: return self._config_assignment_set

def _attributeLevel (attribute_list):
	for attribute_level, attribute_item in enumerate(attribute_list):
		if not attribute_item: continue
		return attribute_level