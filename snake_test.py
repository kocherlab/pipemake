import os
import re
import sys
import logging

from collections import defaultdict

class SnakeFileIO ():
	def __init__ (self, smk_filename, **kwargs):

		# Confirm the file exists
		if not os.path.isfile(smk_filename):
			raise IOError (f'Unable to open: {smk_filename}. Please confirm the file exists.')

		# Assign the basic arguments
		self.filename = smk_filename
		self._indent_style = None
		self._exclude_rules = ['all']
		self._file_rules = []
		
		# Assign the indent style
		self._assignIndent()

		# Parse the snakefile
		self._parseFile()

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

				if smk_list[0] and not (smk_list[0].startswith('rule') and smk_line.rstrip().endswith(':')):
					within_rule_block = False

				# Check if the line is the start of a rule
				if smk_list[0].startswith('rule') and smk_line.rstrip().endswith(':'):

					within_rule_block = True

					# Check if a previous rule block was found
					if rule_block:
						#parseRule(rule_block)
						rule_block = ''
					
				if within_rule_block and smk_line.strip(): rule_block += smk_line

			# Parse the final rule block
			if rule_block: SnakeRuleIO.read(rule_block, indent_style=self._indent_style)
			
	@classmethod
	def open (cls, *args, **kwargs):
		return cls(*args, **kwargs)

class SnakeRuleIO ():
	def __init__ (self, rule_list, indent_style = None, **kwargs):

		# Assign the rule argument to populate
		self.rule_name = rule_list[0].split()[1][: -1]
		self._indent_style = indent_style
		self._rule_attributes = defaultdict(list)
		self._rule_config_params = defaultdict(list)
		self._rule_resource_params = defaultdict(lambda: defaultdict(str))
		
		# Assign the fixed arguments
		self._resource_attributes = ['threads', 'resources']

		self._parseRule(rule_list[1:])

	def _parseRule (self, rule_list):

		def processAttribute (rule_attribute, rule_value_str):

			if rule_value_str.endswith(','): rule_value_str = rule_value_str[:-1]

			# Check if the rule attribute is a resource attribute
			if rule_attribute in self._resource_attributes:

				if rule_attribute == 'threads': rule_resource = 'threads'
				elif rule_attribute == 'resources':

					if self._attributeLevel(rule_value_str.split(self._indent_style)) != 2: 
						raise Exception (f'Config resource error: {rule_value_str}')
					
					# Assign the resource name and value
					rule_resource, rule_value_str = [_v.strip() for _v in rule_value_str.strip().split('=', 1)]

				# Check for a simple integer assignment
				try:
					int(rule_value_str)
					self._rule_resource_params[self.rule_name][rule_resource] = int(rule_value_str)
					self._rule_attributes[rule_attribute].append(f"{rule_resource}=config['resource']['{self.rule_name}']['{rule_resource}']")
				
				except:
					
					# Find all occurences of the word resources followed by one or more sets of square brackets ending with a colon
					for rule_resource_match in re.finditer(r'resources(\[(.*?)\])+\:', rule_value_str):

						# Get the resource str
						rule_resource_str = rule_resource_match.group(0)

						# Split the param str by either square brackets using regex
						rule_resource_list = [_p for _p in re.split(r'\[|\]', rule_resource_str[:-1].replace('"','').replace("'",'')) if _p][1:]
						
						# Check that the resource list is a single argument
						if len(rule_resource_list) != 1: raise Exception (f'Config resource error: {rule_resource_list}')

						# Assign the resource argument
						rule_resource_arg = rule_resource_list[0]

						# Get text after rule_resourcem_str in rule_value_str
						rule_resource_substr = rule_value_str[rule_value_str.find(rule_resource_str) + len(rule_resource_str):]

						# Get the match between curly brackets without returning the brackets
						rule_resource_value_match = re.search(r' ?{(.*?)}', rule_resource_substr)

						# Assign the match str and value
						rule_resource_value_str = rule_resource_value_match.group(0)
						rule_resource_value = rule_resource_value_match.group(1)

						# Assign the rule resource category, argument, and value
						self._rule_resource_params[self.rule_name][rule_resource_arg] = int(rule_resource_value)

						# Create string for resource assignment replacement
						rule_resource_assignment_str = f"config['resource']['{self.rule_name}']['{rule_resource_arg}']"

						# Convert the the text
						rule_value_str = rule_value_str.replace(rule_resource_str + rule_resource_value_str, rule_resource_assignment_str)

						# Replace the resource assignment with the string
						self._rule_attributes[rule_attribute].append(f'{rule_resource}={rule_value_str}')
						
						# Assign the rule resource category, argument, and value
						self._rule_resource_params[self.rule_name][rule_resource] = rule_resource_value
				
				return

			# Find all occurences of the word config follpwed by one or more sets of square brackets
			for rule_param in re.finditer(r'config(\[(.*?)\])+', rule_value_str):

				# Get the param str
				rule_param_str = rule_param.group(0)

				# Split the param str by either square brackets using regex
				rule_param_list = [_p for _p in re.split(r'\[|\]', rule_param_str.replace('"','').replace("'",'')) if _p][1:]

				# Check if the param has a value by position
				rule_param_w_value_pos = rule_value_str.find(f'{rule_param_str}:')

				# Check if the param does not have a value, i.e. config argument
				if rule_param_w_value_pos == -1:
					if len(rule_param_list) != 2: raise Exception (f'Config param error: {rule_param_str}')

					# Assign the rule config category and argument
					rule_config_category, rule_config_arg = rule_param_list
					self._rule_config_params[rule_config_category].append(rule_config_arg)

					print(rule_config_category, rule_config_arg)
				
				# Check if the param has a value, i.e. resource argument
				else:

					if len(rule_param_list) != 3: raise Exception (f'Config param error: {rule_param_str}')

					# Assign the rule config category and argument
					_, rule_config_category, rule_config_arg = rule_param_list

					# Get text after rule_param_str in rule_value_str
					rule_param_substr = rule_value_str[rule_param_w_value_pos + len(f'{rule_param_str}:'):]

					# Get the match between curly brackets without returning the brackets
					rule_value_match = re.search(r'{(.*?)}', rule_param_substr)

					# Assign the match str and value
					rule_value_str = rule_value_match.group(0)
					rule_value = rule_value_match.group(1)

					# Assign the rule resource category, argument, and value
					self._rule_resource_params[rule_config_category][rule_config_arg] = rule_value

					print(rule_config_category, rule_config_arg, rule_value)
				
			
		# Assign the current rule attribute
		rule_attribute = ''

		# Loop the rule block and parse the attributes
		for rule_line in rule_list:
			if not rule_line.strip(): continue

			# Split the line by the indent style
			rule_line_list = rule_line.rstrip().split(self._indent_style)

			# Assign the attribute level
			attribute_level = self._attributeLevel(rule_line_list)

			# Check if the line is a rule attribute (without an assigned value)
			if attribute_level == 1 and rule_line.rstrip().endswith(':'):
				rule_attribute = rule_line_list[1][:-1]
			
			# Check if the line is a rule attribute (with an assigned value)
			elif attribute_level == 1 and not rule_line.rstrip().endswith(':'):
				
				# Split the attribute, confirm a value is assigned
				split_attribute = [_a.strip() for _a in rule_line_list[1].split(':', 1)]
				if len(split_attribute) != 2: raise Exception (f'Rule attribute error: {rule_line_list[1]}')
				processAttribute(*split_attribute)

			# Check if the line is a atrribute value line
			elif attribute_level == 2:

				# Assign the value to the current rule attribute
				#self._rule_attributes[rule_attribute] += rule_line

				processAttribute(rule_attribute, rule_line)

			# Raise an error if the attribute level is not 1 or 2
			else: raise Exception (f'Attribute level error: {rule_line}')

		'''	
				
		# Loop the rule attribute values
		for rule_values_str in self._rule_attributes.values():

			# Find all occurences of the word config follwed by one or more sets of square brackets
			for rule_param in re.finditer(r'config(\[(.*?)\])+', rule_values_str):

				# Get the param str
				rule_param_str = rule_param.group(0)

				# Split the param str by either square brackets using regex
				rule_param_list = [_p for _p in re.split(r'\[|\]', rule_param_str.replace('"','').replace("'",'')) if _p][1:]

				# Check if the param has a value by position
				rule_param_w_value_pos = rule_values_str.find(f'{rule_param_str}:')

				# Check if the param does not have a value, i.e. config argument
				if rule_param_w_value_pos == -1:
					if len(rule_param_list) != 2: raise Exception (f'Config param error: {rule_param_str}')

					# Assign the rule config category and argument
					rule_config_category, rule_config_arg = rule_param_list
					self._rule_config_params[rule_config_category].append(rule_config_arg)
				
				# Check if the param has a value, i.e. resource argument
				else:

					if len(rule_param_list) != 3: raise Exception (f'Config param error: {rule_param_str}')

					# Assign the rule config category and argument
					_, rule_config_category, rule_config_arg = rule_param_list

					# Get text after rule_param_str in rule_values_str
					rule_param_substr = rule_values_str[rule_param_w_value_pos + len(f'{rule_param_str}:'):]

					# Get the match between curly brackets without returning the brackets
					rule_value_match = re.search(r'{(.*?)}', rule_param_substr)

					# Assign the match str and value
					rule_value_str = rule_value_match.group(0)
					rule_value = rule_value_match.group(1)

					# Assign the rule resource category, argument, and value
					self._rule_resource_params[rule_config_category][rule_config_arg] = rule_value

					print(rule_config_category, rule_config_arg, rule_value)
		'''
					




		'''

				# Check if there is a next character
				if rule_config_param_end == len(rule_values_str): self._rule_config_params_wo_values[rule_config_param_str] 
				
				# Check if the next character is a colon
				if rule_values_str[rule_config_param_end] == ':':
					self._rule_config_params_wo_values[rule_config_param] = None
				else:
					self._rule_config_params_w_values[rule_config_param] = rule_values_str[rule_config_param_end: ]
				print(len(rule_values_str), rule_config_param_end, rule_values_str)

		'''
		#print(rule_values_str[rule_config_param_end])

		print('--END--')

	@classmethod
	def read (cls, rule_str, *args, **kwargs):
		return cls(rule_str.splitlines(), *args, **kwargs)
	
	@staticmethod
	def _attributeLevel (attribute_list):
		for attribute_level, attribute_item in enumerate(attribute_list):
			if not attribute_item: continue
			return attribute_level
	



SnakeFileIO.open(sys.argv[1])


pipeline_module_block = False
dict_w_defaults = False
block_name = 'rule'

def parseRule (rule_block_str):

	# Split the rule block into a list
	rule_block_list = rule_block_str.splitlines()

	# Create a dict to store the rule attributes
	rule_attributes = defaultdict(str)

	# Assign the rule name
	rule_name = rule_block_list[0].split()[1][: -1]

	# Assign the current rule attribute
	rule_attribute = ''

	# Loop the rule block and parse the attributes
	for rule_line in rule_block_list[1:]:
		if not rule_line.strip(): continue

		# Split the line by the indent style
		rule_line_list = rule_line.rstrip().split('\t')

		# Assign the attribute level
		attribute_level = attributeLevel(rule_line_list)

		# Check if the line is a rule attribute (without an assigned value)
		if attribute_level == 1 and rule_line.rstrip().endswith(':'):
			rule_attribute = rule_line_list[1][:-1]
		
		# Check if the line is a rule attribute (with an assigned value)
		elif attribute_level == 1 and not rule_line.rstrip().endswith(':'):
			
			# Split the attribute, confirm a value is assigned
			split_attribute = [_a.strip() for _a in rule_line_list[1].split(':')]
			if len(split_attribute) != 2: raise Exception (f'Rule attribute error: {rule_line_list[1]}')
			rule_attributes[split_attribute[0]] = split_attribute[1]

		# Check if the line is a atrribute value line
		elif attribute_level == 2:

			# Assign the value to the current rule attribute
			rule_attributes[rule_attribute] += rule_line

		# Raise an error if the attribute level is not 1 or 2
		else: raise Exception (f'Attribute level error: {rule_line}')
			
	# Loop the rule attribute values
	for rule_values_str in rule_attributes.values():

		# Find all occurences of the word config follwed by one or more sets of square brackets
		rule_values_config = re.finditer(r'config(\[(.*?)\])+', rule_values_str)



		#rule_values_config = re.findall(r'config\[(.*?)\]', rule_values_str)
		

		#rule_values_config = re.findall(r'config(\[(.*?)\])+', rule_values_str)
		for i in rule_values_config:
			print(i.group(0))
		#print(rule_values_str)


		
	#	print(rule_name, rule_attribute, rule_value)



# Create the param dict
if not dict_w_defaults: param_dict = defaultdict(list)
else: param_dict = defaultdict(lambda: defaultdict(str))

sys.exit()

# Create bool to check if in the config block
within_rule_block = False

rule_block = ''

# Parse the snakefile
with open(filename) as smk_file:
	for smk_line in smk_file:

		# Skip blank lines
		if not smk_line.strip(): continue

		# Split the line by the indent style (rule, attribute, param)
		smk_list = smk_line.rstrip().split('\t')

		if smk_list[0] and not (smk_list[0].startswith('rule') and smk_line.rstrip().endswith(':')):
			within_rule_block = False

		# Check if the line is the start of a rule
		if smk_list[0].startswith('rule') and smk_line.rstrip().endswith(':'):

			within_rule_block = True

			# Check if a previous rule block was found
			if rule_block:
				parseRule(rule_block)
				rule_block = ''
			
		if within_rule_block and smk_line.strip(): rule_block += smk_line

	# Parse the final rule block
	if rule_block: parseRule(rule_block)

	'''
		continue

		# Check if within a rule
		if smk_list[0]:
			if not pipeline_module_block and block_name in smk_list[0]: within_rule_block = True
			elif pipeline_module_block and smk_list[0].startswith('module') and block_name in smk_list[0]: within_rule_block = True
			elif within_rule_block: break
			else: within_rule_block = False
			continue

		# Check if within the config block
		if within_rule_block:

			#print(within_rule_block, smk_list)

			continue

			# Loop the snakefile line by level
			for smk_level, smk_item in enumerate(smk_list):
				if not smk_item: continue

				print(smk_level, smk_item)

				# Assign the param attribute
				if smk_item.endswith(':'):
					param_attribute_name = smk_item[:-1].strip()
					param_attribute_level = smk_level

				# Assign the param
				else:
					if (smk_level - param_attribute_level) != 1:
						#print(smk_level, smk_item, param_attribute_level)
						raise Exception (f'Parameter assignment error')

					# Assign the default value, if possible
					if ':' not in smk_item or '://' in smk_item: smk_value = None
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
	'''