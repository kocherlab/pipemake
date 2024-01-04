import os
import re
import sys

from collections import defaultdict

pipeline_module_block = False
dict_w_defaults = False
block_name = 'rule'

def parseRule (rule_block_str):

	def attributeLevel (attribute_list):
		for attribute_level, attribute_item in enumerate(attribute_list):
			if not attribute_item: continue
			return attribute_level

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
		print(rule_values_str)


		
	#	print(rule_name, rule_attribute, rule_value)



# Create the param dict
if not dict_w_defaults: param_dict = defaultdict(list)
else: param_dict = defaultdict(lambda: defaultdict(str))

filename = sys.argv[1]

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