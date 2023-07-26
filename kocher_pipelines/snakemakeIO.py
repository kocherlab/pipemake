import os
from collections import defaultdict

def parseSnakemakefile (filename):

	# Create list the config and input items
	config_items = defaultdict(list)
	input_items = []

	# Set the block bools to false
	config_block = False
	rule_all_block = False
	rule_found = False
	input_found = False

	# Parse the snakemake file
	with open(filename, 'r') as snakemake_file:
		for snakemake_line in snakemake_file:

			# Confirm the config comment block
			if snakemake_line.startswith('"""config'): 
				config_block = True
				continue
			elif snakemake_line.startswith('config"""'): 
				config_block = False
				continue

			# Parse the config block
			if config_block:
				if ':' not in snakemake_line: config_items['root'].append(snakemake_line.strip())
				else:
					if snakemake_line.count(':') > 1: raise Exception(f'Unable to parse config line: {snakemake_line.strip()}')
					config_group, config_item = (s_i.strip() for s_i in snakemake_line.split(':'))
					config_items[config_group].append(config_item)
			
			# Confirm the rule all comment block
			if snakemake_line.startswith('"""rule-all'): 
				rule_all_block = True
				continue
			elif snakemake_line.startswith('rule-all"""'): 
				rule_all_block = False
				rule_found = False
				input_found = False
				continue

			# Confirm the rule all block is correctly formatted, and store the input
			if rule_all_block and 'rule all:' in snakemake_line: rule_found = True
			elif rule_all_block and rule_found and 'input:' in snakemake_line: input_found = True
			elif rule_all_block and rule_found and input_found: input_items.append(snakemake_line.strip())

	if not input_items: 
		raise Exception(f'Unable to assgin input from: {filename}. Please check the input block is correctly formatted')

	return config_items, input_items

def copySnakemakeFile (filename, out_prefix = '', out_filename = '', out_dir = '', keep_blocks = False):

	# Set the rule block to false
	config_block = False
	rule_all_block = False

	if not out_prefix and not out_filename:
		raise Exception(f'No output method given for snakefile')

	# Create the output file
	if out_prefix: out_path = f'{out_prefix}.smk'
	else: out_path = out_filename

	# Update the path with the output directory, if given
	if out_dir: out_path = os.path.join(out_dir, out_path)

	# Create the output file
	snakemake_output_file = open(out_path, 'w')

	# If keeping the blocks, store the config filename and workdir
	if keep_blocks:

		if out_filename: raise Exception(f'out_filename cannnot be used alongside keep_blocks')

		# Include the config filename
		snakemake_output_file.write(f"configfile: '{out_prefix}.yml'\n\n")

		# Assign the work directory
		snakemake_output_file.write("workdir: config['workdir']\n\n")

	# Parse the snakemake file
	with open(filename, 'r') as snakemake_file:
		for snakemake_line in snakemake_file:
			
			# Confirm the config comment block
			if snakemake_line.startswith('"""config'): 
				config_block = True
				continue
			elif snakemake_line.startswith('config"""'): 
				config_block = False
				if not keep_blocks: next(snakemake_file)
				continue

			# Cofirm the rule all comment block
			if snakemake_line.startswith('"""rule-all'): 
				rule_all_block = True
				continue
			elif snakemake_line.startswith('rule-all"""'): 
				rule_all_block = False
				if not keep_blocks: next(snakemake_file)
				continue

			# Write the contents of the file
			if config_block and keep_blocks:
				config_entry = snakemake_line.strip()
				if ':' not in config_entry: snakemake_output_file.write(f"{config_entry} = config['{config_entry}']\n")
				else:
					config_key, config_item = (_i.strip() for _i in config_entry.split(':'))
					snakemake_output_file.write(f"{config_item} = config['{config_key}']['{config_item}']\n")
			elif rule_all_block and keep_blocks: snakemake_output_file.write(snakemake_line)
			elif not config_block and not rule_all_block: snakemake_output_file.write(snakemake_line)

def createSnakemakeFile (job_prefix, input_files, module_dir, modules_to_include):

	# Check the number of input files
	len_input_files = len(input_files)

	snakemake_output_file = open(f'{job_prefix}.smk', 'w')

	# Include the config filename
	snakemake_output_file.write(f"configfile: '{job_prefix}.yml'\n\n")

	# Assign the work directory
	snakemake_output_file.write("workdir: config['workdir']\n\n")

	# Assign rule all:
	snakemake_output_file.write('rule all:\n')
	snakemake_output_file.write('\tinput:\n')
	
	# Assign the input files within rule all:
	for input_pos in range(len_input_files):
		snakemake_output_file.write(f'\t\t{input_files[input_pos]}')
		if input_pos + 1 == len_input_files: snakemake_output_file.write('\n')
		else: snakemake_output_file.write(',\n')

	# Create blank line between input and modules
	snakemake_output_file.write('\n')
	
	# Create the input block of modules
	for module_to_include in modules_to_include:
		snakemake_output_file.write(f'include: "{os.path.join(module_dir, module_to_include)}"\n')
