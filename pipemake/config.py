import os
import copy
import yaml

from collections import defaultdict

from pipemake.processIO import standardizeInput, returnSamples, returnPaths
#from pipemake.snakemakeIO import SnakePipelineIO

def loadPipelineConfigs (directory):

	# Create dicts to store the relevant yaml data
	pipeline_config_args = {}
	pipeline_setup = {}
	pipeline_cmd_line = {}
	pipeline_snakefiles = defaultdict(list)
	
	# Check the config directory exists
	config_dir = os.path.join(directory, 'configs')
	if not os.path.isdir(config_dir): raise Exception(f'Unable to find pipeline directory: {config_dir}')

	# Loop the configs
	for config_filename in os.listdir(config_dir):
		config_file = os.path.join(config_dir, config_filename)

		# Confirm the path is a file
		if not os.path.isfile(config_file): continue

		with open(config_file, "r") as config_stream:
			try: config_yaml = yaml.safe_load(config_stream)
			except: raise Exception('Error opening YAML')

			# Check pipeline_names for repeats, and if found report error
			if config_yaml['pipeline'] in pipeline_config_args:
				raise Exception(f"Pipeline {config_yaml['pipeline']} has already been assigned. Please check pipeline configs")

			pipeline_config_args[config_yaml['pipeline']] = config_yaml['parser']
			pipeline_setup[config_yaml['pipeline']] = config_yaml['setup']
			pipeline_cmd_line[config_yaml['pipeline']] = config_yaml['command-line']
			pipeline_snakefiles[config_yaml['pipeline']] = config_yaml['snakefiles']
		
	return pipeline_config_args, pipeline_setup, pipeline_cmd_line, pipeline_snakefiles

def processPipelineSetup (pipeline_setup, pipeline_args):

	# Create list to store setup args
	process_dict = {}

	for setup_name, setup_methods in pipeline_setup.items():
		for method_args in setup_methods.values():
			
			# Assign the arg groups
			input_args = method_args['input']
			
			# Check for missing arguments
			for input_arg in input_args['args']:
				if input_arg.replace('-', '_') not in pipeline_args: raise Exception(f'Setup argument {input_arg} not found among pipeline argument')

			# Confirm expected args were specified
			method_missing_args = [_a for _a in input_args['args'] if not pipeline_args[_a.replace('-', '_')]]

			# Skip if missing arguments
			if method_missing_args: continue

			if 'standardize' in method_args:

				# Create a dict of the standardize args
				std_args = copy.deepcopy(method_args['standardize'])

				# Check if the work_dir is specified
				if 'work_dir' in pipeline_args and pipeline_args['work_dir']:
					std_args['args']['work_dir'] = pipeline_args['work_dir']

				# Update the arguments
				for std_arg, arg_params in std_args['args'].items():
					if not isinstance(arg_params, str): continue
					std_args['args'][std_arg] = arg_params.replace('-', '_').format(**pipeline_args)

				# Standardize the input
				standardizeInput(**std_args)

				# Assign the method paths
				method_paths = returnPaths(**std_args)

				# Check for method paths
				if len(method_paths) > 0: 
					
					# Check if the args have already been created
					if 'singularity-args' not in process_dict:
						process_dict['singularity-args'] = defaultdict(list)

					# Return any paths that need to be accounted for
					process_dict['singularity-args']['bind'].extend(method_paths)

			if 'samples' in method_args:

				# Create a dict of the standardize args
				samples_args = copy.deepcopy(method_args['samples'])

				# Update the arguments
				for samples_arg, arg_params in samples_args['args'].items():
					if not isinstance(arg_params, str): continue
					samples_args['args'][samples_arg] = arg_params.replace('-', '_').format(**pipeline_args)

				# Assign the samples from the method
				method_samples = returnSamples(**samples_args)

				# Confirm the samples are not already assigned
				if method_samples and 'samples' in process_dict: raise Exception(f'Samples already assigned')

				# Standardize the input
				process_dict['samples'] = method_samples

	return process_dict

def processPipelineCmdLine (pipeline_cmd_line, pipeline_args):

	# Create list to store the command line arguments
	cmd_line_list = ['snakemake']

	# Process the command line
	for cmd_arg, cmd_value in pipeline_cmd_line.items():

		# Add the basic args
		if isinstance(cmd_value, bool): cmd_line_list.append(f'--{cmd_arg}')
		else: cmd_line_list.extend([f'--{cmd_arg}', cmd_value])
		
		# Check if using singularity
		if cmd_arg == 'use-singularity' and cmd_value and 'singularity-args' in pipeline_args:

			# Assign the singularity args
			singularity_args_list = []
			for singularity_arg, singularity_value in pipeline_args['singularity-args'].items():
				singularity_args_list.append(f'--{singularity_arg} ' + ','.join(singularity_value))

			# Add the singularity args
			singularity_args_str = ','.join(singularity_args_list)
			cmd_line_list.extend(['--singularity-args', f'"{singularity_args_str}"'])

	return ' '.join(map(str,cmd_line_list))
