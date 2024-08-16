#!/usr/bin/env python

import os
import sys
import random
import string
import argparse
import datetime

from pydoc import locate

from pipemake.logger import *
from pipemake.config import *
from pipemake.snakemakeIO import SnakePipelineIO

def jobRandomString (num_chars = 4):
	global random_string

	# Create a random string, if not done already
	if not random_string: random_string = ''.join(random.choices(string.digits + string.ascii_uppercase, k = num_chars))
	
	return random_string

def jobTimeStamp ():
	global time_stamp

	# Create a time stamp, if not done already
	if not time_stamp: time_stamp = datetime.datetime.now().strftime("%Y-%m-%d")

	return time_stamp

def pipeline_parser (config_parser_pipelines):

	# Change default behavior on error
	class MyParser(argparse.ArgumentParser): 
		def error(self, message):
			sys.stderr.write('error: %s\n\n' % message)
			self.print_help()
			sys.exit(2)

	# Remove metavar line in parser
	class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
		def _format_action(self, action):
			parts = super(argparse.RawDescriptionHelpFormatter, self)._format_action(action)
			if action.nargs == argparse.PARSER: parts = "\n".join(parts.split("\n")[1:])
			return parts

	def confirmDir ():
		'''Custom action to confirm directory exists'''
		class customAction(argparse.Action):
			def __call__(self, parser, args, value, option_string=None):
				if not os.path.isdir(value):
					raise IOError(f'Unable to find directory: {value}')
				setattr(args, self.dest, value)
		return customAction

	def confirmFile ():
		'''Custom action to confirm file exists'''
		class customAction(argparse.Action):
			def __call__(self, parser, args, value, option_string=None):
				if not os.path.isfile(value):
					raise IOError(f'Unable to find file: {value}')
				setattr(args, self.dest, value)
		return customAction
	
	# Create the pipeline selection parser
	pipeline_parser = MyParser(formatter_class = SubcommandHelpFormatter)
	pipeline_parser._positionals.title = "Pipeline selection argument options (positional)"
		
	# Create the subparsers
	pipeline_subparsers = pipeline_parser.add_subparsers(dest = 'pipeline', required = True)

	# Assign the arguments for each pipeline
	for pipeline_name, pipeline_params in config_parser_pipelines.items():

		# Create the subparser
		pipeline_subparser = pipeline_subparsers.add_parser(pipeline_name, help = pipeline_params['help'], add_help = False)

		pipeline_arg_groups = {}

		# Add the required argument group
		pipeline_arg_groups['required'] = pipeline_subparser.add_argument_group(f"{pipeline_name} required arguments")

		# Loop the pipeline groups
		for pipeline_arg_group in pipeline_params['arg-groups']:

			# Create the argument group, if not the basic subparsers
			if pipeline_arg_group != 'basic':

				# Check if the group already exists, then raise an exception if it does
				if pipeline_arg_group in pipeline_arg_groups: raise Exception(f'Cannot create group {pipeline_arg_group}. Group already exists.')

				# Create the argument group
				pipeline_arg_groups[pipeline_arg_group] = pipeline_subparser.add_argument_group(f"{pipeline_name} {pipeline_arg_group} arguments")

		# Add the optional argument group
		pipeline_arg_groups['optional'] = pipeline_subparser.add_argument_group(f"{pipeline_name} optional arguments")

		# Loop the pipeline mutually exclusive groups
		for pipeline_arg_group, group_args in pipeline_params['arg-groups'].items():
			if 'mutually-exclusive-groups' not in group_args: continue

			# Loop the mutually exclusive groups
			for me_group_name, me_group_args in group_args['mutually-exclusive-groups'].items():

				# Confirm the mutually exclusive group has not been created
				if me_group_name in pipeline_arg_groups: raise Exception(f'Cannot create mutually exclusive group {me_group_name}. Group already exists.')

				# Check if the group is required, if not given
				if 'required' not in me_group_args: me_group_args['required'] = False

				# Create the mutually exclusive group
				if pipeline_arg_group == 'basic' and me_group_args['required']:
					pipeline_arg_groups[me_group_name] = pipeline_arg_groups['required'].add_mutually_exclusive_group(required = me_group_args['required'])
				elif pipeline_arg_group == 'basic': 
					pipeline_arg_groups[me_group_name] = pipeline_arg_groups['optional'].add_mutually_exclusive_group(required = me_group_args['required'])
				elif pipeline_arg_group != 'basic' and pipeline_arg_group in pipeline_arg_groups:
					pipeline_arg_groups[me_group_name] = pipeline_arg_groups['optional'].add_mutually_exclusive_group(required = me_group_args['required'])
				else: raise Exception(f'Unable to assign mutually exclusive group: {me_group_name}')
		
		# Loop the pipeline arguments
		for pipeline_arg_group, group_args in pipeline_params['arg-groups'].items():

			'''2. Add the arguments to the argument subparser'''

			# Loop the arguments in the group
			for pipeline_arg_name, arg_args in group_args['args'].items():

				# Set the datatypes
				if 'type' in arg_args: 
					try: arg_args['type'] = locate(arg_args['type'])
					except: raise Exception(f"Unable to locate type: {arg_args['type']}")
	
				# Configure the action parameter, if specified
				if 'action' in arg_args:
					if arg_args['action'] == 'confirmDir': arg_args['action'] = confirmDir()
					elif arg_args['action'] == 'confirmFile': arg_args['action'] = confirmFile()

				# Configure the default params, if specified
				if 'default' in arg_args and isinstance(arg_args['default'], dict):
					default_str = arg_args['default']['str']
					if 'suffix' in arg_args['default']:
						if isinstance(arg_args['default']['suffix'], str): arg_args['default']['suffix'] = [arg_args['default']['suffix']]
						for suffix in arg_args['default']['suffix']:

							# Add underscores if needed
							if default_str: default_str += '_'

							# Process suffix strings
							if isinstance(suffix, str): default_str += suffix

							# Process suffix dicts
							elif isinstance(suffix, dict):

								# Convert the suffix dict to lowercase keys
								suffix_dict =  {_k.lower():_v for _k, _v in suffix.items()}

								if 'function' not in suffix_dict: raise Exception(f"Suffix dict not supported: {suffix_dict}")

								if suffix['function'] == 'jobTimeStamp':  default_str += jobTimeStamp()
								elif suffix['function'] == 'jobRandomString':  default_str += jobRandomString()
								else: raise Exception(f"Function not supported: {suffix['function']}")
							
					arg_args['default'] = default_str

				# Assign the argument to a mutually exclusive group, if applicable
				if 'mutually-exclusive' in arg_args:
					me_group = arg_args['mutually-exclusive']
					del arg_args['mutually-exclusive']
					pipeline_arg_groups[me_group].add_argument(f'--{pipeline_arg_name}', **arg_args)

				# Assign the argument to it respective group, if applicable
				elif pipeline_arg_group != 'basic' and pipeline_arg_group in pipeline_arg_groups: 
					pipeline_arg_groups[pipeline_arg_group].add_argument(f'--{pipeline_arg_name}', **arg_args)

				# Assign the argument to the required group, if basic and required
				elif pipeline_arg_group == 'basic' and 'required' in arg_args: 
					pipeline_arg_groups['required'].add_argument(f'--{pipeline_arg_name}', **arg_args)

				# Assign the argument to the optional group, if basic and not required
				elif pipeline_arg_group == 'basic': 
					pipeline_arg_groups['optional'].add_argument(f'--{pipeline_arg_name}', **arg_args)
				
				# Report an error if the group is not supported
				else:
					raise Exception(f'Argument group not supported: {pipeline_arg_group}')
				
		# Add the common optional arguments, but at the end of the list
		pipeline_arg_groups['optional'].add_argument('--workflow-prefix', help = 'Assign the workflow prefix', type = str, default = f'Workflow_{jobTimeStamp()}_{jobRandomString()}')
		pipeline_arg_groups['optional'].add_argument('--work-dir', help = 'Assign the working directory for snakemake. If not provided, the current directory will be used.', type = str)
		pipeline_arg_groups['optional'].add_argument('--scale-threads', help = 'Scale the threads for each task', type = float, default = 1.0)
		pipeline_arg_groups['optional'].add_argument('--scale-mem', help = 'Scale the memory (RAM) for each task', type = float, default = 1.0)
		pipeline_arg_groups['optional'].add_argument('--resource-yml', help = 'Create a seperate resource yaml', action = 'store_true')
		pipeline_arg_groups['optional'].add_argument('--singularity-dir', help = 'Assign different directory of singularity images', type = str)
		pipeline_arg_groups['optional'].add_argument('-h', '--help', action = 'help', help = 'show this help message and exit')

	return vars(pipeline_parser.parse_args())

# Create variables to store randomString and timeStamp
random_string = None
time_stamp = None

def main():

	# Assign the pipeline directory
	if os.environ.get('PM_SNAKEMAKE_DIR'): pipeline_storage_dir = os.environ.get('PM_SNAKEMAKE_DIR')
	else: pipeline_storage_dir = 'pipelines'

	# Confirm the pipeline directory exists
	if not os.path.isdir(pipeline_storage_dir): raise Exception(f'Unable to find pipeline directory: {pipeline_storage_dir}')

	# Loads the configs
	pipeline_config_args, pipeline_setup, pipeline_cmd_line, pipeline_snakefiles = loadPipelineConfigs(pipeline_storage_dir)

	# Parse the aguments from the configs
	pipeline_args = pipeline_parser(pipeline_config_args)

	# Assign the pipeline directory to an environment variable, if found
	if os.environ.get('PM_SINGULARITY_DIR') and not pipeline_args['singularity_dir']: pipeline_args['singularity_dir'] = os.environ.get('PM_SINGULARITY_DIR')
	
	# Update the pipeline args with the pipeline directory
	pipeline_args['pipeline_storage_dir'] = pipeline_storage_dir

	# Check that a pipeline was assigned
	if not pipeline_args['pipeline']: raise Exception(f'No pipeline specified')

	# Assign the pipeline job directory
	pipeline_args['pipeline_job_dir'] = os.path.join(pipeline_args['workflow_prefix'], f'pipemake')

	# Check if the pipeline job directory should be updated
	if pipeline_args['work_dir']: #and not os.path.exists(pipeline_args['work_dir']): os.makedirs(pipeline_args['work_dir'])
		pipeline_args['pipeline_job_dir'] = os.path.join(pipeline_args['work_dir'], pipeline_args['pipeline_job_dir'])

	# Create the pipeline job directory
	if not os.path.exists(pipeline_args['pipeline_job_dir']): os.makedirs(pipeline_args['pipeline_job_dir'])

	# Start logger and log the arguments
	startLogger(os.path.join(pipeline_args['pipeline_job_dir'], f'pipeline.log'))
	logArgDict(pipeline_args, omit = ['pipeline_job_dir'])

	# Process the pipeline setup
	setup_arg_dict = processPipelineSetup(pipeline_setup[pipeline_args['pipeline']], pipeline_args)

	# Update the pipeline args if the setup created new args
	if setup_arg_dict: pipeline_args.update(setup_arg_dict)

	# Create the snakemake pipeline
	snakemake_pipeline = SnakePipelineIO.open(**pipeline_args)

	# Add the snakemake modules to the pipeline
	for smkm_filename in pipeline_snakefiles[pipeline_args['pipeline']]:
		snakemake_pipeline.addModule(smkm_filename)

	# Build the singularity containers
	snakemake_pipeline.buildSingularityContainers()

	# Create the snakemake config file
	snakemake_pipeline.writeConfig(pipeline_args)

	# Create the snakemake piepline file
	snakemake_pipeline.writePipeline()

	# Close the snakemake pipeline
	snakemake_pipeline.close()

	# Create the command line
	print(processPipelineCmdLine(pipeline_cmd_line[pipeline_args['pipeline']], pipeline_args))

if __name__ == '__main__':
	main()