#!/usr/bin/env python

import os
import sys
import yaml
import random
import string
import logging
import argparse
import datetime

from pydoc import locate
from collections import defaultdict

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

		# Add the argument categories: required, optional, and config-assigned
		pipeline_required = pipeline_subparser.add_argument_group(f"{pipeline_name} required arguments")
		pipeline_optional = pipeline_subparser.add_argument_group(f"{pipeline_name} optional arguments")
		pipeline_arg_categories = {}

		# Add argument groups, for special cases such as mutually exclusive arguments
		pipeline_arg_groups = {}

		# Create additional categories, if needed
		for pipeline_arg_category in pipeline_params['args']:

			# Skip the default subparser, as args are stored in required and optional
			if pipeline_arg_category == 'default': continue

			# Create the subparser if it doesn't exist
			if pipeline_arg_category not in pipeline_arg_categories:
				pipeline_arg_categories[pipeline_arg_category] = pipeline_subparser.add_argument_group(f"{pipeline_name} {pipeline_arg_category}")

		# Check for pipeline groups
		if 'groups' in pipeline_params:

			# Loop the pipeline groups
			for pipeline_group, group_info in pipeline_params['groups'].items():

				# Confirm a type has been specified
				if 'type' not in group_info:
					raise Exception(f'Cannot create group {pipeline_group}. No type given.')

				# Check if a category was specified, then assign the group to the category
				if 'category' in group_info['args']:
					if group_info['args']['category'] not in pipeline_arg_categories:
						raise Exception(f'Cannot create group {pipeline_group}. Category {group_info["args"]["category"]} does not exist.')
					group_parser = pipeline_arg_categories[group_info['args']['category']]

				# Assign to required if no category was specified and the group is required
				elif 'required' in group_info['args']: group_parser = pipeline_required

				# Assign to optional if no category was specified and the group is optional
				else: group_parser = pipeline_optional

				# Check if the group type is mutually exclusive, then add the group
				if group_info['type'] == 'mutually_exclusive':
					pipeline_arg_groups[pipeline_group] = group_parser.add_mutually_exclusive_group(**group_info['args'])

				# Raise exception if group type is not supported
				else: raise Exception(f'Group type not supported: {group_info["type"]}')

		# Loop the pipeline arguments
		for pipeline_arg_category, pipeline_args in pipeline_params['args'].items():
			for pipeline_arg, arg_params in pipeline_args.items():
				
				'''
				1. Store the argument in the appropriate parser
					a) If the argument is in a special group, add it to the group
					b) If the argument is in a defined category, add it to the category
					c) If the argument is required, add it to the required parser
					d) Otherwise, add it to the optional parser
				'''
				if 'group' in arg_params:
					argument_parser = pipeline_arg_groups[arg_params['group']]
					del arg_params['group']
				elif pipeline_arg_category != 'params': 
					argument_parser = pipeline_arg_categories[pipeline_arg_category]
				elif 'required' in arg_params: argument_parser = pipeline_required
				else: argument_parser = pipeline_optional
				
				# Set the datatypes
				if 'type' in arg_params: 
					try: arg_params['type'] = locate(arg_params['type'])
					except: raise Exception(f"Unable to locate type: {arg_params['type']}")
	
				# Configure the action parameter, if specified
				if 'action' in arg_params:
					if arg_params['action'] == 'confirmDir': arg_params['action'] = confirmDir()
					elif arg_params['action'] == 'confirmFile': arg_params['action'] = confirmFile()

				# Configure the default params, if specified
				if 'default' in arg_params and not isinstance(arg_params['default'], str):
					default_str = arg_params['default']['str']
					if 'suffix' in arg_params['default']:
						if isinstance(arg_params['default']['suffix'], str): arg_params['default']['suffix'] = [arg_params['default']['suffix']]
						for suffix in arg_params['default']['suffix']:

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
							
					arg_params['default'] = default_str

				# Assign the argument to the parser
				argument_parser.add_argument(f'--{pipeline_arg}', **arg_params)

		# Add the common optional arguments, but at the end of the list
		pipeline_optional.add_argument('--scale-threads', help = 'Scale the threads for each task', type = float, default = 1.0)
		pipeline_optional.add_argument('--scale-mem', help = 'Scale the memory (RAM) for each task', type = float, default = 1.0)
		pipeline_optional.add_argument('--resource-yml', help = 'Create a seperate resource yaml', action = 'store_true')
		pipeline_optional.add_argument('--singularity-dir', help = 'Assign different directory of singularity images', type = str, default = '/Genomics/argo/users/aewebb/.local/images/')
		pipeline_optional.add_argument('-h', '--help', action = 'help', help = 'show this help message and exit')

	return vars(pipeline_parser.parse_args())

# Create variables to store randomString and timeStamp
random_string = None
time_stamp = None

def main():

	# Assign the pipeline directory
	if os.environ.get('KPDIR'): pipeline_storage_dir = os.environ.get('KPDIR')
	else: pipeline_storage_dir = 'pipelines'

	# Confirm the pipeline directory exists
	if not os.path.isdir(pipeline_storage_dir): raise Exception(f'Unable to find pipeline directory: {pipeline_storage_dir}')

	# Loads the configs
	pipeline_config_args, pipeline_setup, pipeline_cmd_line, pipeline_snakefiles = loadPipelineConfigs(pipeline_storage_dir)

	# Parse the aguments from the configs
	pipeline_args = pipeline_parser(pipeline_config_args)

	# Update the pipeline args with the pipeline directory
	pipeline_args['pipeline_storage_dir'] = pipeline_storage_dir

	# Check that a pipeline was assigned
	if not pipeline_args['pipeline']: raise Exception(f'No pipeline specified')

	# Create the working directory
	if not os.path.exists(pipeline_args['work_dir']): os.makedirs(pipeline_args['work_dir'])

	# Create the pipeline job directory
	pipeline_args['pipeline_job_dir'] = os.path.join(pipeline_args['work_dir'], f'.pipeline')
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