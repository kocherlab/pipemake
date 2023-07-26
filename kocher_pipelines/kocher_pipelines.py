#!/usr/bin/env python

import os
import sys
import yaml
import random
import string
import logging
import argparse
import datetime

from collections import defaultdict

from kocher_pipelines.logger import *
from kocher_pipelines.config import loadPipelineConfigs, processPipelineSetup, processPipelineCmdLine, processSnakemakeModules

class MyDumper(yaml.Dumper):
	def increase_indent(self, flow=False, indentless=False):
		return super(MyDumper, self).increase_indent(flow, False)

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
	#pipeline_parser.add_argument('--singularity-dir', help = 'Assign different directory of singularity images', type = str, default = '/Genomics/argo/users/aewebb/.local/images/')
	
	# Create the subparsers
	pipeline_subparsers = pipeline_parser.add_subparsers(dest = 'pipeline', required = True)

	# Assign the arguments for each pipeline
	for pipeline_name, pipeline_params in config_parser_pipelines.items():

		# Create the sub parser and groups
		pipeline_subparser = pipeline_subparsers.add_parser(pipeline_name, help = pipeline_params['help'], add_help = False)
		pipeline_required = pipeline_subparser.add_argument_group(f"{pipeline_name} required arguments")
		pipeline_optional = pipeline_subparser.add_argument_group(f"{pipeline_name} optional arguments")
		pipeline_groups = {}

		# Check for pipeline groups
		if 'groups' in pipeline_params:

			# Loop the pipeline groups
			for pipeline_group, group_info in pipeline_params['groups'].items():

				# Confirm a type has been specified
				if 'type' not in group_info:
					raise Exception(f'Cannot create group {pipeline_group}. No type given.')

				if 'required' in group_info['args']: group_parser = pipeline_required
				else: group_parser = pipeline_optional

				# Assign the group subparser
				if group_info['type'] == 'mutually_exclusive':

					pipeline_groups[pipeline_group] = group_parser.add_mutually_exclusive_group(**group_info['args'])

		# Loop the pipeline arguments
		for pipeline_arg, arg_params in pipeline_params['args'].items():

			# Assign the argument to a group
			if 'group' in arg_params:
				pipeline_group = pipeline_groups[arg_params['group']]
				del arg_params['group']
			elif 'required' in arg_params: pipeline_group = pipeline_required
			else: pipeline_group = pipeline_optional
			
			# Set the datatypes
			if 'type' in arg_params: 
				if arg_params['type'] == 'str': arg_params['type'] = str
				else: print(arg_params['type'])

			# Configure the action parameter
			if 'action' in arg_params:
				if arg_params['action'] == 'confirmDir': arg_params['action'] = confirmDir()
				elif arg_params['action'] == 'confirmFile': arg_params['action'] = confirmFile()

			# Configure the default parameter
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

			# Assign the argument
			pipeline_group.add_argument(f'--{pipeline_arg}', **arg_params)

		# Add the help argument back, but at the end of the list
		pipeline_optional.add_argument('-h', '--help', action = 'help', help = 'show this help message and exit')

	return vars(pipeline_parser.parse_args())

# Create variables to store randomString and timeStamp
random_string = None
time_stamp = None

def main():

	# Assign the pipeline directory
	if os.environ.get('KPDIR'): pipeline_dir = os.environ.get('KPDIR')
	else: pipeline_dir = 'Snakemake'

	# Confirm the pipeline directory exists
	if not os.path.isdir(pipeline_dir): raise Exception(f'Unable to find pipeline directory: {pipeline_dir}')

	# Loads the cofigs
	pipeline_config_args, pipeline_setup, pipeline_cmd_line, pipeline_snakefiles = loadPipelineConfigs(pipeline_dir)

	# Parse the aguments from the configs
	pipeline_args = pipeline_parser(pipeline_config_args)

	# Check that a pipeline was assigned
	if not pipeline_args['pipeline']: raise Exception(f'No pipeline specified')

	# Create the working directory
	if not os.path.exists(pipeline_args['work_dir']): os.makedirs(pipeline_args['work_dir'])

	# Create the pipeline config directory
	pipeline_args['pipeline_config_dir'] = os.path.join(pipeline_args['work_dir'], f'.pipeline')
	if not os.path.exists(pipeline_args['pipeline_config_dir']): os.makedirs(pipeline_args['pipeline_config_dir'])

	# Start logger and log the arguments
	startLogger(os.path.join(pipeline_args['pipeline_config_dir'], f'pipeline.log'))
	logArgDict(pipeline_args, omit = ['pipeline_config_dir'])

	# Process the pipeline setup
	setup_arg_dict = processPipelineSetup(pipeline_setup[pipeline_args['pipeline']], pipeline_args)

	# Update the pipeline args if the setup created new args
	if setup_arg_dict: pipeline_args.update(setup_arg_dict)

	# Process the snakemake modules
	snakemake_config = processSnakemakeModules(pipeline_snakefiles[pipeline_args['pipeline']], pipeline_dir, pipeline_args)

	# Build YAML
	yaml_dict = {'workdir': pipeline_args['work_dir']}

	for group, config_list in snakemake_config.items():
		if group == 'root':
			for config in config_list:
				yaml_dict[config] = pipeline_args[config.replace('-', '_')]
		else:
			group_dict = {}
			for config in config_list:
				group_dict[config] = pipeline_args[config.replace('-', '_')]
			yaml_dict[group] = group_dict

	# Create the YAML file
	yaml_file = open(f"{pipeline_args['snakemake_job_prefix']}.yml", 'w')
	yaml_file.write(yaml.dump(yaml_dict, Dumper = MyDumper, sort_keys = False))
	yaml_file.close()

	# Create the command line
	print(processPipelineCmdLine(pipeline_cmd_line[pipeline_args['pipeline']], pipeline_args))

if __name__ == '__main__':
	main()