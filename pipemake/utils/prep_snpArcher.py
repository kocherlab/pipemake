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

#from kocher_pipelines.logger import *
#from kocher_pipelines.config import loadPipelineConfigs, processPipelineSetup, processPipelineCmdLine, processSnakemakeModules

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

def utilParser ():

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

	util_parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, add_help = False)
	util_required = util_parser.add_argument_group('General arguments')
	util_optional = util_parser.add_argument_group('Optional general arguments')
	util_variants = util_parser.add_argument_group('Variant calling arguments')
	util_calling = util_parser.add_argument_group('Calling sites arguments')
	util_coverage = util_parser.add_argument_group('Coverage filtering arguments')
	util_qcmodule = util_parser.add_argument_group('QC module arguments')
	util_postprocessing = util_parser.add_argument_group('Postprocessing module arguments')
	util_other = util_parser.add_argument_group('Other arguments')
	
	# General arguments
	util_input = util_required.add_mutually_exclusive_group(required = True)
	util_input.add_argument('--fastq-wildcard', help = 'Wildcard statement to represent input FASTQs', type = str)
	util_input.add_argument('--fastq-table', help = 'Table with FASTQs filenames and other relevant data', type = str, action = confirmFile())
	util_required.add_argument('--assembly-fasta', help = 'Filename of the assembly file', type = str, action = confirmFile(), required = True)
	util_required.add_argument('--out-prefix', dest = 'final_prefix', help = 'Prefix of the output files', type = str, required = True)

	# Optional arguments
	util_optional.add_argument('--assembly-name', help = 'Name of the assembly. Inferred from --assembly-fasta if not specified', type = str)
	util_optional.add_argument('--species-name', help = 'Name of the species (i.e. scientific name)', type = str)
	util_optional.add_argument('--bio-project', help = 'BioProject ID. If multiple please use --fastq-table', type = str, default = '')

	# Variant calling options
	util_variants.add_argument('--minP', help = 'GATK --min-pruning. Recommendation: 1 for low coverage (<10x), 2 for high coverage (>10x)', type = int, default = 2)
	util_variants.add_argument('--minD', help = 'GATK --min-dangling-branch-length. Recommendation: 1 for low coverage (<10x), 4 for high coverage (>10x)', type = int, default = 4)

	# Site calling options
	util_calling.add_argument('--mappability-min', help = 'Minimum mappability score for genomic regions to be callable sites.', type = int, default = 1)
	util_calling.add_argument('--mappability-k', help = 'Kmer size for mappability calculations', type = int, default = 150)
	util_calling.add_argument('--mappability-merge', help = 'Max distance (in bp) to merge passing regions by mappability', type = int, default = 100)
	util_calling.add_argument('--cov-merge', help = 'Max distance (in bp) to merge passing regions by coverage', type = int, default = 100)

	# Coverage options
	util_coverage.add_argument('--cov-filter', help = 'Include coverage thresholds in the callable sites', action = 'store_true')
	util_coverage.add_argument('--cov-threshold-lower', help = 'Lower coverage threshold', type = int, default = 1)
	util_coverage.add_argument('--cov-threshold-upper', help = 'Upper coverage threshold', type = int, default = 50000)
	util_coverage.add_argument('--cov-threshold-stdev', help = 'Standard deviation threshold to consider callable', type = int)
	util_coverage.add_argument('--cov-threshold-rev', help = 'Scaling factor for coverage threshold', type = int)

	# QC module options
	util_qcmodule.add_argument('--nClusters', help = 'Number of clusters for PCA', type = int, default = 3)
	util_qcmodule.add_argument('--min-depth', help = 'Average depth threshold to exlcude samples from QC analysis', type = int, default = 2)
	
	# Postprocessing options
	util_postprocessing.add_argument('--contig-size', help = "Contig size threshold for excluding SNPs from the ‘clean’ VCF", type = int, default = 10000)
	util_postprocessing.add_argument('--maf', help = "MAF threshold for excluding SNPs from the ‘clean’ VCF", type = int, default = 0.1)
	util_postprocessing.add_argument('--missingness', help = "Missingness threshold for excluding SNPs from the ‘clean’ VCF", type = int, default = 0.75)
	util_postprocessing.add_argument('--scaffolds-to-exclude', help = "List of scaffolds/contigs to exclude from the 'clean' VCF", type = str, nargs = '+')

	# Other options
	util_other.add_argument('-h', '--help', help = 'Show this help message and exit', action = 'help')

	return vars(util_parser.parse_args())

# Create variables to store randomString and timeStamp
random_string = None
time_stamp = None

def main():

	# Parse the aguments from the configs
	util_args = utilParser()

	# Set unspecified arguments
	util_args['intervals'] = False
	util_args['sentieon'] = False
	util_args['sentieon_lic'] = ''
	util_args['remote_reads'] = False
	util_args['remote_reads_prefx'] = ''
	util_args['bigtmp'] = ''
	
	# Check if the input method uses a fastq wildcard
	if util_args['fastq_wildcard']:

		# Warn user is no species is specified
		if not util_args['species_name']:
			util_args['species_name'] = f'Species_{jobRandomString()}'
			print(f"No species name given, defaulting to: {util_args['species_name']}")

		#'bio_project': None,

	# Check if the input method uses a fastq table
	elif util_args['fastq_table']:

	else: raise Exception (f'No FASTQ input method specified')

	# Process the assembly fasta
	print(util_args['assembly_fasta'])
	print(util_args['assembly_name'])


	: None, 
	'': '../processIO.py', 
	'assembly_fasta': '../processIO.py', 
	



	'''

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

	'''

if __name__ == '__main__':
	main()