import logging
import sys

def startLogger (log_filename = None, filemode = 'w'):

	# Close any old loggers
	for handler in logging.root.handlers[:]:
		handler.close()
		logging.root.removeHandler(handler)

	# Star the logger
	if log_filename: logging.basicConfig(filename = log_filename, filemode = filemode, level = 'INFO', format = '%(asctime)s - %(levelname)s: %(message)s', datefmt = '%Y-%m-%d %H:%M:%S')
	else: logging.basicConfig(stream = sys.stdout, level = 'INFO', format = '%(asctime)s - %(levelname)s: %(message)s', datefmt = '%Y-%m-%d %H:%M:%S')

	# Start logging to stdout
	stdout_log = logging.StreamHandler()

	# Assign the stdout logging level
	stdout_log.setLevel(logging.WARNING)

	# Define the stdout format
	console_format = logging.Formatter('%(funcName)s - %(levelname)s: %(message)s')

	# Assign the format
	stdout_log.setFormatter(console_format)

	# Add the stdout handler
	logging.getLogger('').addHandler(stdout_log)

	# Update the exception handler to log the exception
	def expHandler(etype,val,tb):

		# Log the error
		logging.error("%s" % (val), exc_info=(etype,val,tb))

	# Update the exception hook
	sys.excepthook = expHandler

def logArgDict (arg_dict, print_undefined = False, omit = []):

	# Loop the arguments
	for arg, value in arg_dict.items():

		# Skip arg if in omit list
		if arg in omit: continue

		# Report only defined arguments, unless print_undefined is True
		if value is None and not print_undefined: continue

		# Log the argument
		logging.info('Argument %s: %s' % (arg, value))