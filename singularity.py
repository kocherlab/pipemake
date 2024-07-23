import os
import subprocess

def checkSingularity ():
	try: subprocess.run(['singularity', '--version'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	except FileNotFoundError: raise FileNotFoundError('Singularity is not installed on this system')

class Singularity:
	def __init__ (self, url = '', singularity_path = ''):

		# Check if singularity is installed
		checkSingularity()

		# Set the singularity attributes
		self._url = url
		self._singularity_path = singularity_path
		self._image_filename = ''

		# If the url is provided, process it
		if self._url: self._image_filename = self.processURL()

	def processURL (self):

		# Set the image name and tag
		image_name, image_tag = self._url.rsplit('/', 1)[-1].split(':')

		# Set the image name
		image_filename = f"{image_name}.{image_tag}.sif"

		# Set the image path if the singularity path is provided
		if self._singularity_path: image_filename = os.path.join(self._singularity_path, image_filename)

		# Check if the image already exists
		if os.path.isfile(image_filename): return image_filename

		# Download the image using singularity
		singularity_process = subprocess.Popen(['singularity', 'pull', image_filename, self._url], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

		# Get the stdout and stderr
		_, singularity_stderr = singularity_process.communicate()

		# Check for errors
		if 'FATAL' in singularity_stderr.decode(): raise Exception(singularity_stderr.decode().strip())

		# Confirm the image was downloaded
		if not os.path.isfile(image_filename): raise Exception('Image was not downloaded')

		# Return the image filename
		return image_filename

	def returnPath (self):
		if not self._image_filename: raise Exception('Image file not found')
		return os.path.abspath(self._image_filename)
	
	@classmethod
	def fromURL (cls, url, singularity_path = ''):
		return cls(url = url, singularity_path = singularity_path)
	


test = Singularity.fromURL('docker://aewebb/pipemake_utils:v0.1.27', 'toy_files')
print(test.returnPath())