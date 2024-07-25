import os
import logging
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
		if self._url: self.processURL()

	def processURL (self):

		# Set the image name and tag
		image_name, image_tag = self._url.rsplit('/', 1)[-1].split(':')

		# Set the image name
		image_filename = f"{image_name}.{image_tag}.sif"

		# Set the image path if the singularity path is provided
		if self._singularity_path: image_filename = os.path.join(self._singularity_path, image_filename)

		# Set the image filename
		self._image_filename = image_filename

	def download (self):
		
		# Check if the image already exists
		if os.path.isfile(self._image_filename):
			logging.info(f"Image already exists at: {self._image_filename}")
			print(f"Image already exists at: {self._image_filename}")
			return self._image_filename
		
		# Check if the singularity path exists
		if self._singularity_path and not os.path.isdir(self._singularity_path): os.makedirs(self._singularity_path)

		logging.info(f"Downloading image: {self._url}")
		print(f"Downloading image: {self._url}")

		# Download the image using singularity
		singularity_process = subprocess.Popen(['singularity', 'pull', self._image_filename, self._url], stdout = subprocess.PIPE, stderr = subprocess.PIPE)

		# Get the stdout and stderr
		_, singularity_stderr = singularity_process.communicate()

		# Check for errors
		if 'FATAL' in singularity_stderr.decode(): raise Exception(singularity_stderr.decode().strip())

		# Confirm the image was downloaded
		if not os.path.isfile(self._image_filename): raise Exception('Image was not downloaded')

		logging.info(f"Image created at: {self._image_filename}")
		print(f"Image created at: {self._image_filename}")
	
	def returnPath (self):
		return self._image_filename

	def updateContainer (self, str):
		return str.replace(self._url, os.path.abspath(self._image_filename))
	
	@classmethod
	def fromURL (cls, url, singularity_path = ''):
		return cls(url = url, singularity_path = singularity_path)
