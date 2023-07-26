import io
import os

from setuptools import setup

import kocher_pipelines

with io.open("README.rst", "rt", encoding="utf8") as f:
    readme = f.read()

# Read requirements.txt
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = lib_folder + "/requirements.txt"
requirements = []
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        requirements = f.read().splitlines()

setup(
    name=kocher_pipelines.__name__,
    version=kocher_pipelines.__version__,
    project_urls={
        "Documentation": kocher_pipelines.__docs__,
        "Code": kocher_pipelines.__code__,
        "Issue tracker": kocher_pipelines.__issue__,
    },
    license=kocher_pipelines.__license__,
    url=kocher_pipelines.__url__,
    description=kocher_pipelines.__summary__,
    long_description_content_type="text/x-rst",
    long_description=readme,
    packages=["kocher_pipelines"],
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "kocher-pipelines=kocher_pipelines.kocher_pipelines:main",
        ],
    },
    python_requires=">=3.7",
)
