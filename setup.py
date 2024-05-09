import io
import os

from setuptools import setup

import pipemake

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
    name=pipemake.__name__,
    version=pipemake.__version__,
    project_urls={
        "Documentation": pipemake.__docs__,
        "Code": pipemake.__code__,
        "Issue tracker": pipemake.__issue__,
    },
    license=pipemake.__license__,
    url=pipemake.__url__,
    description=pipemake.__summary__,
    long_description_content_type="text/x-rst",
    long_description=readme,
    packages=["pipemake"],
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "pipemake=pipemake.pipemake:main",
        ],
    },
    python_requires=">=3.7",
)
