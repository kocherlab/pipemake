.. _installation:

##################################
Pipemake Installation Instructions
##################################

*****
Conda
*****

.. attention::

    We recommend installing Pipemake using Conda. This is the simplest and most reliable installation method to ensure that all dependencies are installed correctly.

To install Pipemake in an environment with `Snakemake <https://snakemake.readthedocs.i/>`_, please run the following commands:

.. code-block:: bash

    conda install -n base -c conda-forge mamba
    mamba create -c conda-forge -c bioconda -c kocherlab -n pipemake pipemake

.. hint::

    Please note that a version of Conda (Anaconda/Miniconda) must be installed prior to running the above snippet. We recommend Miniconda, a lightweight version of Anaconda that does not include unnecessary packages.

    Miniconda can be installed from the `Anaconda website <https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links>`_. Please follow the provided instructions for your operating system.

***
pip
***

While Pipemake may also be installed using pip, Snakemake will have limited functionality and this method is not recommended. Instead we only recommend using pip if Snakemake is already installed.

.. code-block:: bash

    pip install pipemake

.. note::

    Requires: Python >=3.11


