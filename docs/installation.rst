.. _installation:

#########################
Installation Instructions
#########################

We recommend installing pipemake using mamba. This is the simplest and most reliable installation method to ensure that all dependencies are installed correctly.

To obtain mamba, we currently recommend using Miniforge, a lightweight package that includes both conda and mamba. Miniforge may be installed from the `Miniforge GitHub page <https://github.com/conda-forge/miniforge>`_.

*****
mamba
*****

To install pipemake with all currently available pipelines, you may run the following command:

.. code-block:: bash

    mamba create -c conda-forge -c bioconda -c kocherlab -n pipemake pipemake

.. note::
    
    The pipelines directory may be found within the share directory of the conda environment.

If you wish to maintain the pipelines in a separate directory, you may install the pipemake without the pipelines using the following command:

.. code-block:: bash

    mamba create -c conda-forge -c bioconda -c kocherlab -n pipemake pipemake-minimal

You may then used the environmental variable `PM_SNAKEMAKE_DIR` to specify the location of the pipelines directory. For example:

.. code-block:: bash

    export PM_SNAKEMAKE_DIR=/path/to/pipelines

***
pip
***
.. caution::
    
    While pipemake may also be installed using pip, Snakemake will have limited functionality and this method is not recommended. Instead we only recommend using pip if an environment with Snakemake already exists.

.. code-block:: bash

    pip install pipemake

.. note::

    We recommend Snakemake 8 or greater and require python 3.8 or greater.

    `Snakemake <https://snakemake.readthedocs.i/>`_,


