.. _installation:

#########################
Installation Instructions
#########################

We recommend installing pipemake using mamba. This is the simplest and most reliable installation method to ensure that all dependencies are installed correctly.

To obtain mamba, we currently recommend using Miniforge, a lightweight package that includes both conda and mamba. Miniforge may be installed from the `Miniforge GitHub page <https://github.com/conda-forge/miniforge>`_.

*****
mamba
*****

To install pipemake with all currently available pipelines you may run the following command:

.. code-block:: bash

    mamba create -c conda-forge -c bioconda -c kocherlab -n pipemake pipemake

.. note::
    
    The pipelines directory may be found within the share directory of the conda environment.

If you wish to maintain the pipelines in a separate directory, you may install the pipemake without the pipelines using the following command:

.. code-block:: bash

    mamba create -c conda-forge -c bioconda -c kocherlab -n pipemake pipemake-minimal

The pipelines directory may then be specified using an environmental variable, as described in :ref:`pipelines-directory`.

***
pip
***
.. caution::
    
    While pipemake may also be installed using pip, Snakemake will have limited functionality and this method is not recommended. Instead we only recommend using pip if an environment with Snakemake already exists.

.. code-block:: bash

    pip install pipemake

.. important::

    The pip installation does not include the pipelines. You will need to set the environmental variable `PM_SNAKEMAKE_DIR` to the location of the pipelines directory, as described in :ref:`pipelines-directory`.

.. note::

    We recommend `Snakemake <https://snakemake.readthedocs.io/>`_ 8 or greater and require python 3.10 or greater.


.. _environment-variables:

*********************
Environment variables
*********************

pipemake uses the following environmental variables to locate shared resources:

.. list-table::
    :header-rows: 1
    :widths: 25 75

    * - Variable
      - Description
    * - `PM_SNAKEMAKE_DIR`
      - Location of the pipelines directory. See :ref:`pipelines-directory`.
    * - `PM_SINGULARITY_DIR`
      - Directory used to store Singularity containers. See :ref:`singularity-containers`.

.. _pipelines-directory:

Pipelines directory
===================

The environmental variable `PM_SNAKEMAKE_DIR` may be used to specify the location of the pipelines directory. If it is not set, pipemake looks for a ``pipelines`` directory in the current working directory. Setting it is therefore necessary when pipemake is installed without the pipelines (e.g. ``pipemake-minimal`` or pip).

An ideal method is to clone the `pipemake GitHub repository <https://github.com/kocherlab/pipemake>`_ and point `PM_SNAKEMAKE_DIR` to its ``pipelines`` directory. This allows the pipelines to be updated (e.g. with ``git pull``) independently of the pipemake installation. For example:

.. code-block:: bash

    git clone https://github.com/kocherlab/pipemake.git /path/to/pipemake
    export PM_SNAKEMAKE_DIR=/path/to/pipemake/pipelines

This approach may also be used to share a single pipelines directory. A single pipelines directory may be maintained in a common location, and each user may point `PM_SNAKEMAKE_DIR` to it (e.g. in their ``.bashrc``). This ensures all users run the same version of the pipelines, and updates need only be made once.

.. _singularity-containers:

Singularity containers
======================

pipemake also includes the option to store Singularity containers in a predefined directory. This is useful for groups that wish to maintain a single set of containers. This can be done by setting the environmental variable `PM_SINGULARITY_DIR` to the desired directory. For example:

.. code-block:: bash

    export PM_SINGULARITY_DIR=/path/to/singularity

If no directory is specified, pipemake will display a warning, as this may result in redundant containers being stored for each run.

.. note::

    It's also possible to use the argument `--singularity-dir` when running the `pipemake` command to specify the desired directory. If both are provided, `--singularity-dir` takes precedence over `PM_SINGULARITY_DIR`.