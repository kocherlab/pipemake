.. only:: html
    
    |Stable version| |Documentation| |github ci| |Coverage| |conda| |Conda Upload| |PyPI Upload| |LICENSE|

    .. |Stable version| image:: https://img.shields.io/github/v/release/kocherlab/pipemake?label=stable
        :target: https://github.com/kocherlab/pipemake/releases/
        :alt: Stable version

    .. |Documentation| image::
        https://readthedocs.org/projects/pipemake/badge/?version=latest
        :target: https://pipemake.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

    .. |github ci| image::
        https://github.com/kocherlab/pipemake/actions/workflows/ci.yml/badge.svg?branch=main
        :target: https://github.com/kocherlab/pipemake/actions/workflows/ci.yml
        :alt: Continuous integration status

    .. |Coverage| image::
        https://codecov.io/gh/kocherlab/pipemake/branch/main/graph/badge.svg
        :target: https://codecov.io/gh/kocherlab/pipemake
        :alt: Coverage

    .. |conda| image::
        https://anaconda.org/kocherlab/pipemake/badges/version.svg
        :target: https://anaconda.org/kocherlab/pipemake

    .. |Conda Upload| image::
        https://github.com/kocherlab/pipemake/actions/workflows/upload_conda.yml/badge.svg
        :target: https://github.com/kocherlab/pipemake/actions/workflows/upload_conda.yml

    .. |PyPI Upload| image::
        https://github.com/kocherlab/pipemake/actions/workflows/python-publish.yml/badge.svg
        :target: https://github.com/kocherlab/pipemake/actions/workflows/python-publish.yml

    .. |LICENSE| image::
        https://anaconda.org/kocherlab/pipemake/badges/license.svg
        :target: https://github.com/kocherlab/pipemake/blob/main/LICENSE

********
pipemake
********
pipemake is a lightweight, flexible, and easy-to-use tool for creating and managing `Snakemake <https://snakemake.readthedocs.i/>`_ pipelines. It was designed with three primary goals: 

1. Offer a collection of curated, customizable genomic analysis pipelines for researchers seeking to rapidly integrate Snakemake-based workflows into their research.
2. Optimize computational efficiency and reproducibility by fully operating in the Snakemake ecosystem.
3. Streamline development by creating a flexible platform with swappable pipelines that easily reuse previously written Snakemake code.

.. toctree::
    :hidden:
    :caption: pipemake Documentation

    installation
    pipelines
    usage
    filetypes
    contact

================
Getting pipemake
================
------------
Conda
------------

.. code-block:: bash

    conda install -n base -c conda-forge mamba
    mamba create -c conda-forge -c bioconda -c kocherlab -n pipemake pipemake


======
Issues
======

1. Check the `docs <https://pipemake.rtfd.io/>`_.
2. Search the `issues on GitHub <https://github.com/kocherlab/pipemake/issues>`_ or open a new one.

============
Contributors
============

* **Andrew Webb**, Department of Ecology and Evolutionary Biology, Lewis-Sigler Institute for Integrative Genomics, Princeton University, Howard Hughes Medical Institute, Chevy Chase, MD, USA
* **Scott Wolf**, Research Computing, Princeton University, Princeton, NJ, USA
* **Ian M Traniello**, Department of Ecology and Evolutionary Biology and Lewis-Sigler Institute for Integrative Genomics, Princeton University
* **Sarah Kocher**, Department of Ecology and Evolutionary Biology, Lewis-Sigler Institute for Integrative Genomics, Princeton University, Howard Hughes Medical Institute, Chevy Chase, MD, USA

=======
License
=======

Pipemake is licensed under the MIT license. See the `LICENSE <https://github.com/kocherlab/pipemake/blob/main/LICENSE>`_ file for details.
