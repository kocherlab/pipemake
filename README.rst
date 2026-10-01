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
pipemake is a lightweight, flexible, and easy-to-use tool for creating and managing `Snakemake <https://snakemake.readthedocs.io/>`_ pipelines. It makes Snakemake's reproducibility safeguards automatic, reusable, and accessible to researchers who would otherwise struggle to implement them. It was designed with four primary goals:

1. Offer a collection of curated, customizable genomic analysis pipelines for researchers seeking to rapidly integrate Snakemake-based workflows into their research.
2. Build safeguards into every pipeline as enforced defaults, including validated command-line arguments, mandatory containerization, and standardized record-keeping in a workflow directory.
3. Optimize computational efficiency and reproducibility by fully operating in the Snakemake ecosystem. pipemake generates standard Snakemake files, so a workflow can be run with Snakemake alone.
4. Streamline development and sharing by creating a flexible platform with swappable pipelines that easily reuse previously written Snakemake code, allowing a group to maintain and share a single collection of pipelines.

pipemake is most useful for analyses that need to be repeated, shared, or reproduced by someone other than the original designer, such as in research groups whose members have different levels of computational expertise. It is not intended to replace community-curated pipelines, such as those from `nf-core <https://nf-co.re/>`_, where one already exists for the analysis.

================
Getting pipemake
================

-----
mamba
-----

.. code-block:: bash

   mamba create -c conda-forge -c bioconda kocherlab::pipemake

For more information, see the `installation instructions <https://pipemake.readthedocs.io/en/latest/installation.html>`_.

======
Issues
======

1. Check the `docs <https://pipemake.rtfd.io/>`_.
2. Search the `issues on GitHub <https://github.com/kocherlab/pipemake/issues>`_ or open a new one.

============
Contributors
============

* **Andrew Webb**, Department of Integrative Biology, University of California Berkeley, Berkeley, CA, USA, Howard Hughes Medical Institute, Chevy Chase, MD, USA
* **Scott Wolf**, Department of Integrative Biology, University of California Berkeley, Berkeley, CA, USA
* **Ian M Traniello**, Department of Ecology and Evolutionary Biology and Lewis-Sigler Institute for Integrative Genomics, Princeton University, Princeton, NJ, USA
* **Sarah Kocher**, Department of Integrative Biology, University of California Berkeley, Berkeley, CA, USA, Howard Hughes Medical Institute, Chevy Chase, MD, USA

=======
License
=======

Pipemake is licensed under the MIT license. See the `LICENSE <https://github.com/kocherlab/pipemake/blob/main/LICENSE>`_ file for details.
