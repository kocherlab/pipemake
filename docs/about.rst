.. _about:

#####
About
#####

*****************
What is pipemake?
*****************

The goal of pipemake is to provide a lightweight, flexible, and easy-to-use tool for managing and creating `Snakemake <https://snakemake.readthedocs.io/>`_ pipelines. Pipelines are written in Snakemake, to easily optimize computational efficiency and reproducibility.

For the majority of users, pipemake provides a collection of curated and customizable genomic analysis pipelines that can be easily integrated into research projects. Running most pipelines simply requires users to specify the desired pipeline and the appropriate input files. 

Pipemake will then automatically create a workflow directory, which stores all the necessary files and configurations for running the pipeline. The workflow directories then can be executed using Snakemake, which will handle the execution of all pipeline steps.

.. image:: _static/Pipemake_workflow_figure_v3.jpg
   :alt: The workflow directory contains all the Snakemake files, configuration files, input files, and output files.

Workflow directories were devised as a method to simplify record keeping and reproducibility. They enable users to always have a complete record of the pipeline run, including the exact Snakemake rules used, the input files, and the output files generated. This is particularly useful for sharing results with collaborators or for future reference.

******************
Creating Pipelines
******************

For users who want to create their own pipelines, pipemake provides a framework that allows for simplified development of new pipelines.

pipemake uses configurable YAML files to define pipeline, including:

* The command-line arguments for the pipeline
* The process to correctly store the input files (i.e. naming conventions, if the input files are compressed, etc.)
* The required snakemake files
* The snakemake rules to link if their input/output files do not match

While creating configurable YAML may take more time then simply combining existing Snakemake rules, it allows for pipelines that be easily reused and modified.

For a detailed guide on how to create pipelines (see :ref:`create`).
