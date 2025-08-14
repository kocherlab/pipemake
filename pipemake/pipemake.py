#!/usr/bin/env python

import os
import shutil

from pipemake.logger import startLogger, logArgDict
from pipemake.parser import PipelineParser
from pipemake.snakemakeIO import SnakePipelineIO
from pipemake.pipelineIO import ConfigPipelinesIO


def main():
    # Assign the pipeline directory
    if os.environ.get("PM_SNAKEMAKE_DIR"):
        pipeline_storage_dir = os.environ.get("PM_SNAKEMAKE_DIR")
    else:
        pipeline_storage_dir = "pipelines"

    # Confirm the pipeline directory exists
    if not os.path.isdir(pipeline_storage_dir):
        raise Exception(f"Unable to find pipeline directory: {pipeline_storage_dir}")

    # Load the pipeline configs
    pipeline_configs = ConfigPipelinesIO.fromDirectory(pipeline_storage_dir)

    # Parse the aguments from the configs
    pipeline_parser = PipelineParser.create(pipeline_configs)
    pipeline_args = pipeline_parser.returnArgs()

    # Assign the pipeline directory to an environment variable, if found
    if os.environ.get("PM_SINGULARITY_DIR") and not pipeline_args["singularity_dir"]:
        pipeline_args["singularity_dir"] = os.environ.get("PM_SINGULARITY_DIR")

    # Update the pipeline args with the pipeline directory
    pipeline_args["pipeline_storage_dir"] = pipeline_storage_dir

    # Check that a pipeline was assigned
    if not pipeline_args["pipeline"]:
        raise Exception("No pipeline specified")

    # Assign the pipeline job directory
    pipeline_args["pipeline_job_dir"] = "pipemake"

    # Check for a previous workflow directory
    if os.path.exists(pipeline_args["workflow_dir"]) and pipeline_args["overwrite"]:
        shutil.rmtree(pipeline_args["workflow_dir"])
    elif os.path.exists(pipeline_args["workflow_dir"]):
        raise Exception(
            f"Workflow directory already exists: {pipeline_args['workflow_dir']}"
        )

    # Create the pipeline job directory
    if not os.path.exists(
        os.path.join(pipeline_args["workflow_dir"], pipeline_args["pipeline_job_dir"])
    ):
        os.makedirs(
            os.path.join(
                pipeline_args["workflow_dir"], pipeline_args["pipeline_job_dir"]
            )
        )

    # Start logger and log the arguments
    startLogger(
        os.path.join(
            pipeline_args["workflow_dir"],
            pipeline_args["pipeline_job_dir"],
            "pipeline.log",
        )
    )
    logArgDict(pipeline_args, omit=["pipeline_job_dir"])

    # Assign the pipeline config
    pipeline_config = pipeline_configs[pipeline_args["pipeline"]]

    # Process the pipeline setup
    pipeline_config.setupPipeline(pipeline_args)

    # Add the samples to the pipeline args
    pipeline_args.update(pipeline_config.samples)

    # Add the setup pipeline args to the pipeline args
    pipeline_args.update(pipeline_config.setup_pipeline_args)

    # Create the snakemake pipeline
    snakemake_pipeline = SnakePipelineIO.open(**pipeline_args)

    # Add the snakemake modules to the pipeline
    for smkm_filename in pipeline_config.snakefiles:
        snakemake_pipeline.addModule(smkm_filename)

    # Add the snakemake links to the pipeline
    for smkm_links in pipeline_config.snakelinks:
        snakemake_pipeline.addSnakeLink(smkm_links)

    # Build the singularity containers
    snakemake_pipeline.buildSingularityContainers()

    # Create the snakemake config file
    snakemake_pipeline.writeConfig(pipeline_args)

    # Create the snakemake piepline file
    snakemake_pipeline.writePipeline()

    # Close the snakemake pipeline
    snakemake_pipeline.close()

    # Print the singularity help message
    pipeline_config.helpMessage(pipeline_args)


if __name__ == "__main__":
    main()
