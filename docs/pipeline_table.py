import os
import yaml

# Assign the pipelines table filename
pipeline_table_filename = os.path.join("assets", "pipelines.csv")

# Open the pipelines table file
with open(pipeline_table_filename, "w") as table_file:
    # Write the header to the pipelines table file
    table_file.write("Pipeline,Description\n")

    # Assign the pipeline directory
    pipeline_dir = os.path.join("..", "pipelines", "configs")

    # Loop the pipeline files in the pipeline directory
    for pipeline_filename in os.listdir(pipeline_dir):
        if not pipeline_filename.endswith(".yml") and not pipeline_filename.endswith(
            ".yaml"
        ):
            continue

        # Assign the pipeline file path
        pipeline_path = os.path.join(pipeline_dir, pipeline_filename)

        # Load the pipeline file
        with open(pipeline_path, "r") as file:
            pipeline_data = yaml.safe_load(file)

        # Write the pipeline data to the pipelines table file
        table_file.write(
            f'{pipeline_data["pipeline"]},{pipeline_data["parser"]["help"]}\n'
        )
