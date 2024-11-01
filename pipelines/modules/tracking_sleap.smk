rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["tracking_sleap_dir"],
                "{video_id}_sleap_tracked.slp",
            ),
            video_id=config["video_ids"],
        ),


rule run_sleap:
    input:
        input_video=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["tracking_videos_dir"],
            "{video_id}.mp4",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["tracking_sleap_dir"],
            "{video_id}_sleap_tracked.slp",
        ),
    params:
        centroid_model=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["tracking_model_dir"],
            "centroid_model",
        ),
        instance_model=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["tracking_model_dir"],
            "instance_model",
        ),
        tracking_batch_size=config["tracking_batch_size"],
        tracking_tracker=config["tracking_tracker"],
        tracking_max_tracking=config["tracking_max_tracking"],
    singularity:
        "docker://swwolf/sleap:latest"
    shell:
        """
        sleap-track \
        {input.input_video} \
        -m {params.centroid_model} \
        -m {params.instance_model} \
        --output {output} \
        --verbosity json \
        --batch-size {params.tracking_batch_size} \
        --tracking.tracker {params.tracking_tracker} \
        --tracking.max_tracking {params.tracking_max_tracking} \
        """
