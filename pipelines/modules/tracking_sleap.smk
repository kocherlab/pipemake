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
        centroid_model=directory(
            config["paths"]["workflow_prefix"],
            config["paths"]["model_dir"],
            "SLEAP",
            "centroid_model",
        ),
        instance_model=directory(
            config["paths"]["workflow_prefix"],
            config["paths"]["model_dir"],
            "SLEAP",
            "instance_model",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["tracking_sleap_dir"],
            "{video_id}_sleap_tracked.slp",
        ),
    params:
        tracking_batch_size=config["tracking_batch_size"],
        tracking_tracker=config["tracking_tracker"],
        tracking_max_tracking=config["tracking_max_tracking"],
        tracking_max_tracks=config["tracking_max_tracks"],
    singularity:
        "docker://swwolf/sleap:latest"
    shell:
        """
        sleap-track \
        {input.input_video} \
        -m {input.centroid_model} \
        -m {input.instance_model} \
        --output {output} \
        --verbosity json \
        --batch-size {params.tracking_batch_size} \
        --tracking.tracker {params.tracking_tracker} \
        --tracking.max_tracking {params.tracking_max_tracking} \
        --tracking.max_tracks {params.tracking_max_tracks} \
        """
