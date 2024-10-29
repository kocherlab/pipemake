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
        tracking_tracker=config["tracking_tracker"],
        tracking_similarity=config["tracking_similarity"],
        tracking_max_tracking=config["tracking_max_tracking"],
        tracking_min_new_track_points=config["tracking_min_new_track_points"],
        tracking_min_match_points=config["tracking_min_match_points"],
        tracking_post_connect_single_breaks=config[
            "tracking_post_connect_single_breaks"
        ],
        tracking_target_instance_count=config["tracking_target_instance_count"],
    singularity:
        "docker://swwolf/sleap:latest"
    shell:
        """
        sleap-track \
        {input.input_video} \
        -m {input.centroid_model} \
        -m {input.instance_model} \
        --output "{output}" \
        --verbosity json \
        --tracking.tracker {params.tracking_tracker} \
        --tracking.similarity {params.tracking_similarity} \
        --tracking.max_tracking {params.tracking_max_tracking} \
        --tracking.min_new_track_points {params.tracking_min_new_track_points} \
        --tracking.min_match_points {params.tracking_min_match_points} \
        --tracking.post_connect_single_breaks {params.tracking_post_connect_single_breaks} \
        --tracking.target_instance_count {params.tracking_target_instance_count} \
        """
