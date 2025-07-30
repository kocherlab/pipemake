rule all:
    input:
        expand(
            "Tracking/SLEAP/{video_id}_sleap_tracked.slp", video_id=config["video_ids"]
        ),


rule run_sleap:
    input:
        input_video="Tracking/Videos/{video_id}.mp4",
    output:
        "Tracking/SLEAP/{video_id}_sleap_tracked.slp",
    params:
        centroid_model="Tracking/Models/centroid_model",
        instance_model="Tracking/Models/instance_model",
        tracking_batch_size=config["tracking_batch_size"],
        tracking_tracker=config["tracking_tracker"],
        tracking_similarity=(
            f'--tracking.similarity {config["tracking_similarity"]}'
            if "tracking_similarity" in config
            else ""
        ),
        tracking_post_connect_single_breaks=(
            f'--tracking.post_connect_single_breaks {config["tracking_post_connect_single_breaks"]}'
            if "tracking_post_connect_single_breaks" in config
            else ""
        ),
        tracking_pre_cull_to_target=(
            f'--tracking.pre_cull_to_target {config["tracking_pre_cull_to_target"]}'
            if "tracking_pre_cull_to_target" in config
            else ""
        ),
        tracking_target_instance_count=(
            f'--tracking.target_instance_count {config["tracking_target_instance_count"]}'
            if "tracking_target_instance_count" in config
            else ""
        ),
        tracking_clean_instance_count=(
            f'--tracking.clean_instance_count {config["tracking_clean_instance_count"]}'
            if "tracking_clean_instance_count" in config
            else ""
        ),
        tracking_max_tracking=(
            f"--tracking.max_tracking" if config["tracking_max_tracking"] else ""
        ),
        tracking_max_tracks=(
            f'--tracking.max_tracks {config["tracking_max_tracks"]}'
            if "tracking_max_tracks" in config
            else ""
        ),
    singularity:
        "docker://aewebb/sleap:v1.2.8"
    shell:
        """
        sleap-track \
        {input.input_video} \
        --output {output} \
        -m {params.centroid_model} \
        -m {params.instance_model} \
        --verbosity json \
        --batch-size {params.tracking_batch_size} \
        --tracking.tracker {params.tracking_tracker} \
        {params.tracking_similarity} \
        {params.tracking_post_connect_single_breaks} \
        {params.tracking_pre_cull_to_target} \
        {params.tracking_target_instance_count} \
        {params.tracking_clean_instance_count} \
        {params.tracking_max_tracking} \
        {params.tracking_max_tracks}
        """
