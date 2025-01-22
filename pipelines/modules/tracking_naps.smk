rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["tracking_naps_dir"],
                "{video_id}_naps_tracked.slp",
            ),
            video_id=config["video_ids"],
        ),


rule run_naps:
    input:
        sleap_output=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["tracking_sleap_dir"],
            "{video_id}_sleap_tracked.slp",
        ),
        input_video=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["tracking_videos_dir"],
            "{video_id}.mp4",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["tracking_naps_dir"],
            "{video_id}_naps_tracked.slp",
        ),
    params:
        output_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["tracking_naps_dir"],
            "{video_id}_naps_tracked",
        ),
        start_frame=config["start_frame"],
        end_frame=config["end_frame"],
        tag_node_name=config["tag_node_name"],
        aruco_marker_set=config["aruco_marker_set"],
        aruco_error_correction_rate=config["aruco_error_correction_rate"],
        aruco_adaptive_thresh_constant=config["aruco_adaptive_thresh_constant"],
        aruco_adaptive_thresh_win_size_max=config["aruco_adaptive_thresh_win_size_max"],
        aruco_adaptive_thresh_win_size_step=config[
            "aruco_adaptive_thresh_win_size_step"
        ],
        aruco_adaptive_thresh_win_size_min=config["aruco_adaptive_thresh_win_size_min"],
        half_rolling_window_size=config["half_rolling_window_size"],
        aruco_crop_size=config["aruco_crop_size"],
    singularity:
        "docker://swwolf/naps:latest"
    resources:
        mem_mb=32000,
    threads: 1
    shell:
        """
        naps-track \
        --slp-path {input.sleap_output} \
        --video-path {input.input_video} \
        --tag-node-name {params.tag_node_name} \
        --start-frame {params.start_frame} \
        --end-frame {params.end_frame} \
        --aruco-marker-set {params.aruco_marker_set} \
        --aruco-error-correction-rate {params.aruco_error_correction_rate} \
        --aruco-adaptive-thresh-constant {params.aruco_adaptive_thresh_constant} \
        --aruco-adaptive-thresh-win-size-max {params.aruco_adaptive_thresh_win_size_max} \
        --aruco-adaptive-thresh-win-size-step {params.aruco_adaptive_thresh_win_size_step} \
        --aruco-adaptive-thresh-win-size-min {params.aruco_adaptive_thresh_win_size_min} \
        --half-rolling-window-size {params.half_rolling_window_size} \
        --aruco-crop-size {params.aruco_crop_size} \
        --output {params.output_prefix}
        """

