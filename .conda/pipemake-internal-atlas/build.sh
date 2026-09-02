#!/usr/bin/env bash

image_dir="/project/beenome100_collab/conda_envs"

mkdir -p "${PREFIX}/etc/conda/activate.d"
echo "export PM_SINGULARITY_DIR=${image_dir}" >> "${PREFIX}/etc/conda/activate.d/${PKG_NAME}_activate.sh"

mkdir -p "${PREFIX}/etc/conda/deactivate.d"
echo "unset PM_SINGULARITY_DIR" >> "${PREFIX}/etc/conda/deactivate.d/${PKG_NAME}_deactivate.sh"

python -m pip install --no-deps --ignore-installed .