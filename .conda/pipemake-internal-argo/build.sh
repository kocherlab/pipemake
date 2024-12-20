#!/usr/bin/env bash

pipeline_dir="/Genomics/kocherlab/lab/Pipelines/pipemake/pipelines"
image_dir="/Genomics/kocherlab/lab/Pipelines/images"

mkdir -p "${PREFIX}/etc/conda/activate.d"
echo "export PM_SNAKEMAKE_DIR=${pipeline_dir}" > "${PREFIX}/etc/conda/activate.d/${PKG_NAME}_activate.sh"
echo "export PM_SINGULARITY_DIR=${image_dir}" >> "${PREFIX}/etc/conda/activate.d/${PKG_NAME}_activate.sh"

mkdir -p "${PREFIX}/etc/conda/deactivate.d"
echo "unset PM_SNAKEMAKE_DIR" > "${PREFIX}/etc/conda/deactivate.d/${PKG_NAME}_deactivate.sh"
echo "unset PM_SINGULARITY_DIR" >> "${PREFIX}/etc/conda/deactivate.d/${PKG_NAME}_deactivate.sh"

python -m pip install --no-deps --ignore-installed .