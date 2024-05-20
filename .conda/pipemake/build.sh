#!/usr/bin/env bash

pipeline_dir="${PREFIX}/share/pipelines"
mkdir -p $pipeline_dir
cp -r pipelines/* $pipeline_dir

mkdir -p "${PREFIX}/etc/conda/activate.d"
echo "export KPDIR=${pipeline_dir}" > "${PREFIX}/etc/conda/activate.d/${PKG_NAME}_activate.sh"

mkdir -p "${PREFIX}/etc/conda/deactivate.d"
echo "unset KPDIR=${pipeline_dir}" > "${PREFIX}/etc/conda/deactivate.d/${PKG_NAME}_deactivate.sh"

python -m pip install --no-deps --ignore-installed .