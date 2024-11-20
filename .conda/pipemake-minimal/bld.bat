set pipeline_dir=%CONDA_PREFIX%\share\pipelines

mkdir %CONDA_PREFIX%\etc\conda\activate.d
echo set PM_SNAKEMAKE_DIR=%pipeline_dir% > %PREFIX%\etc\conda\activate.d\%PKG_NAME%_activate.sh

mkdir %CONDA_PREFIX%\etc\conda\deactivate.d
echo set PM_SNAKEMAKE_DIR= > %PREFIX%\etc\conda\deactivate.d\%PKG_NAME%_deactivate.sh

python -m pip install --no-deps --ignore-installed .
