pipeline_dir=%CONDA_PREFIX%/share/pipelines
mkdir %pipeline_dir%
xcopy pipelines %pipeline_dir% /E

mkdir %CONDA_PREFIX%\etc\conda\activate.d
echo set KPDIR=%pipeline_dir% > %PREFIX%\etc\conda\activate.d\%PKG_NAME%_activate.sh

mkdir %CONDA_PREFIX%\etc\conda\deactivate.d
echo set KPDIR= > %PREFIX%\etc\conda\deactivate.d\%PKG_NAME%_deactivate.sh

python -m pip install --no-deps --ignore-installed .
