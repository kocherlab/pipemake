rule all:
	input:
		os.path.join(config['paths']['models_dir'], f"{config['species']}.{config['model_name']}.ind.txt"),
		os.path.join(config['paths']['models_dir'], f"{config['species']}.{config['model_name']}.ind.log")

rule model_name_ind_file:
	input:
		os.path.join(config['paths']['models_dir'], f"{config['species']}.model")
	output:
		os.path.join(config['paths']['models_dir'], f"{config['species']}.{config['model_name']}.ind.txt"),
		os.path.join(config['paths']['models_dir'], f"{config['species']}.{config['model_name']}.ind.log")
	params:
		out_prefix=os.path.join(config['paths']['models_dir'], f"{config['species']}.{config['model_name']}"),
		model_name=config['model_name']
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/pipemake_utils.sif"
	shell:
		"model-inds --model-file {input} --model-name {params.model_name} --out-prefix {params.out_prefix}"