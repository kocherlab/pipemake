rule all:
	input:
		expand(os.path.join(config['paths']['reseq_popgen_dir'], 'Fst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.fst.summary"), model=config['models']),
		expand(os.path.join(config['paths']['reseq_popgen_dir'], 'Fst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.fst.log"), model=config['models'])

rule reseq_model_fst_phenotype_file:
	input:
		fam_file=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.fam")
		model_file=os.path.join(config['paths']['models_dir'], f"{config['species']}.model")
	output:
		os.path.join(config['paths']['models_dir'], 'Fst', f"{config['species']}.{{model}}.pheno.txt"),
		os.path.join(config['paths']['models_dir'], 'Fst', f"{config['species']}.{{model}}.pheno.log")
	params:
		out_prefix=os.path.join(config['paths']['models_dir'], 'Fst', f"{config['species']}.{{model}}")
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/pipemake_utils.sif"
	shell:
		"ped-phenotype-file --fam {input.fam_file} --model-file {input.model_file} --model-name {wildcards.model} --out-prefix {params.out_prefix} --pheno-header {wildcards.model}"

rule reseq_model_calc_fst_plink:
	input:
		bed_file=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bed"),
		bim_file=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bim"),
		fam_file=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.fam"),
		pheno_file=os.path.join(config['paths']['models_dir'], 'Fst', f"{config['species']}.{{model}}.pheno.txt"),
		ind_file=os.path.join(config['paths']['models_dir'], f"{config['species']}.{{model}}.ind.txt")
	output:
		os.path.join(config['paths']['reseq_popgen_dir'], 'Fst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.fst.summary"),
		os.path.join(config['paths']['reseq_popgen_dir'], 'Fst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.fst.log")
	params:
		bed_prefix=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered"),
		fst_prefix=os.path.join(config['paths']['reseq_popgen_dir'], 'Fst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered"),
		fst_method=config['fst_method']
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/pipemake_utils.sif"
	shell:
		"plink2 --bfile {params.bed_prefix} --pheno {input.pheno_file} --keep {input.ind_file} --fst {wildcards.model} report-variants method={params.fst_method} --allow-extra-chr --out {params.fst_prefix}"