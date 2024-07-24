rule all:
	input:
		expand(os.path.join(config['paths']['reseq_popgen_dir'], f"{config['species']}.{{model}}.report"), model=config['models'])

rule reseq_model_fst_phenotype_file:
	input:
		fam_file=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.fam"),
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
		"docker://aewebb/pipemake_utils:v0.1.27"
	shell:
		"ped-phenotype-file --fam {input.fam_file} --model-file {input.model_file} --model-name {wildcards.model} --out-prefix {params.out_prefix} --out-format plink2 --pheno-header {wildcards.model}"

checkpoint reseq_model_calc_fst_plink:
	input:
		bed_file=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bed"),
		bim_file=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bim"),
		fam_file=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.fam"),
		pheno_file=os.path.join(config['paths']['models_dir'], 'Fst', f"{config['species']}.{{model}}.pheno.txt"),
		ind_file=os.path.join(config['paths']['models_dir'], f"{config['species']}.{{model}}.ind.txt")
	output:
		fst_dir=directory(os.path.join(config['paths']['reseq_popgen_dir'], 'Fst', '{model}'))
	params:
		bed_prefix=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered"),
		fst_prefix=os.path.join(config['paths']['reseq_popgen_dir'], 'Fst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered"),
		fst_method=config['fst_method']
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"docker://aewebb/plink2:20240418"
	shell:
		"""
		mkdir -p {output.fst_dir}
		plink2 --bfile {params.bed_prefix} --pheno {input.pheno_file} --keep {input.ind_file} --fst {wildcards.model} report-variants method={params.fst_method} --allow-extra-chr --out {params.fst_prefix}
		"""

def get_fst_files (wildcards):
	checkpoint_output = checkpoints.reseq_model_calc_fst_plink.get(**wildcards).output['fst_dir']
	return expand(os.path.join(checkpoint_output, f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.var"),
				  pair = glob_wildcards(os.path.join(checkpoint_output, f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.var")).pair)

rule reseq_model_fst_tmp_report:
	input:
		get_fst_files
	output:
		temp(os.path.join(config['paths']['reseq_popgen_dir'], f"{config['species']}.{{model}}.report"))
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"echo {input} > {output}"