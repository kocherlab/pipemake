rule all:
	input:
		expand(os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'ZFst', f"{config['species']}.{{model}}.report"), model=config['models'])

rule reseq_model_calc_zfst_pipemake:
	input:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'Fst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.var")
	output:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'ZFst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.tsv")
	params:
		out_prefix=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'ZFst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst"),
		fst_method=config['fst_method']
	singularity:
		"docker://aewebb/pipemake_utils:v0.1.27"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"""
		fst_method={params.fst_method}
		fst_col=${{fst_method^^}}
		z-normalize --input-file {input} --out-prefix {params.out_prefix} --normalize-col ${{fst_col}}_FST
		"""

rule plot_zfst_pipemake:
	input:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'ZFst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.tsv")
	output:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'ZFst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.manhattan.png")
	params:
		out_prefix=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'ZFst', '{model}', f"{config['species']}_{config['assembly_version']}.{{pair}}.fst")
	singularity:
		"docker://aewebb/pipemake_utils:v0.1.27"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"manhattan-plot --input-file {input} --chrom-col #CHROM --pos-col POS --stat-col 'Z(HUDSON_FST)' --out-prefix {params.out_prefix}"

def get_zfst_files (wildcards):
	checkpoint_output = checkpoints.reseq_model_calc_fst_plink.get(**wildcards).output['fst_dir']
	return expand(os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'ZFst', '{model}', f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.tsv"),
				  model = checkpoint_output.split(os.sep)[-1],
				  pair = glob_wildcards(os.path.join(config['paths']['workflow_prefix'], checkpoint_output, f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.var")).pair)

rule reseq_model_zfst_tmp_report:
	input:
		get_zfst_files
	output:
		temp(os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'ZFst', f"{config['species']}.{{model}}.report"))
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"echo {input} > {output}"