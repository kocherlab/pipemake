rule all:
	input:
		expand(os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.lmm.filtered.assoc.txt"), model=config['models'])

rule gemma_model_bed:
	input:
		bed_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bed"),
		bim_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bim"),
		fam_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.fam"),
		ind_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['models_dir'], f"{config['species']}.{{model}}.ind.txt")
	output:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.bed"),
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.bim"),
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.fam")
	params:
		bed_prefix=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned"),
		out_prefix=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned")
	resources:
		mem_mb=8000
	threads: 1
	singularity:
		"docker://aewebb/plink2:20240418"
	shell:
		"plink2 --bfile {params.bed_prefix} --keep {input.ind_file} --allow-extra-chr --make-bed --out {params.out_prefix}"

rule gemma_model_phenotype_file:
	input:
		fam_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.fam"),
		model_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['models_dir'], f"{config['species']}.model")
	output:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['models_dir'], 'GEMMA', f"{config['species']}.{{model}}.pheno.txt"),
		os.path.join(config['paths']['workflow_prefix'], config['paths']['models_dir'], 'GEMMA', f"{config['species']}.{{model}}.pheno.log")
	params:
		out_prefix=os.path.join(config['paths']['workflow_prefix'], config['paths']['models_dir'], 'GEMMA', f"{config['species']}.{{model}}")
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"docker://aewebb/pipemake_utils:v0.1.27"
	shell:
		"ped-phenotype-file --fam {input.fam_file} --model-file {input.model_file} --model-name {wildcards.model} --numeric --out-prefix {params.out_prefix}"

rule run_gemma_gk:
	input:
		bed_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.bed"),
		bim_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.bim"),
		fam_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.fam"),
		pheno_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['models_dir'], 'GEMMA', f"{config['species']}.{{model}}.pheno.txt")
	output:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.gk.cXX.txt"),
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.gk.log.txt")
	params:
		bed_prefix=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned"),
		out_prefix=f"{config['species']}_{config['assembly_version']}.pruned.gk",
		out_dir=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}'),
		kinship_matrix = config['kinship_matrix']
	resources:
		mem_mb=8000
	threads: 1
	singularity:
		"docker://aewebb/gemma:v0.98.5"
	shell:
		"gemma -bfile {params.bed_prefix} -p {input.pheno_file} -gk {params.kinship_matrix} -outdir {params.out_dir} -o {params.out_prefix}"

rule run_gemma_lmm:
	input:
		bed_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.bed"),
		bim_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.bim"),
		fam_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.fam"),
		pheno_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['models_dir'], 'GEMMA', f"{config['species']}.{{model}}.pheno.txt"),
		gk_file=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.gk.cXX.txt")
	output:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.lmm.assoc.txt"),
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.lmm.log.txt")
	params:
		bed_prefix=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned"),
		out_prefix=f"{config['species']}_{config['assembly_version']}.pruned.lmm",
		out_dir=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}'),
		linear_model = config['linear_model'],
		maf = config['maf']
	resources:
		mem_mb=8000
	threads: 1
	singularity:
		"docker://aewebb/gemma:v0.98.5"
	shell:
		"gemma -p {input.pheno_file} -bfile {params.bed_prefix} -lmm {params.linear_model} -k {input.gk_file} -maf {params.maf} -outdir {params.out_dir} -o {params.out_prefix}"

rule filter_gemma:
	input:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.lmm.assoc.txt")
	output:
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.lmm.filtered.assoc.txt"),
		os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.lmm.filtered.log")
	params:
		out_prefix=os.path.join(config['paths']['workflow_prefix'], config['paths']['reseq_popgen_dir'], 'GEMMA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.lmm"),
		min_log_pvalue = config['min_log_pvalue']
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"docker://aewebb/pipemake_utils:v0.1.27"
	shell:
		"filter-gemma --gemma-file {input} --min-log-pvalue {params.min_log_pvalue} --out-prefix {params.out_prefix}"