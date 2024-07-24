rule all:
	input:
		expand(os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec"), model=config['models']),
		expand(os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenval"), model=config['models']),
		expand(os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.pdf"), model=config['models'])

rule reseq_calc_freq_plink:
	input:
		bed_file=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bed"),
		bim_file=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bim"),
		fam_file=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.fam"),
		ind_file=os.path.join(config['paths']['models_dir'], f"{config['species']}.{{model}}.ind.txt")
	output:
		os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.afreq")
	params:
		bed_prefix=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned"),
		pca_prefix=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned")
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"docker://aewebb/plink2:20240418"
	shell:
		"plink2 --bfile {params.bed_prefix} --keep {input.ind_file} --freq --allow-extra-chr --out {params.pca_prefix}"

rule reseq_model_calc_pca_plink:
	input:
		bed_file=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bed"),
		bim_file=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bim"),
		fam_file=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.fam"),
		ind_file=os.path.join(config['paths']['models_dir'], f"{config['species']}.{{model}}.ind.txt"),
		freq_file=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.afreq")
	output:
		os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec"),
		os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenval"),
		os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.log")
	params:
		bed_prefix=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned"),
		pca_prefix=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca"),
		pca_count=config['pca_count']
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"docker://aewebb/plink2:20240418"
	shell:
		"plink2 --bfile {params.bed_prefix} --keep {input.ind_file} --read-freq {input.freq_file} --pca {params.pca_count} --allow-extra-chr --out {params.pca_prefix}"

rule reseq_model_plot_pca:
	input:
		eigenvec_file=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec"),
		eigenval_file=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenval"),
		model_file=os.path.join(config['paths']['models_dir'], f"{config['species']}.model")
	output:
		os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.pdf")
	params:
		pca_dir=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}'),
		out_prefix=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca")
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"docker://aewebb/pipemake_utils:v0.1.27"
	shell:
		"plot-pca --pca-dir {params.pca_dir} --model-file {input.model_file} --model-name {wildcards.model} --out-prefix {params.out_prefix}"