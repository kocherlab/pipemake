rule all:
	input:
        expand(os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec"), model=config['models']),
        expand(os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenval"), model=config['models']),

rule reseq_model_calc_pca_plink:
	input:
		bed_file=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bed"),
		bim_file=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bim"),
		fam_file=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.fam"),
        ind_file=os.path.join(config['paths']['models_dir'], f"{config['species']}.{{model}}.ind.txt")
	output:
		os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec"),
        os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenval"),
        os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.log")
	params:
		bed_prefix=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned"),
        pca_prefix=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA'. '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca"),
        pca_count=config['pca_count']
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/pipemake_utils.sif"
	shell:
		"plink2 --bfile {params.bed_prefix} --keep {input.ind_file} --pca {params.pca_count} --allow-extra-chr --out {params.pca_prefix}"

rule reseq_model_plot_pca:
	input:
		eigenvec_file=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec"),
        eigenval_file=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenval"),
        model_file=os.path.join(config['paths']['models_dir'], f"{config['species']}.model")
	output:
		os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec"),
        os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.eigenval"),
        os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca.log")
	params:
        pca_dir=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}'),
        out_prefix=os.path.join(config['paths']['reseq_popgen_dir'], 'PCA', '{model}', f"{config['species']}_{config['assembly_version']}.pruned.pca")
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/pipemake_utils.sif"
	shell:
		"plot-pca --pca-dir {params.bed_prefix} --model-file {input.model_file} --model-name {wildcards.modelt} --out-prefix {params.out_prefix}"