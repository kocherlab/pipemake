rule all:
	input:
		os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bed"),
		os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bim"),
		os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.fam")

rule reseq_indep_pairwise_plink:
	input:
		os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bed"),
		os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bim"),
		os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.fam")
	output:
		os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.prune.in")
	params:
		bed_prefix=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered"),
		out_prefix=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}"),
		ld_window_size={config['ld_window_size']},
		ld_window_step={config['ld_window_step']},
		ld_threshold={config['ld_threshold']}
	singularity:
		"/Genomics/kocherlab/lab/Pipelines/imagesplink.sif"
	resources:
		mem_mb=8000
	threads: 1
	shell:
		"plink2 --bfile {params.bed_prefix} --indep-pairwise {params.ld_window_size} {params.ld_window_step} {params.ld_threshold} --bad-ld --out {params.out_prefix} --allow-extra-chr --threads {threads}"
	
rule reseq_ld_prune_plink:
	input:
		bed=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bed"),
		bim=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bim"),
		fam=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.fam"),
		prune_in=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.prune.in")
	output:
		os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bed"),
		os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.bim"),
		os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned.fam")
	params:
		input_prefix=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered"),
		output_prefix=os.path.join(config['paths']['reseq_pruned_plink_dir'], f"{config['species']}_{config['assembly_version']}.pruned")
	singularity:
		"/Genomics/kocherlab/lab/Pipelines/imagesplink.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"plink2 --bfile {params.input_prefix} --extract {input.prune_in} --make-bed --out {params.output_prefix} --allow-extra-chr --threads {threads}"