rule all:
	input:
		os.path.join(config['paths']['reseq_filtered_vcf_dir'], f"{config['species']}_{config['assembly_version']}.filtered.vcf.gz")

rule convert_flitered_vcf_to_plink:
	input:
		os.path.join(config['paths']['reseq_filtered_vcf_dir'], f"{config['species']}_{config['assembly_version']}.filtered.vcf.gz")
	output:
		os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bed"),
		os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.bim"),
		os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered.fam")
	params:
		out_prefix=os.path.join(config['paths']['reseq_filtered_plink_dir'], f"{config['species']}_{config['assembly_version']}.filtered")
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/plink.sif"
	shell:
		"plink2 --vcf {input} --make-bed --out {params.out_prefix} --allow-extra-chr --set-missing-var-ids @:# --threads {threads}"