rule all:
	input:
		os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}_{config['assembly_version']}.filtered.missing_data.vcf.gz")

checkpoint pop_ind_file:
	input:
		os.path.join(config['paths']['models_dir'], f"{config['species']}.model")
	output:
		temp(directory(os.path.join(config['paths']['models_dir'], 'BCFtools')))
	params:
		model_name=config['model_name']
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/pipemake_utils.sif"
	shell:
		"model-pop-files --model {input} --model-category {params.model_name} --out-dir {output}"

rule pop_vcf_bcftools:
	input:
		vcf_file=os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}_{config['assembly_version']}.vcf.gz"),
		pop_file=os.path.join(config['paths']['models_dir'], 'BCFtools', "{pop}.pop")
	output:
		temp(os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}_{config['assembly_version']}.filtered.missing_data.{{pop}}.vcf.gz"))
	params:
		missing_cutoff=config['missing_cutoff']
	resources:
		mem_mb=8000
	threads: 1
	singularity:
		"/home/aewebb/kocher_POP.sif"
	shell:
		"""
		bcftools view --samples-file {input.pop_file} {input.vcf_file} | bcftools view -i 'F_MISSING<{params.missing_cutoff}' --output-type z --output-file {output}
		bcftools index {output}
		"""

rule isec_pop_vcfs_bcftools:
	input:
		expand(os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}_{config['assembly_version']}.filtered.missing_data.{{pop}}.vcf.gz"))
	output:
		temp(os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}_{config['assembly_version']}.filtered.missing_data.sites"))
	params:
        pop_count=lambda wildcards, input: len(input)
	resources:
		mem_mb=16000
	threads: 4
	singularity:
		"/home/aewebb/kocher_POP.sif"
	shell:
		"bcftools isec {input} -n={params.pop_count} | cut -f1,2 > {output}"

rule filter_missing_data_vcf_bcftools:
	input:
		vcf_file=os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}_{config['assembly_version']}.vcf.gz"),
		sites_file=os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}_{config['assembly_version']}.filtered.missing_data.sites")
	output:
		os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}_{config['assembly_version']}.filtered.missing_data.vcf.gz")
	resources:
		mem_mb=8000
	threads: 4
	singularity:
		"/home/aewebb/kocher_POP.sif"
	shell:
		"bcftools view -R {input.sites_file} --output-type z --output-file {output} --threads {threads} {input.vcf_file}"
