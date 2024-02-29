rule all:
	input:
		expand(os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.lmm.filtered.assoc.txt"), category=config['categories'])

rule category_ind_file:
	input:
		os.path.join(config['paths']['models'], f"{config['species']}.model")
	output:
		temp(os.path.join(config['paths']['models'], 'GEMMA', f"{config['species']}.{{category}}.ind.txt")),
		temp(os.path.join(config['paths']['models'], 'GEMMA', f"{config['species']}.{{category}}.ind.log"))
	params:
		out_prefix=os.path.join(config['paths']['models'], 'GEMMA', f"{config['species']}.{{category}}")
	shell:
		"model-category-inds --model {input} --model-category {wildcards.category} --out-prefix {params.out_prefix}"

rule gemma_category_bed:
	input:
		bed_file=os.path.join(config['paths']['reseq_ped_dir'], f"{config['species']}_{config['assembly_version']}.bed"),
		bim_file=os.path.join(config['paths']['reseq_ped_dir'], f"{config['species']}_{config['assembly_version']}.bim"),
		fam_file=os.path.join(config['paths']['reseq_ped_dir'], f"{config['species']}_{config['assembly_version']}.fam"),
		ind_file=os.path.join(config['paths']['models'], 'GEMMA', f"{config['species']}.{{category}}.ind.txt")
	output:
		temp(os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.bed")),
		temp(os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.bim")),
		temp(os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.fam"))
	params:
		bed_prefix=os.path.join(config['paths']['reseq_gwas_dir'], f"{config['species']}_{config['assembly_version']}"),
		out_prefix=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}")
	shell:
		"plink --bfile {input.bed_prefix} --keep {input.ind_file} --allow-extra-chr --make-bed --out {params.out_prefix}"

rule category_phenotype_file:
	input:
		fam_file=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.fam"),
		model_file=os.path.join(config['paths']['models'], f"{config['species']}.model")
	output:
		temp(os.path.join(config['paths']['models'], 'GEMMA', f"{config['species']}.{{category}}.pheno.txt")),
		temp(os.path.join(config['paths']['models'], 'GEMMA', f"{config['species']}.{{category}}.pheno.log"))
	params:
		out_prefix=os.path.join(config['paths']['models'], f"{config['species']}.{{category}}")
	shell:
		"ped-phenotype-file --fam {input.fam_file} --model {input.model_file} --model-category {wildcards.category} --out-prefix {params.out_prefix}"

rule run_gemma_gk:
	input:
		bed_file=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.bed"),
		bim_file=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.bim"),
		fam_file=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.fam"),
		pheno_file=os.path.join(config['paths']['models'], 'GEMMA', f"{config['species']}.{{category}}.pheno.txt")
	output:
		os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.gk.cXX.txt"),
		os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.gk.log.txt")
	params:
		bed_prefix=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}"),
		out_prefix=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.gk"),
		kinship_matrix = config['kinship_matrix']
	shell:
		"gemma -bfile {params.bed_prefix} -p {input.pheno_file} -gk {params.kinship_matrix} -o {params.out_prefix}"

rule run_gemma_lmm:
	input:
		bed_file=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.bed"),
		bim_file=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.bim"),
		fam_file=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.fam"),
		pheno_file=os.path.join(config['paths']['models'], 'GEMMA', f"{config['species']}.{{category}}.pheno.txt"),
		gk_file=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.gk.cXX.txt")
	output:
		os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.lmm.assoc.txt"),
		os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.lmm.log.txt")
	params:
		bed_prefix=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}"),
		out_prefix=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.lmm"),
		linear_model = config['linear_model'],
		maf = config['maf']
	shell:
		"gemma -bfile {params.bed_prefix} -p {input.pheno_file} -lmm {params.linear_model} -k {input.gk_file} -maf {params.maf} -o {params.out_prefix}"

rule filter_gemma:
	input:
		os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.lmm.assoc.txt")
	output:
		os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.lmm.filtered.assoc.txt"),
		os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.lmm.filtered.log")
	params:
		out_prefix=os.path.join(config['paths']['reseq_gwas_dir'], 'GEMMA', f"{config['species']}_{config['assembly_version']}.{{category}}.lmm"),
		min_log_pvalue = config['min_log_pvalue']
	shell:
		"filter-gemma --gemma-file {input} --min-log-pvalue {params.min_log_pvalue} --out-prefix {params.out_prefix}"