rule all:
	input:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.vcf.gz")

rule vcf_header_bcftools:
	input:
		os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.vcf.gz")
	output:
		temp(os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.header"))
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/bcftools.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"bcftools view -h {input} > {output}"

checkpoint split_by_chrom_bcftools:
	input:
		os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.vcf.gz")
	output:
		temp(directory(os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'SplitByChrom')))
	params:
		out_dir=os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'SplitByChrom'),
		out_prefix=os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'SplitByChrom', '')
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/bcftools.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"""
		mkdir {params.out_dir}
		bcftools index -f {input}
		bcftools index -s {input} | cut -f 1 | while read chrom; do bcftools view --regions $chrom -O z -o {params.out_prefix}${{chrom}}.vcf.gz {input}; done
		"""

rule phase_chroms_shapeit4:
	input:
		vcf=os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'SplitByChrom', '{chrom}.vcf.gz'),
		header=os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.header")
	output:
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.vcf.gz'))
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/shapeit4.sif"
	resources:
		mem_mb=24000
	threads: 12
	shell:
		"""
		date=$(date +'%a %b %H:%M:%S %Y')
		insert_pos=$(wc -l < {input.header})
		awk -v n=$insert_pos -v s="##shapeit4_phaseCommand=--input {input.vcf} --region {wildcards.chrom} --output {output} --thread {threads};  Date=$date" 'NR == n {{print s}} {{print}}' {input.header} > {input.header}.tmp && mv {input.header}.tmp {input.header}
		bcftools index -f {input.vcf}
		shapeit4 --input {input.vcf} --region {wildcards.chrom} --output {output} --thread {threads}
		"""

def aggregate_phased (wildcards):
	checkpoint_output = checkpoints.split_by_chrom_bcftools.get(**wildcards).output[0]
	return expand(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.vcf.gz'), chrom = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.vcf.gz")).chrom)

rule cat_phased_chroms_bcftools:
	input:
		vcfs=aggregate_phased,
		header=os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.header")
	output:
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.shapeit_header.vcf.gz"))
	params:
		unphased_split_dir=os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'SplitByChrom')
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/bcftools.sif"
	resources:
		mem_mb=8000
	threads: 4
	shell:
		"""
		date=$(date +'%a %b %H:%M:%S %Y')
		insert_pos=$(wc -l < {input.header})
		awk -v n=$insert_pos -v s="##bcftools_concatCommand=concat --threads {threads} -O z -o {output} {input.vcfs};  Date=$date" 'NR == n {{print s}} {{print}}' {input.header} > {input.header}.tmp && mv {input.header}.tmp {input.header}
		bcftools concat --threads {threads} -O z -o {output} {input.vcfs}
		rm -r {params.unphased_split_dir}
		"""

rule replace_shapit_header_bcftools:
	input:
		shapeit_vcf=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.shapeit_header.vcf.gz"),
		header=os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.header")
	output:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.vcf.gz")
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/bcftools.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"""
		date=$(date +'%a %b %H:%M:%S %Y')
		insert_pos=$(wc -l < {input.header})
		awk -v n=$insert_pos -v s="##bcftools_viewCommand=view {input.shapeit_vcf} | reheader {input.header};  Date=$date" 'NR == n {{print s}} {{print}}' {input.header} > {input.header}.tmp && mv {input.header}.tmp {input.header}
		bcftools view {input.shapeit_vcf} | bcftools reheader -h {input.header} | bcftools view -O z -o {output}
		"""