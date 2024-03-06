rule all:
	input:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.vcf.gz")

checkpoint split_by_chrom_bcftools:
	input:
		os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.vcf.gz")
	output:
		temp(directory(os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'SplitByChrom')))
	params:
		out_dir=os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'SplitByChrom'),
		out_prefix=os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'SplitByChrom', '')
	singularity:
		"/home/aewebb/kocherPOP.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"""
		mkdir {params.out_dir}
		bcftools index -f {input}
		bcftools index -s {input} | cut -f 1 | while read chrom; do bcftools view -O z -o {params.out_prefix}${{chrom}}.vcf.gz {input} '${{chrom}}'; done
		"""

rule phase_chroms_shapeit2:
	input:
		os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'SplitByChrom', '{chrom}.vcf.gz')
	output:
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.haps')),
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.sample')),
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.phase.snp.mm')),
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.phase.ind.mm')),
		log=os.path.join(config['paths']['reseq_vcf_phased_dir'], '{chrom}.phase.log')
	params:
		haps_prefix=os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}')
	singularity:
		"/home/aewebb/kocherPOP.sif"
	resources:
		mem_mb=8000
	threads: 1
	shell:
		"shapeit -V {input} -O {params.haps_prefix} --output-log {output.log}"

rule convert_chroms_shapeit2:
	input:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.haps'),
		os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.sample')
	output:
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.shapeit_header.vcf.gz')),
		log=os.path.join(config['paths']['reseq_vcf_phased_dir'], '{chrom}.convert.log')
	params:
		haps_prefix=os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}'),
		shapeit_vcf=os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.shapeit_header.vcf')
	singularity:
		"/home/aewebb/kocherPOP.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"""
		shapeit -convert --input-haps {params.haps_prefix} --output-vcf {params.shapeit_vcf} --output-log {output.log}
		bgzip {params.shapeit_vcf}
		"""

def aggregate_phased (wildcards):
	checkpoint_output = checkpoints.split_by_chrom_bcftools.get(**wildcards).output[0]
	return expand(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.shapeit_header.vcf.gz'), chrom = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.shapeit_header.vcf.gz")).chrom)

rule cat_phased_chroms_bcftools:
	input:
		aggregate_phased
	output:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.shapeit_header.vcf.gz")
	singularity:
		"/home/aewebb/kocherPOP.sif"
	resources:
		mem_mb=8000
	threads: 4
	shell:
		"""
		bcftools concat --threads {threads} -O z -o {output} {input}
		"""

rule bcftools_create_header:
	input:
		os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.vcf.gz")
	output:
		temp(os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.header"))
	params:
		haps_output=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.haps"),
		shapeit_vcf=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.shapeit_header.vcf.gz")
	shell:
		"""
		bcftools view -h {input} > {output}
		date=$(date +'%a %b %H:%M:%S %Y')
		insert_pos=$(wc -l < {output})
		awk -v n=$insert_pos -v s="##shapeit_covertCommand=-convert --input-haps {params.haps_output} --output-vcf {params.shapeit_vcf}; Date=$date" 'NR == n {{print s}} {{print}}' {output} > {output}.tmp && mv {output}.tmp {output}
		awk -v n=$insert_pos -v s="##shapeit_phaseCommand=-V {input} -O {params.haps_output};  Date=$date" 'NR == n {{print s}} {{print}}' {output} > {output}.tmp && mv {output}.tmp {output}
		awk -v n=$insert_pos -v s="##bcftools_pipelineCommand=view {params.shapeit_vcf} | reheader {output};  Date=$date" 'NR == n {{print s}} {{print}}' {output} > {output}.tmp && mv {output}.tmp {output}
		"""

rule bcftools_replace_shapit_header:
	input:
		shapeit_vcf=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.shapeit_header.vcf.gz"),
		header=os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['species']}.header")
	output:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.vcf.gz")
	shell:
		"bcftools view {input.shapeit_vcf} | bcftools reheader -h {input.header} | bcftools view -O z -o {output}"