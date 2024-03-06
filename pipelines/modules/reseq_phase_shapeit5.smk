rule all:
	input:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.vcf.gz")

checkpoint split_by_chrom_bcftools:
	input:
		os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['vcf_prefix']}.vcf.gz")
	output:
		temp(directory(os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'Split')))
	params:
		out_dir=os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'Split'),
		out_prefix=os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'Split', '')
	singularity:
		"/home/aewebb/kocherPOP.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"""
		mkdir {params.out_dir}
		bcftools index -f {input}
		bcftools index -s {input} | cut -f 1 | while read chrom; do bcftools view -O b -o {params.out_prefix}${{chrom}}.bcf {input} '${{chrom}}'; done
		"""
		
rule phase_chroms_shapeit5:
	input:
		os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'Split', '{chrom}.bcf')
	output:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], 'Split', '{chrom}.bcf')
	singularity:
		"/home/aewebb/kocherPOP.sif"
	resources:
		mem_mb=24000
	threads: 12
	shell:
		"""
		bcftools index -f {input}
		SHAPEIT5_phase_common --input {input} --region {wildcards.chrom} --output {output} --thread {threads}
		"""

def aggregate_phased (wildcards):
	checkpoint_output = checkpoints.split_by_chrom_bcftools.get(**wildcards).output[0]
	return expand(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'Split', '{chrom}.bcf'), chrom = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.bcf")).chrom)

rule cat_phased_chroms_bcftools:
	input:
		aggregate_phased
	output:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.vcf.gz")
	params:
		tmp_dir=os.path.join(config['paths']['reseq_vcf_unphased_dir'], 'Split')
	singularity:
		"/home/aewebb/kocherPOP.sif"
	resources:
		mem_mb=8000
	threads: 4
	shell:
		"""
		bcftools concat --threads {threads} -O z -o {output} {input}
		rm -r {params.tmp_dir}
		"""