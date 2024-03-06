rule all:
	input:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['species']}.vcf.gz")

rule reseq_phased_map_plink:
	input:
		os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.vcf.gz')
	output:
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.map')),
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.ped')),
		temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.log')),
	params:
		out_prefix=os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}')
	singularity:
		"/home/aewebb/plink.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"plink2 --vcf {input} --export ped --out {params.out_prefix} --set-all-var-ids @:# --allow-extra-chr"

rule reseq_ihs_selscan:
	input:
		vcf=os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.vcf.gz'),
		map=os.path.join(config['paths']['reseq_vcf_phased_dir'], 'SplitByChrom', '{chrom}.map')
	output:
		temp(os.path.join(config['paths']['reseq_popgen_anlyses'], 'ihs', '{chrom}.ihs.out')),
		temp(os.path.join(config['paths']['reseq_popgen_anlyses'], 'ihs', '{chrom}.ihs.log'))
		
	params:
		out_prefix=os.path.join(config['paths']['reseq_popgen_anlyses'], 'ihs', '{chrom}'),
		maf = config['maf']
	singularity:
		"/home/aewebb/selscan.sif"
	resources:
		mem_mb=24000
	threads: 12
	shell:
		"selscan --ihs --vcf {input.vcf} --map {input.map} --pmap --maf {params.maf} --threads {threads} --out {params.out_prefix}"

def aggregate_ihs_reseq (wildcards):
	checkpoint_output = checkpoints.reseq_split_unphased_bcftools.get(**wildcards).output[0]
	return {'ihs': expand(os.path.join(config['paths']['reseq_popgen_anlyses'], 'ihs', '{chrom}.ihs.out'), chrom = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.vcf.gz")).chrom),
			'log': expand(os.path.join(config['paths']['reseq_popgen_anlyses'], 'ihs', '{chrom}.ihs.log'), chrom = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.vcf.gz")).chrom)}

rule reseq_cat_phased_bcftools:
	input:
		unpack(aggregate_ihs_reseq)
	output:
		os.path.join(config['paths']['reseq_popgen_anlyses'], 'ihs', f"{config['species']}.ihs.out")
		os.path.join(config['paths']['reseq_popgen_anlyses'], 'ihs', f"{config['species']}.ihs.log")
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"""
		cat {input.ihs} > {output.ihs}
		cat {input.log} > {output.log}
		"""