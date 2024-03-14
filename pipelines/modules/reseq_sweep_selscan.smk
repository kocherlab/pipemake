rule all:
	input:
		os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{config['species']}.ihs.out"),
		os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{config['species']}.ihs.out.log"),
		os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{config['species']}.ihs.norm"),
		os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{config['species']}.ihs.norm.log")

rule reseq_phased_map_plink:
	input:
		os.path.join(config['paths']['reseq_phased_vcf_dir'], 'SplitByChrom', '{chrom}.vcf.gz')
	output:
		temp(os.path.join(config['paths']['reseq_phased_vcf_dir'], 'SplitByChrom', '{chrom}.map')),
		temp(os.path.join(config['paths']['reseq_phased_vcf_dir'], 'SplitByChrom', '{chrom}.ped')),
		temp(os.path.join(config['paths']['reseq_phased_vcf_dir'], 'SplitByChrom', '{chrom}.log'))
	params:
		out_prefix=os.path.join(config['paths']['reseq_phased_vcf_dir'], 'SplitByChrom', '{chrom}')
	singularity:
		"/home/aewebb/plink.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"plink2 --vcf {input} --export ped --out {params.out_prefix} --set-all-var-ids @:# --allow-extra-chr"

rule reseq_ihs_selscan:
	input:
		vcf=os.path.join(config['paths']['reseq_phased_vcf_dir'], 'SplitByChrom', '{chrom}.vcf.gz'),
		map=os.path.join(config['paths']['reseq_phased_vcf_dir'], 'SplitByChrom', '{chrom}.map')
	output:
		temp(os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', '{chrom}.ihs.out')),
		temp(os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', '{chrom}.ihs.log'))
	params:
		out_prefix=os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', '{chrom}'),
		maf = config['maf']
	singularity:
		"/home/aewebb/selscan.sif"
	resources:
		mem_mb=24000
	threads: 12
	shell:
		"selscan --ihs --ihs-detail --vcf {input.vcf} --map {input.map} --pmap --maf {params.maf} --threads {threads} --out {params.out_prefix}"

rule reseq_ihs_normalize_norm:
	input:
		os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', '{chrom}.ihs.out')
	output:
		temp(os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{{chrom}}.ihs.out.{config['bins']}bins.norm")),
		temp(os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{{chrom}}.ihs.out.{config['bins']}bins.log"))
	params:
		out_prefix=os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', '{chrom}.ihs.out'),
		bins = config['bins']
	singularity:
		"/home/aewebb/selscan.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"norm --ihs --files {input} --bins {params.bins} 2> {params.out_prefix}.{params.bins}bins.log"

def aggregate_ihs_reseq (wildcards):
	checkpoint_output = checkpoints.reseq_split_unphased_bcftools.get(**wildcards).output[0]
	return {'scan_ihs': expand(os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', '{chrom}.ihs.out'), chrom = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.vcf.gz")).chrom),
			'scan_log': expand(os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', '{chrom}.ihs.log'), chrom = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.vcf.gz")).chrom),
			'norm_ihs': expand(os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{{chrom}}.ihs.out.{config['bins']}bins.norm"), chrom = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.vcf.gz")).chrom),
			'norm_log': expand(os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{{chrom}}.ihs.out.{config['bins']}bins.log"), chrom = glob_wildcards(os.path.join(checkpoint_output, "{chrom}.vcf.gz")).chrom)}

rule reseq_cat_ihs_bash:
	input:
		unpack(aggregate_ihs_reseq)
	output:
		scan_ihs=os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{config['species']}.ihs.out"),
		scan_log=os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{config['species']}.ihs.out.log"),
		norm_ihs=os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{config['species']}.ihs.norm"),
		norm_log=os.path.join(config['paths']['reseq_popgen_dir'], 'ihs', f"{config['species']}.ihs.norm.log")
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"""
		cat {input.scan_ihs} > {output.scan_ihs}
		cat {input.norm_ihs} > {output.norm_ihs}
		cat {input.scan_log} > {output.scan_log}
		cat {input.norm_log} > {output.norm_log}
		"""