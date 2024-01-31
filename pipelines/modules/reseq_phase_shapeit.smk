rule all:
    input:
        os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.vcf.gz")

# Rule to phase the unphased VCF using Shapeit
rule shapeit_phase_reseq:
    input:
        os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['vcf_prefix']}.vcf.gz")
    output:
        temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.haps")),
        temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.sample")),
        temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.phase.snp.mm")),
        temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.phase.ind.mm")),
        log=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.phase.log")
    params:
        haps_prefix=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}")
    shell:
        "shapeit -V {input} -O {params.haps_prefix} --output-log {output.log}"

# Rule to convert the phased HAPS file back to VCF using Shapeit convert
rule shapeit_convert_reseq_to_vcf:
    input:
        os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.haps"),
        os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.sample")
    output:
        temp(os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.shapeit_header.vcf.gz")),
        log=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.convert.log")
    params:
        haps_prefix=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}"),
        shapeit_vcf=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.shapeit_header.vcf"),
    shell:
        """
        shapeit -convert --input-haps {params.haps_prefix} --output-vcf {params.shapeit_vcf} --output-log {output.log}
        bgzip {params.shapeit_vcf}
        """

# Rule to create extract the header from the original VCF
rule bcftools_create_header:
    input:
        os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['vcf_prefix']}.vcf.gz")
    output:
        temp(os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['vcf_prefix']}.header"))
    params:
        haps_output=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.haps"),
        shapeit_vcf=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.shapeit_header.vcf.gz")
    shell:
        """
        bcftools view -h {input} > {output}
        date=$(date +'%a %b %H:%M:%S %Y')
        insert_pos=$(wc -l < {output})
        awk -v n=$insert_pos -v s="##shapeit_covertCommand=-convert --input-haps {params.haps_output} --output-vcf {params.shapeit_vcf}; Date=$date" 'NR == n {{print s}} {{print}}' {output} > {output}.tmp && mv {output}.tmp {output}
        awk -v n=$insert_pos -v s="##shapeit_phaseCommand=-V {input} -O {params.haps_output};  Date=$date" 'NR == n {{print s}} {{print}}' {output} > {output}.tmp && mv {output}.tmp {output}
        awk -v n=$insert_pos -v s="##bcftools_pipelineCommand=view {params.shapeit_vcf} | reheader {output};  Date=$date" 'NR == n {{print s}} {{print}}' {output} > {output}.tmp && mv {output}.tmp {output}
        """

# Rule 2: Use the header file to replace the header of another VCF
rule bcftools_replace_shapit_header:
    input:
        shapeit_vcf=os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.shapeit_header.vcf.gz"),
        header=os.path.join(config['paths']['reseq_vcf_unphased_dir'], f"{config['vcf_prefix']}.header")
    output:
        os.path.join(config['paths']['reseq_vcf_phased_dir'], f"{config['vcf_prefix']}.vcf.gz")
    shell:
        "bcftools view {input.shapeit_vcf} | bcftools reheader -h {input.header} | bcftools view -O z -o {output}"