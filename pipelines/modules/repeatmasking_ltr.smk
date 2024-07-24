rule all:
	input:
		os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa.masked")

rule repeat_modeler_r1:
	input:
		os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa")
	output:
		os.path.join(config['paths']['repeatmodeler_dir'], 'R1', 'Families', f"{config['species']}_{config['assembly_version']}-families.fa")
	params:
		repeatmodeler_db=os.path.join(config['paths']['repeatmodeler_dir'], 'R1', 'DB', f"{config['species']}_{config['assembly_version']}_DB"),
		db_dir=os.path.join(config['paths']['repeatmodeler_dir'], 'R1', 'DB'),
		wd_dir=os.path.join(config['paths']['repeatmodeler_dir'], 'R1', 'WorkingDirectory')
	threads: 20
	singularity: 
		"docker://dfam/tetools:v1.88.5"
	shell:
		r"""
		mkdir {params.db_dir} &&
		mkdir {params.wd_dir} &&
		BuildDatabase -name {params.repeatmodeler_db} {input} &&
		RepeatModeler -database {params.repeatmodeler_db} -threads {threads} &&
		rm_working_dir=$(grep 'Using output directory' {params.repeatmodeler_db}-rmod.log | grep -o '[^/]*$') &&
		tar -czf {params.wd_dir}/RepeatModeler_R1_WD.tar.gz $rm_working_dir &&
		rm -rf $rm_working_dir &&
		mv {params.repeatmodeler_db}-families.fa {output} &&
		mv {params.repeatmodeler_db}-families.stk {params.wd_dir} &&
		mv {params.repeatmodeler_db}-rmod.log {params.wd_dir}
		"""

rule repeat_masker_r1:
	input:
		assembly=os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa"),
		families=os.path.join(config['paths']['repeatmodeler_dir'], 'R1', 'Families', f"{config['species']}_{config['assembly_version']}-families.fa")
	output:
		os.path.join(config['paths']['repeatmodeler_dir'], 'R1', 'MaskedAssembly', f"{config['species']}_{config['assembly_version']}.fa.masked")
	params:
		mask_dir=os.path.join(config['paths']['repeatmodeler_dir'], 'R1', 'MaskedAssembly')
	threads: 20
	singularity: 
		"docker://dfam/tetools:v1.88.5"
	shell:
		"RepeatMasker -par {threads} -dir {params.mask_dir} -lib {input.families} {input.assembly}"

rule repeat_modeler_r2:
	input:
		os.path.join(config['paths']['repeatmodeler_dir'], 'R1', 'MaskedAssembly', f"{config['species']}_{config['assembly_version']}.fa.masked")
	output:
		directory(os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'DB')),
		directory(os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'WorkingDirectory')),
		os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'Families', f"{config['species']}_{config['assembly_version']}-families.fa")
	params:
		repeatmodeler_db=os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'DB', f"{config['species']}_{config['assembly_version']}_DB"),
		db_dir=os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'DB'),
		wd_dir=os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'WorkingDirectory')
	threads: 20
	singularity: 
		"docker://dfam/tetools:v1.88.5"
	shell:
		r"""
		mkdir {params.db_dir} &&
		mkdir {params.wd_dir} &&
		BuildDatabase -name {params.repeatmodeler_db} {input} && 
		RepeatModeler -database {params.repeatmodeler_db} -threads {threads} -LTRStruct &&
		rm_working_dir=$(grep 'Using output directory' {params.repeatmodeler_db}-rmod.log | grep -o '[^/]*$') &&
		tar -czf {params.wd_dir}/RepeatModeler_R1_WD.tar.gz $rm_working_dir &&
		rm -rf $rm_working_dir &&
		mv {params.repeatmodeler_db}-families.fa {output} &&
		mv {params.repeatmodeler_db}-families.stk {params.wd_dir} &&
		mv {params.repeatmodeler_db}-rmod.log {params.wd_dir}
		"""

rule repeat_masker_r2:
	input:
		assembly=os.path.join(config['paths']['repeatmodeler_dir'], 'R1', 'MaskedAssembly', f"{config['species']}_{config['assembly_version']}.fa.masked"),
		families=os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'Families', f"{config['species']}_{config['assembly_version']}-families.fa")
	output:
		os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'MaskedAssembly', f"{config['species']}_{config['assembly_version']}.fa.masked.masked")
	params:
		mask_dir=os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'MaskedAssembly')
	threads: 20
	singularity: 
		"docker://dfam/tetools:v1.88.5"
	shell:
		"RepeatMasker -par {threads} -dir {params.mask_dir} -lib {input.families} {input.assembly}"

rule softmask_r2:
	input:
		assembly=os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa"),
		masked_assembly=os.path.join(config['paths']['repeatmodeler_dir'], 'R2', 'MaskedAssembly', f"{config['species']}_{config['assembly_version']}.fa.masked.masked")
	output:
		os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa.masked")
	singularity: 
		"docker://aewebb/pipemake_utils:v0.1.27"
	shell:
		"softmask.py {input.assembly} {input.masked_assembly} {output}"