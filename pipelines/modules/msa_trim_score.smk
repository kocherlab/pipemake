rule all:
    input:
        expand("Stats/MSA/{sample}.trim_fraction.txt", sample=config["samples"]),


rule clipkit_trim_fraction:
    input:
        "logs/clipkit/{sample}.log",
    output:
        "Stats/MSA/{sample}.trim_fraction.txt",
    log:
        "logs/clipkit/{sample}.trim_fraction.log",
    run:
        import re

        with open(log[0], "w") as logf:
            with open(input[0]) as fh:
                content = fh.read()

            match = re.search(
                r"Percentage of alignment trimmed:\s*([0-9]+\.?[0-9]*)%",
                content,
            )
            if not match:
                msg = f"Could not find 'Percentage of alignment trimmed' in {input[0]}"
                logf.write(msg + "\n")
                raise ValueError(msg)

            fraction = float(match.group(1)) / 100
            logf.write(f"Parsed percentage: {match.group(1)}%\n")
            logf.write(f"Converted fraction: {fraction:.2f}\n")

            with open(output[0], "w") as out:
                out.write(f"{fraction:.2f}\n")
