configfile: "include/rules/config.yaml"

include: "include/rules/fasterqdump.rule"

rule all:
    input:
        expand("01_raw/done__{srr}_dump", srr=config['srr'])


rule prefetch:
    output:
        "01_raw/.prefetch/sra/{srr}.sra"
    params:
        "{srr} --max-size 50GB -O 01_raw"
    log:
        "01_raw/.prefetch/sra/{srr}.log"
    conda:
        "yamls/sra-tools.yaml"
    shell:
        """
        prefetch {params} > {log} 2>&1 && touch {output}
        """

rule fastqdump:
    input:
        "01_raw/.prefetch/sra/{srr}.sra"
    output:
        touch("01_raw/done__{srr}_dump")
    params:
        args = "-S -O 01_raw/ -t 01_raw/",
        id_srr = "{srr}"
    log:
        "01_raw/{srr}.log"
    conda:
        "yamls/sra-tools.yaml"
    shell:
        """
        fasterq-dump {params.args} {params.id_srr} > {log} 2>&1
        """
