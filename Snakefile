configfile : "config.yaml"

rule all:
    input:
        expand("{accession}_reads/{accession}_1.fastq", accession= config["accession"]),
        expand("{accession}_reads/{accession}_2.fastq", accession= config["accession"]),
        expand("{accession}_STAR/{accession}_Log.final.out", accession= config["accession"])

rule get_SRA_by_accession:
    """
    Retrieve FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "{accession}_reads/{accession}_1.fastq",
        "{accession}_reads/{accession}_2.fastq"
    params:
        args = "--split-files --progress --details",
        accession = "{accession}"
    log:
        "{accession}_reads/{accession}.log"
    conda:
        "sra_env.yml"
    shell:
        'mkdir -p {params.accession}_reads && '
        'fasterq-dump {params.args} {params.accession} -O {params.accession}_reads'


rule STARindex_Down:
        output:
            "humanINDEX/Genome",
            "humanINDEX/Homo_sapiens.GRCh38.99.gtf",
            "humanINDEX/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
            "humanINDEX/SA",
            directory("humanINDEX")
        conda:
            "sra_env.yml"
        log:
            "wget.log"
        shell:
            'mkdir {output} && \
            wget -r -np -nH -nc --cut-dirs=10 -P {output} https://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/ > {log} 2>&1 |& tee {log}'


rule STARsolo:
    """
    STARsolo alignment
    """
        input:
            R1="{accession}_reads/{accession}_1.fastq",
            R2="{accession}_reads/{accession}_2.fastq",
            refdir="humanINDEX/"
        params:
            outdir = "{accession}_STAR",
            prefix = "{accession}_"
            #extra = "{extra_star_params}"
        output:
            '{accession}_STAR/{accession}_Log.final.out'
        threads:
            20 # set the maximum number of available cores
        log:
            "{accession}_STAR/STARsolo_{accession}.log"
        conda:
            "star_env.yml"
        shell:
            'mkdir -p {params.outdir} && '
            'STAR --genomeDir {input.refdir} \
            --readFilesIn {input.R1},{input.R2} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.outdir}/{params.prefix}'
