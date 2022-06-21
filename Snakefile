configfile : "config.yaml"

rule all:
    input:
        expand("{accession}_reads/{accession}_1.fastq", accession= config["accession"]),
        expand("{accession}_reads/{accession}_2.fastq", accession= config["accession"]),
        expand("{accession}_STAR/", accession= config["accession"])

rule get_SRA_by_accession:
    """
    Retrieve FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "{accession}_reads/{accession}_1.fastq",
        "{accession}_reads/{accession}_2.fastq"
    params:
        args = "--split-files --progress",
        accession = "{accession}"
    log:
        "{accession}_reads/{accession}.log"
    conda:
        "caaps_env.yml"
    shell:
        'mkdir -p {params.accession}_reads && '
        'fasterq-dump {params.args} {params.accession} -O reads > {log} 2>&1 |& tee {log}'


rule STARindex_Down:
        output:
            "humanINDEX/Genome",
            "humanINDEX/Homo_sapiens.GRCh38.99.gtf",
            "humanINDEX/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
            "humanINDEX/SA",
            directory("humanINDEX")
        conda:
            "caaps_env.yml"
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
            prefix = "{accession}"
            #extra = "{extra_star_params}"
        output:
            '{accession}_STAR/'
        threads:
            20 # set the maximum number of available cores
        log:
            "{accession}_STAR/STARsolo_{accession}.log"
        conda:
            "caaps_env.yml"
        shell:
            'mkdir -p {params.outdir} && '
            'STAR --genomeDir {input.refdir} \
            --readFilesIn {input.R1},{input.R2} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.outdir}/{params.prefix} > {log} 2>&1 |& tee {log}'

rule UMI_dedup:
        '''
PCR duplicates Removal, from scAPA;
"PCR duplicates are removed using UMI-tools dedup. As UMI tools dedup
requires that each line in the BAM file has a molecular barcode tag, Dropseq tools is first used to filter the BAMs, leaving only reads for which cell ranger counts produced corrected molecular barcode tag.Then UMI tools is ran with “method=unique” so that cellranger corrected
molecular barcodes are used"
	'''
        input:
            Bam="{accession}_STAR/Aligned.bam"
        params:
            outdir = "{accession}_dedup",
            prefix = "{accession}"
            #extra = "{extra_star_params}"
        output:
            'dedup.Aligned.bam'
        threads:
            20 # set the maximum number of available cores
        log:
            "{accession}_dedup/UMItools_{accession}.log"
        conda:
            "caaps_env.yml"
        shell:
            'mkdir -p {params.outdir} && '
            'FilterBAM TAG_RETAIN=UB I={input.bam} O= UB.Aligned.bam && \
            umi_tools dedup -I UB.Aligned.bam -S {output} --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB'
