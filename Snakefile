configfile : "~/CAAPS/config.yaml"

print (config['accession'])
print (config['read_no'])

#You can specify the wildcards with the method written below
wildcard_constraints:
    accession = config["accession"]


rule all:
    input:
        #["~/CAAPS/SRR8325947_1.fastq.gz",
         #"~/CAAPS/SRR8325947_2.fastq.gz"]
        expand("~/CAAPS/{accession}_{read_no}.fastq.gz", read_no= config["read_no"], accession= config["accession"])
'''
rule prefetch:
    output:
        "~/CAAPS/{accession}.sra"
    params:
        "{accession} --max-size 50GB -O 01_raw"
    log:
        "~/CAAPS/{accession}.log"
    shell:
        """
        prefetch {params} > {log} 2>&1 && touch {output}
        """

'''

accession = config['accession']

rule get_SRA_by_accession:
    """
    Retrieve FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "~/CAAPS/{accession}_1.fastq.gz",
        "~/CAAPS/{accession}_2.fastq.gz"
    params:
        conda_env = config['conda_env'],
        args = "--split-files --threads 4 --mem 2048MB --progress --details",
    log:
        "~/CAAPS/{accession}.log"
    shell:
        """
        fasterq-dump {params.args} {accession}  > {output}
        """

rule STARindex:
    """
    Creating a index of human genome to align
    we can later on put a wget link for downloading an already existing one
    https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
    refdata-gex-GRCh38-2020-A/genes/genes.gtf
    refdata-gex-GRCh38-2020-A/fasta/genome.fa
    10x uses a filtered gtf, need to download from 10x or ->
    #
    If you want to use your own GTF (e.g. newer version of ENSEMBL or GENCODE), you can generate the "filtered" GTF file using 10X's mkref tool:
    https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
    #
    """
        input:
            fa = '~/CAAPS/genome.fa', # provide human genome FASTA
            gtf = '~/CAAPS/genes.gtf' # provide human genome GTF
        output:
            directory('humanINDEX') # Human Index Folder
        threads: 20 # set the maximum number of available cores
        shell:
            'mkdir {output} && '
            'STAR --runThreadN {threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--genomeFastaFiles {input.fa} '
            '--sjdbGTFfile {input.gtf} '
            '--sjdbOverhang 100'

rule STAR:
    """
    STARsolo alignment
    """
        input:
            R1="~/CAAPS/{accession}_1.fastq.gz",
            R2="~/CAAPS/{accession}_2.fastq.gz"
            refdir = directory('humanINDEX')
        params:
            outdir = '{sample}_STAR'
            whitelist = '~/CAAPS/10x_V3_whitelist.txt'
        output:
            '{sample}_STAR/'
        threads: 20 # set the maximum number of available cores
        prefix: 's1'
        shell:
            'mkdir {params.outdir} && ' # snakemake had problems finding output files with --outFileNamePrefix, so I used this approach instead
            'cd {params.outdir} && '
            'STAR --genomeDir {input.refdir} \
                --soloType CB_UMI_Simple \
                --soloCBwhitelist {params.whitelist} \
                --soloUMIlen 12 \
                --readFilesIn {input.R1}, {input.R2} \
                --runThreadN {threads} \
                --outFileNamePrefix {prefix} \
                --outSAMtype BAM SortedByCoordinate \
                --outReadsUnmapped elp1_s2_Unmapped \
                --twopassMode Basic \
                --chimSegmentMin 20 \
                --readFilesCommand zcat \
                --clipAdapterType CellRanger4 \
                --outFilterScoreMin 30 \
                --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
                --soloUMIfiltering MultiGeneUMI_CR \
                --soloUMIdedup 1MM_CR '
