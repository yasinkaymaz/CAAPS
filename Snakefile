configfile : "config.yaml"
accession= config["accession"]

print("Accession number is:", accession)

rule all:
    input: 
        expand("~/CAAPS/{accession}_1.fastq.gz", accession= config["accession"]),
        expand("~/CAAPS/{accession}_2.fastq.gz", accession= config["accession"]),
        expand("{accession}_STAR/", accession= config["accession"])
        	
rule get_SRA_by_accession:
    """
    Retrieve FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "~/CAAPS/{accession}_1.fastq.gz",
        "~/CAAPS/{accession}_2.fastq.gz"
    params:
        args = "-S --progress"
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
            fa = '~/CAAPS/genomefile/genome.fa', # provide human genome FASTA
            gtf = '~/CAAPS/genomefile/genes.gtf' # provide human genome GTF
        output:
            directory("CAAPS/humanINDEX") # Human Index Folder
        threads: 20 # set the maximum number of available cores
        shell:
            'mkdir {output} && '
            'STAR --runThreadN {threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--genomeFastaFiles {input.fa} '
            '--sjdbGTFfile {input.gtf} '
            '--sjdbOverhang 100'
# a pre-made genome can be found in https://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/
#It needs STAR 2.7.1 to be run
rule STAR:
    """
    STARsolo alignment
    """
        input:
            R1="~/CAAPS/{accession}_1.fastq.gz",
            R2="~/CAAPS/{accession}_2.fastq.gz",
            refdir="~/CAAPS/humanINDEX"
        params:
            outdir = '{accession}_STAR',
            whitelist = '~/CAAPS/10x_V3_whitelist.txt',
            prefix = 's1'
        output:
            '{accession}_STAR/'
        threads: 20 # set the maximum number of available cores
        
        shell:
            'mkdir {params.outdir} && ' # snakemake had problems finding output files with --outFileNamePrefix, so I used this approach instead
            'cd {params.outdir} && '
            'STAR --genomeDir {input.refdir} \
                --soloType CB_UMI_Simple \
                --soloCBwhitelist {params.whitelist} \
                --soloUMIlen 12 \
                --readFilesIn {input.R1}, {input.R2} \
                --runThreadN {threads} \
                --outFileNamePrefix {params.prefix} \
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
