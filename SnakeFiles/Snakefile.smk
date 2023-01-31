configfile : "config.yaml"
valueType = ['Count', 'Percent', 'CPM']

rule all:
    input:
#         #expand("{accession}_reads/{accession}_1.fastq.gz", accession= config["accession"]),
#         #expand("{accession}_reads/{accession}_2.fastq.gz", accession= config["accession"]),
#         #expand("{accession}_STAR/{accession}Aligned.out.bam", accession= config["accession"]),
        expand('{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_HieRFIT_cellTypes.list', accession= config["accession"])

rule get_SRA_by_accession:
    """
    Retrieve FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        '{accession}_reads/{accession}_1.fastq.gz',
        '{accession}_reads/{accession}_2.fastq.gz',
    params:
        args = "--split-files --progress",
        accession = "{accession}"
    log:
        "{accession}_reads/{accession}.log"
#    conda:
#        "caaps_env.yml"
    shell:
        #'mkdir -p {params.accession}_reads && '
        'fasterq-dump {params.args} {params.accession} -O {params.accession}_reads > {log} 2>&1 | tee {log} && \
        gzip {params.accession}_reads/*'


rule STARindex_create:
        output:
            "humanSTARindex/Genome",
            "humanSTARindex/Homo_sapiens.GRCh38.99.gtf",
            "humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
            "humanSTARindex/SA",
            directory("humanSTARindex")
#        conda:
#            "caaps_env.yml"
        log:
            "wget.log"
        threads:
            4
        shell:
            #'mkdir {output} && \
            'STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir humanSTARindex \
            --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
            --sjdbGTFfile Homo_sapiens.GRCh38.99.gtf'

rule STARsolo:
    """
    STARsolo alignment
    """
        input:
            R1="{accession}_reads/{accession}_1.fastq.gz",
            R2="{accession}_reads/{accession}_2.fastq.gz",
            refdir="humanSTARindex/"
        params:
            outdir = "{accession}_STAR",
            prefix = "{accession}"
        output:
            temp('{accession}_STAR/{accession}Aligned.out.bam'),
            '{accession}_STAR/{accession}Solo.out/Gene/filtered/matrix.mtx'
        threads:
            2 # set the maximum number of available cores
        log:
            "{accession}_STAR/STARsolo_{accession}.log"
        shell:
            'STAR --genomeDir {input.refdir} \
            --readFilesIn {input.R2} {input.R1} \
            --readFilesCommand gunzip -c \
            --soloType Droplet \
            --outSAMattributes CR CY UR UY \
            --outSAMtype BAM Unsorted \
            --runThreadN {threads} \
            --soloCBwhitelist {input.refdir}/737K-august-2016.txt \
            --outFileNamePrefix {params.outdir}/{params.prefix} > {log} 2>&1 | tee {log}'

rule GzipMtx:
    """
    Zip to prepare for Seurat.
    """
        input:
            '{accession}_STAR/{accession}Solo.out/Gene/filtered/matrix.mtx'
        output:
            '{accession}_STAR/{accession}Solo.out/Gene/filtered/matrix.mtx.gz'
        params:
            prefix = "{accession}"
        shell:
            'gzip {params.prefix}_STAR/{params.prefix}Solo.out/Gene/filtered/*'


rule StandardSeurat:
    """
    To process the count matrix and obtain cell clusters.
    """
        input:
            StarMatrix='{accession}_STAR/{accession}Solo.out/Gene/filtered/matrix.mtx.gz'
        output:
            seuratObj='{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_seurat.rds',
            CellClusters= '{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_seurat_clusters.txt'
        params:
            RunDir='{accession}_STAR/{accession}Solo.out/Gene/filtered/',
            RunName='{accession}'
        shell:
            """
            bash -c '
                source /usr/local/anaconda3/etc/profile.d/conda.sh
                conda activate CAAPS
                Rscript ~/Dropbox/codes/CAAPS/scripts/MakeSeuratObject.R {params.RunDir} {params.RunName}
                '
            """

rule HieRFIT:
    input:
        seuratObj='{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_seurat.rds'
    output:
        CellTypes='{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_HieRFIT_cellTypes.txt'
    params:
        HieRmod="cancerHiermod.hm.RDS"
    shell:
        """
        bash -c '
            source /usr/local/anaconda3/etc/profile.d/conda.sh
            conda activate CAAPS
            Rscript ~/Dropbox/codes/CAAPS/scripts/HieRFIT.R {params.HieRmod} {input.seuratObj} {output}
            '
            """
        
rule CreateCellTypeList:
    input:
        CellClusters='{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_HieRFIT_cellTypes.txt'
    output:
        CellClustersList='{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_HieRFIT_cellTypes.list'
    shell:
        "cut -f2 {input.CellClusters}|sort |uniq > {output.CellClustersList}"

# rule CreateClustersList:
#     input:
#         CellClusters='{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_seurat_clusters.txt'
#     output:
#         CellClustersList='{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_seurat_clusters.list'
#     shell:
#         "cut -f2 {input.CellClusters}|sort |uniq > {output.CellClustersList}"


#snakemake --cores 2 -s SnakeFile.smk --reason -n
#snakemake --cores 6 -s SnakeFile.smk -n --dag | dot -Tpdf > dag.1.pdf