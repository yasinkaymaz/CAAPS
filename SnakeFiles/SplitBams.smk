configfile : "config.yaml"

rule all:
    input:
        expand("{accession}_STAR/sinto_{accession}.log", accession=config["accession"])


rule ProcessBam:
        input:
            bam='{accession}_STAR/{accession}Aligned.out.bam'
        output:
            Pbam='{accession}_STAR/{accession}Aligned.sorted.bam'
        shell:
            'samtools sort {input.bam} -o {output.Pbam} && \
             samtools index {output.Pbam} '


rule SplitBams2Clusters:
    """
    Spliting read alignments based on cell type clusters.
    """
        input:
            CellClusters='{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_HieRFIT_cellTypes.txt',
            CellClustersList='{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_HieRFIT_cellTypes.list',
            Pbam='{accession}_STAR/{accession}Aligned.sorted.bam'
        output:
            #expand("{accession}_STAR/{accession}Aligned.{cls}.bam", accession= config["accession"], cls=read_wildcards('{accession}_STAR/{accession}Solo.out/Gene/filtered/{accession}_seurat_clusters.txt'))
            "{accession}_STAR/sinto_{accession}.log"
        params:
            runDir='{accession}_STAR/',
            srrids='{accession}'
        threads:
            2
        shell:
            """
            bash -c '
                sinto filterbarcodes -b {input.Pbam} \
                -c {input.CellClusters} \
                --barcodetag "CR" \
                --outdir {params.runDir} \
                -p {threads}
                
                for i in `cat {input.CellClustersList}`;
                    do
                    mv {params.runDir}/$i.bam {params.runDir}/{params.srrids}_$i.bam;
                done

                touch {output}
                '
            """


#snakemake --cores 2 -s SplitBams.smk --reason -n
#snakemake --cores 2 -s SplitBams.smk -n --dag | dot -Tpdf > dag.2.pdf