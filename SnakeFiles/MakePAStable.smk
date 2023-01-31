configfile : "config.yaml"
valueType = ['Count', 'Percent', 'CPM']

PasTableList=[]
for srr in config["accession"]:
    for val in valueType:
        PasTableList.append(srr+"_STAR/"+srr+"_mergedPeaks_annotated_"+val+".bed")
print(PasTableList)

rule all:
    input:
        PasTableList
        #expand("{accession}_STAR/{accession}_mergedPeaks_annotated_{val}.bed", accession= config["accession"], val=valueType)


rule dumpSamplesList:
    input:
        files="{accession}_STAR/{accession}_{cls}_Atrimmed_sorted_stranded_read_count.nbpred"
        #lambda wildcards: "{accession}_STAR/{accession}Aligned.{cls}_stranded_read_count.nbpred".format(cls=wildcards.cls)
        #expand("{accession}_STAR/{accession}Aligned.{cls}_stranded_read_count.nbpred", accession= config["accession"], cls=Cluster_outs)
    output:
        "{accession}_STAR/{accession}_{cls}.fastq"
    params:
        "master_samples.list.txt"
    shell:
        """
        touch {output}
        for item in {output};
            do
            echo $item >> {params};
        done
        """


#files = [f'{gene}/assoc/{marker}' for gene in genes for marker in markers[gene]]
#print(config["accession"])
#print([read_wildcards(x) for x in config["accession"] ])

rule MakePAStable:
    """
    Combine all cell type PAS data in a table.
    """
        input:
            #PAStablelist
            SamplesList= "master_samples.list.txt"
            #lambda wildcards: "{accession}_STAR/{accession}Aligned.{cls}_stranded_read_count.nbpred".format(cls=read_wildcards(wildcards.accession))
            #"{accession}_STAR/{accession}Aligned.{cls}_stranded_read_count.bed"
            #lambda wildcards: "{accession}_STAR/{accession}Aligned.{cls}_stranded_read_count.bed".format(cls=read_wildcards(wildcards.accession), accession=wildcards.accession)
            #nbpreds=lambda wildcards: "{accession}_STAR/{accession}Aligned.{cls}_stranded_read_count.nbpred".format(cls=read_wildcards(wildcards.accession), accession=config["accession"])
        output:
            "{accession}_STAR/{accession}_mergedPeaks_annotated_{val}.bed"
            #cntTab=expand("{accession}_STAR/{accession}_mergedPeaks_annotated_{val}.bed", accession= config["accession"], val=valueType)
        params:
            SamplesList= "master_samples.list.txt",
            GTFfile=config["gftFile"],
            endSeqToolsDir=config["endDir"]
        shell:
            """
            bash -c '
                source /usr/local/anaconda3/etc/profile.d/conda.sh
                conda activate CAAPS
                bash ~/Dropbox/codes/CAAPS/scripts/make_PAS_table.v3.sh {params.SamplesList} {params.GTFfile} {params.endSeqToolsDir}
                '
            """
            

#snakemake --cores 2 -s ClustersPAS.smk --reason -n
#snakemake --cores 6 -s ClustersPAS.smk -n --dag | dot -Tpdf > dag.4.pdf