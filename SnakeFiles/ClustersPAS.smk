configfile : "config.yaml"

def read_wildcards(accession):
    file=accession+'_STAR/'+accession+'Solo.out/Gene/filtered/'+accession+'_HieRFIT_cellTypes.list'
    with open(file) as f:
        return [line.strip().split("\t")[0] for line in f]

NBpredList=[]
for srr in config["accession"]:
    for cls in read_wildcards(srr):
        NBpredList.append(srr+'_STAR/'+srr+'_'+cls+'_Atrimmed_sorted_stranded_read_count.nbpred')

rule all:
    input:
        NBpredList

rule ProcessBam:
        input:
            bam="{accession}_STAR/{accession}_{cls}.bam"
        output:
            Pbam="{accession}_STAR/{accession}_{cls}_sorted.bam"
        shell:
            'samtools sort {input.bam} -o {output.Pbam} && \
             samtools index {output.Pbam} '

rule GenomeCov:
    """
    Creating a coverage data through the genome.
    """
        input:
            "{accession}_STAR/{accession}_{cls}_sorted.bam"
        output:
            "{accession}_STAR/{accession}_{cls}_sorted_stranded_read_count.bed"
        params:
            "{accession}_STAR/{accession}_{cls}_sorted_stranded_read_count.tmp"
        shell:
            """
            bash -c '
                source /usr/local/anaconda3/etc/profile.d/conda.sh
                conda activate CAAPS
                python ~/Dropbox/codes/endSeqTools/endSeq_tools.py genomeCov \
                --bam {input} --re 3
                '
             """

rule CorrectChrs:
    input:
        "{accession}_STAR/{accession}_{cls}_sorted_stranded_read_count.bed"
    output:
        "{accession}_STAR/{accession}_{cls}_Atrimmed_sorted_stranded_read_count.bed"
    shell:
        """
        awk '{{print "chr"$0}}' {input} > {output}
        """
rule NBClassifier:
    """
    Running a Naive Bayes classifier to evaluate Internal Priming Events (IPE).
    """
        input:
            "{accession}_STAR/{accession}_{cls}_Atrimmed_sorted_stranded_read_count.bed"
        output:
            "{accession}_STAR/{accession}_{cls}_Atrimmed_sorted_stranded_read_count.nbpred"
        shell:
            """
            bash -c '
                source /usr/local/anaconda3/etc/profile.d/conda.sh
                conda activate CAAPS
                Rscript --vanilla ~/Dropbox/codes/CAAPS/scripts/SetNBpred.1.0.R \
                {input} {output} human \
                '
            """


#snakemake --cores 2 -s ClustersPAS.smk --reason -n
#snakemake --cores 6 -s ClustersPAS.smk -n --dag | dot -Tpdf > dag.3.pdf