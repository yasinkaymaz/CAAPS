configfile : "/home/asus/CAAPS/config.yaml"

print (config['accession'])
print (config['read_no'])

#You can specify the wildcards with the method written below
wildcard_constraints:
    accession = config["accession"]


rule all:
    input:
        #["/home/asus/CAAPS/SRR8325947_1.fastq.gz",
         #"/home/asus/CAAPS/SRR8325947_2.fastq.gz"]
        expand("/home/asus/CAAPS/{accession}_{read_no}.fastq.gz", read_no= config["read_no"], accession= config["accession"])
'''
rule prefetch:
    output:
        "/home/asus/CAAPS/{accession}.sra"
    params:
        "{accession} --max-size 50GB -O 01_raw"
    log:
        "/home/asus/CAAPS/{accession}.log"
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
        "/home/asus/CAAPS/{accession}_1.fastq.gz",
        "/home/asus/CAAPS/{accession}_2.fastq.gz"


    params:
        conda_env = config['conda_env'],
        args = "--split-files --threads 4 --mem 2048MB --progress --details",
    log:
        "/home/asus/CAAPS/{accession}.log"

    shell:
        """
        fasterq-dump {params.args} {accession}  > {output}
        """
