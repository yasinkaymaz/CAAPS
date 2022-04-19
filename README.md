# CAAPS

Current version can:
  -download from sra
  -download STAR_index from dobin's repo
  -Align using STAR (needs >30GB of RAM)

needs conda and snakemake to run

<code> conda install -c bioconda snakemake </code>

to run

<code> snakemake --use-conda -c all </code>

you can adjust the config file to change the accession entries in the given format
<code>
accession:
- SRR1
- SRR2
- etc.
</code>

Currently maintained by BMGlab
