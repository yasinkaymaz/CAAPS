# CAAPS

Current version can:
  -download from sra
  -download STAR_index from dobin's repo
  -Align using STAR (needs >30GB of RAM)
  
Curre

needs conda and snakemake to run
<code> conda install -c bioconda snakemake </code>

to run 
<code> snakemake --use-conda -c all </code>

you can adjust the config file to change the accession entries in the given format

<code> <pre>
accession:
- SRR1
- SRR2
- etc.
</pre>
</code>

Currently maintained by BMGlab
