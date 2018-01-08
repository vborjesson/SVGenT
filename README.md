# SVGenT
Tool for better genotyping of whole genome sequencing (WGS) data. 

![alt text](https://github.com/vborjesson/SVGenT/blob/master/SVGenT.png)

Takes an vcf, bam and tab-file with read depth per 100 bp as input. 

dependencies: 

```
SQLite3
ABySS
SSAKE
BWA
samtools
clustalw
SVGenT.db (look futher down in this readme)
```

## Create database
requirements; SQLite3
Statistics are needed to confirm and support the new predicted genotype and SV-type, therefor SVGenT requires SVGenT.db with GC-content and mappability information. Download mappability scores for reference genome hg19 from rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability and convert the bigwig file to readable bedGraph. Do also download the hg19 reference genome. 

```
python SVGenT_DB.py --hg19 /path/to/refgenome --bed /path/to/bedGraph
```
This will create SVGenT.db SQL database in the SVGenT directory, takes about 2 hours.

## Run SVGenT
```
python SVGenT.py -h 

usage: SVGenT.py [-h] --vcf VCF_IN --bam BAM_IN --tab TAB_IN [--ID ID]
                 [--sam SAM] --bwa_ref BWA_REF

SVGenT takes vcf- and bam-files as input and improve the prediction of
genotype

optional arguments:
  -h, --help         show this help message and exit
  --vcf VCF_IN       Path to vcf-file
  --bam BAM_IN       Path to bam-file
  --tab TAB_IN       Path to tab_file
  --ID ID            sample ID
  --sam SAM          path to sam-file for dry run
  --bwa_ref BWA_REF  Path to reference genome for bwa mem
```
An example:
```
python SVGenT.py --vcf --bam --tab --bwa_ref --ID
```
This will generate a vcf file with new predicted breakpoints, genotyping predictions. SVGenT updates ~ 700 SVs/hour, plus 30 min for developing the database. 





