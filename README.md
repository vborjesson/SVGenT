# SVGenT
Tool for better genotyping of whole genome sequencing (WGS) data. 

![alt text](https://github.com/vborjesson/SVGenT/blob/master/SVGenT.png)

Takes an vcf, bam and tab-file with read depth per 100 bp as input. 

dependencies: 

```
ABySS
SSAKE
velvet
BWA
samtools
clustalw
SVGenT.db (look futher down in this readme)
consensus (https://github.com/J35P312/SplitVision) 
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
python SVGenT.py --vcf --bam --tab --bwa_ref --ID
```
This will generate a vcf file with new predicted breakpoints, genotyping predictions.





