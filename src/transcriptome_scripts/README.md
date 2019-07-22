This folder contains all local scripts written for the analysis of the C. virginica OA exposure experiment. 

Folders:
* `STAR_scripts` : Scripts specifically used for performing the mapping and counting (via RSEM) of the RNA-seq data using the STAR mapping software and RSEM as the primary transcript quanfication approach.
* `salmon_scripts` : Scripts used for the salmon, "reference-free", approach.
	R : R scripts specifically for use on the salmon outputs.
* `software` : Contains the various software used for analysis that I couldn't get installed on the shared software section of the server. These include:
	RSEM :  probabilistic transcript quantification approach used in conjunction with STAR for determine gene and transcript (isoform) level counts.
	gffread : allows for the converstion from gff gene annotation files to gtf
	gffcompare : allows for comparisons between file formats
