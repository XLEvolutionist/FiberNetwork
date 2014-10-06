FiberNetwork
============

WGCNA analysis using Fiber RNA seq data.

#Puropse

This set of scripts takes in RNA seq expression data and constructs co-expression networks from the data. When using these scripts you will need to modify several lines in order that your custom code will work. As such this is more of a record of the processes and code used to generate the data associated with this project and these set of scripts are meant to serve as a starting spot fro your own analysis, not a distribution ready set of scripts.

Brief details of each script can are provided in the file descriptions but it is worth mentioning the logic order in which these operate:

#included scripts and brief descriptions.

[1] Fiber_Norm_RPM.R will take in raw count data and output normalized libraries ans well as library sizes.

[2] WGCNA_construct.R will, trim the gene set. generate the co-expression object and into a directory (you can specify this by modifying the script.)

[3] WGCNA_analysis.R will perform some higher level analysis in order to visualize similarities and differences between co-expression networks.

[4] construct_igraph_networks.R will constrcut iGraph objects that can be plotted using the "plotg.R" script.

[5] edgeList.pl converts a n x n adjecency matrix (where n is the number of genes in the co-expression network) into edgeList format that igraph can understand.
