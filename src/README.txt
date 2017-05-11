-----------------------
| CONVERGENCE PACKAGE |
-----------------------

--------------------------
REQUIREMENT

	The package needs the gsl and the nlopt libraries.

--------------------------
COMPILING

	Just type
	> make all
	to build all the binaries.

--------------------------
DESCRIPTION

	'msd' detects molecular signatures of convergence events by using the convergence index.
	
	'enr' computes GO terms enrichment of a list of genes.

--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	msd - detection of molecular signatures of convergence events
	
SYNOPSIS
	msd [OPTIONS] [FILE TREE] [FILE CONV.] [ALIGN. FILE 1] [ALIGN. FILE 2] ...

DESCRIPTION
	return a table where each line displays the result of the detection of molecular signature
	of each alignment files [ALIGN. FILE i] in Fasta format, with the regard to the tree in Newick format of 
	[FILE TREE] and the character stored as a table in [FILE CONV.].
	
	Options are
	
	-n [FILE]
		set the optimisation options to those contained in [FILE].
	-m [FILE]
		set the evolution model to that contained in [FILE].
	-o [NAME]
		set the name of the output file to [NAME]. By default, it is 'output.txt'.
	-e [VALUE]
		set the threshold under which a site is considered significant.
	-e [VALUE]
		save an extra file named 'significant.txt', which contains all the names of alignments
		with a corrected p-value under [VALUE].
	-t [NUMBER]
		set the maximal number of simultaneous threads to [NUMBER].
	-h
		display help

-------------------------

NAME
	enr - computing the GO terms enrichments of a list of genes.
	
SYNOPSIS
	enr [OPTIONS] [FILE GO] [FILE LIST]
	
DESCRIPTION
	save a file containg the GO terms associated to at least a gene of [FILE LIST] according to [FILE GO].
	
	-o [FILE]
		load the GO descriptors from [FILE].
	-b [FILE LIST]
		set the background set of genes to those of [FILE LIST].
	-h
		display help
