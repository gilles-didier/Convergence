Repertory 'results' contains 3 types of files for all pairs of species of the dataset:

1) Files 'output....csv' are tables with columns 
	#name: the ident of the gene,
	status: signals a possible issue in the p-value computation (see 'warnings..'),
	length: the length of the gene
	unknown percent: the percentage of wildcards or gap in the alignment,
	pvalue brute: binomial p-value,
	pvalue corrected: the preceding p-value corrected by using Benjamini-Hochberg,
	number of convergent sites: 
	convergent sites: the position of the convergent sites.
	
2) Files 'significant....txt' contains the list of the significant genes (i.e. here with a corrected p-value lower than 5%).

3) Files 'warnings....txt' contains the level of confidence of the p-values with a possible issue (cf column 'status' with '*' of the corresponding 'output...' files - all other p-values can be trust at a confidence level of 99.99%)

Results were obtained with a significance threshold for the site (gamma) 0.0001
