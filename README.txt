README

Sannie Loi Bioinformatics Tool for motif function prediction:

1. Divide protein into motif fragments of length N 
2. Search Arabidopsis proteins for matches (allowing 1 mismatch)
3. Test identified proteins for enriched functions (i.e., GO terms) using hypergeometric test
4. Produce list of top functional enrichments for each motif
5. Highlight the most commonly occurring functional enrichments throughout the entire protein
6. Cluster results to remove redundancy
7. Adjust p-values for multiple hypothesis testing (Bonferroni)

Currently I have separate scripts for the steps inorder for easier testing, but after I optimize the program I plan to combine them.

Steps 1&2 --> splitMotif.py
- Function: splitMotif(size of sliding window, user_file, outputfile_name)
- the userfile is the sequence where there are no functions associated with 

Steps 3 --> hypergeometricTest.py
- takes the output from steps 1&2 and runs it through a hypergeometric test. Loops through all the possible terms and motifs, so takes a while to run.
-Function: hypergeometric(userfile, output)

Steps 4 --> frequencyTable.py
- Given a cutoff value (i.e. Top 50), the frequency table produces a list of terms ordered by the most appeared to least within the cutoff range
- Function: freqtable(Output_from_step3_file, output_file, cutoff_value)

Steps 5,6,&7 --> plotMotifComplete.py
Plot motif plots the -log(pvalue) of the terms(provided by user) and the user has the option to either have the graphs as overlays (True) or as subplots (False) given a list of terms, a cutoff value, and a coefficient value. The p-values were corrected with multiple-hypothesis (MH) testing adjustment by multiplying by the # of motifs and # of terms. 
Overlay
-terms that are greater than the given coefficient value are shown in grey
-terms that are less than the given coefficient value are shown in different colours
Subplots
-terms with less than the given coefficient values each have their own subplot
-terms that are greater than the given coefficient values are shown in the legend of the subplot of the similar term
-Red shows the average
- Function: plotmotif(termlist, hyperg_out, output, Boolean (True for overlay, False for subplot), cutoff_value, coefficient value (for clustering))


KNOWN PROBLEMS:
-Hypergeometric test can take a long time to run
	- Potential solution: Pre-compute and store results for all motifs
-Homology bias


NEXT STEPS:
- Speed up the programs
	-pre-compute data
-Clean up some of the codes
-Create a plot where colour changes based on p-value rank 
-Try to get it onto a webserver
