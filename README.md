# TaxSBP

Vitor C. Piro (vitorpiro@gmail.com)


Implementation of the approximation algorithm for the hierarchically structured bin packing problem [1] based on the NCBI Taxonomy database [2].

Dependencies:
-------------
 - binpacking=1.3 (https://pypi.python.org/pypi/binpacking)

Input: 
------
 * sequence information (identifier, length and taxonomic assignment) and taxonomy (nodes.dmp and merged.dmp)

sequence information (*_updated_sequence_accession.txt from https://github.com/pirovc/genome_updater):

    # Archaea - RefSeq - Complete Genomes ** also downloads taxonomy **
	./genome_updater.sh -d "refseq" -g "archaea" -c "all" -l "Complete Genome" -f "genomic.fna.gz,assembly_report.txt" -o refseq_archaea_cg -t 12 -r -a
	
taxonomy:
	
	# Unpack  
	tar xfz refseq_archaea_cg/${date}_taxdump.tar.gz nodes.dmp merged.dmp
	
Output:
-------
 * A tab-separated file with sequence identifier, length, taxonomic assignment and bin identifier

Running:
--------

	# Create new bins
	grep "^A" refseq_archaea_cg/${date}_updated_sequence_accession.txt | cut -f 3,4,5 > input.txt
	python3 taxsbp/TaxSBP.py create -f input.txt -n nodes.dmp -s 2157 -b 20 > bins.txt
	
	# Add new sequences to existing bins
	grep "^A" refseq_archaea_cg/${date}_updated_sequence_accession.txt | cut -f 3,4,5 > added.txt.
	python3 taxsbp/TaxSBP.py add -f added.txt -i bins.txt -n nodes.dmp -m merged.dmp > added_bins.txt
	
	# Remove sequences to existing bins
	grep "^R" refseq_archaea_cg/${date}_updated_sequence_accession.txt | cut -f 3 > removed.txt
	python3 taxsbp/TaxSBP.py remove -f removed.txt -i bins.txt > bins_updated.txt
	

References:
-----------

[1] Codenotti, B., De Marco, G., Leoncini, M., Montangero, M., & Santini, M. (2004). Approximation algorithms for a hierarchically structured bin packing problem. Information Processing Letters, 89(5), 215–221. http://doi.org/10.1016/j.ipl.2003.12.001

[2] Federhen, S. (2012). The NCBI Taxonomy database. Nucleic Acids Research, 40(D1), D136–D143. http://doi.org/10.1093/nar/gkr1178