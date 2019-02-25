# TaxSBP

Vitor C. Piro (vitorpiro@gmail.com)

Implementation of the approximation algorithm for the hierarchically structured bin packing problem [1] based on the NCBI Taxonomy database [2].

Dependencies:
-------------

python >=3.5.0

binpacking=1.3 (https://pypi.python.org/pypi/binpacking)
 
	pip install binpacking==1.3

or install manually into TaxSBP.py folder:
	
	wget https://pypi.python.org/packages/c9/fe/56782753922a195d332d419949f889c1d59cab7b1780db2351bd8b99501c/binpacking-1.3.tar.gz
	tar -xvf binpacking-1.3.tar.gz binpacking-1.3/binpacking/ --strip-components=1
	
Input: 
------
 * A tab-separated file:
	
	`sequence id <tab> sequence length <tab> taxonomic id [ <tab> specialization]`
 
 * nodes.dmp and merged.dmp from NCBI Taxonomy (ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
	
Output:
-------
 * A tab-separated file:

 	`sequence id <tab> sequence length <tab> taxonomic id [ <tab> specialization] <tab> bin id`

Example:
--------

	# From one or more FASTA files, grep accession version identifier
	grep -o "^>\S*" sequences.fna | tr -d ">" > accessions.txt

	# Retrieve length and taxonomic id from NCBI eutils
	scripts/get_len_taxid.sh accessions.txt > acc_len_txid.txt

	# Get taxonomy
	wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
	tar xfz taxdump.tar.gz nodes.dmp merged.dmp

	# Run TaxSBP (clusters of 15M bases)
	python3 TaxSBP.py -f acc_len_txid.txt -n nodes.dmp -m merged.dmp -l 15000000 > bins.txt
	
	# Re-write fasta into separated files 
	scripts/split_fna_bins.sh sequences.fna bins.txt output_folder/

Parameters:
-----------

$ python3 TaxSBP.py -h

	usage: TaxSBP [-h] -f <input_file> -n <nodes_file> [-m <merged_file>]
	              [-b <bins>] [-l <bin_len>] [-a <fragment_len>]
	              [-o <overlap_len>] [-p <pre_cluster>] [-r <bin_exclusive>]
	              [-z <specialization>] [-u <update_file>] [-v]

	optional arguments:
	  -h, --help           show this help message and exit
	  -f <input_file>      Tab-separated with the fields: sequence id <tab>
	                       sequence length <tab> taxonomic id [<tab> group]
	  -n <nodes_file>      nodes.dmp from NCBI Taxonomy
	  -m <merged_file>     merged.dmp from NCBI Taxonomy
	  -b <bins>            Approximate number of bins (estimated by total
	                       length/bin number). Default: 50 [Mutually exclusive -l]
	  -l <bin_len>         Maximum bin length (in bp). Use this parameter insted
	                       of -b to define the number of bins [Mutually exclusive
	                       -b]
	  -a <fragment_len>    Fragment sequences into pieces, output accession will
	                       be modified with positions: ACCESION/start:end
	  -o <overlap_len>     Overlap length between fragments [Only valid with -a]
	  -p <pre_cluster>     Pre-cluster sequences into rank/taxid/specialization,
	                       so they won't be splitted among bins
	                       [none,specialization name,taxid,species,genus,...]
	                       Default: none
	  -r <bin_exclusive>   Make bins rank/taxid/specialization exclusive, so bins
	                       won't have mixed sequences. When the chosen rank is not
	                       present on a sequence lineage, this sequence will be
	                       taxid/group exclusive. [none,specialization
	                       name,taxid,species,genus,...] Default: none
	  -z <specialization>  Specialization name (e.g. assembly, strain). If given,
	                       TaxSBP will cluster entries on a specialized level
	                       after the taxonomic id. The specialization identifier
	                       should be provided as an extra collumn in the
	                       input_file ans should respect the taxonomic hiercharchy
	                       (one taxid -> multiple specializations / one
	                       specialization -> one taxid). Default: ''
	  -u <update_file>     Previously generated files to be updated. Default: ''
	  -v                   show program's version number and exit


References:
-----------

[1] Codenotti, B., De Marco, G., Leoncini, M., Montangero, M., & Santini, M. (2004). Approximation algorithms for a hierarchically structured bin packing problem. Information Processing Letters, 89(5), 215–221. http://doi.org/10.1016/j.ipl.2003.12.001

[2] Federhen, S. (2012). The NCBI Taxonomy database. Nucleic Acids Research, 40(D1), D136–D143. http://doi.org/10.1093/nar/gkr1178
