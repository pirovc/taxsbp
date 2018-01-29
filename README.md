# TaxSBP

Vitor C. Piro (vitorpiro@gmail.com)

Implementation of the approximation algorithm for the hierarchically structured bin packing problem [1] based on the NCBI Taxonomy database [2].

Dependencies:
-------------
binpacking=1.3 (https://pypi.python.org/pypi/binpacking)
 
	pip install binpacking

Input: 
------
 * A tab-separated file with sequence id, sequence length and taxonomic id 
 * nodes.dmp and merged.dmp from NCBI Taxonomy (ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
	
Output:
-------
 * A tab-separated file with sequence id, sequence length, taxonomic id and bin number

Example:
--------

	# From one or more FASTA files, grep accession version identifier
	grep -o "^>\S*" sequences.fna | tr -d ">" > accessions.txt

	# Retrieve length and taxonomic id from NCBI eutils (24 threads)
	<accessions.txt xargs --max-procs=24 -I '{}' bash get_len_taxid.sh '{}' > acc_len_txid.txt

	# Get taxonomy
	wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
	tar xfz taxdump.tar.gz nodes.dmp merged.dmp

	# Run TaxSBP
	python3 TaxSBP.py create -f acc_len_txid.txt -n nodes.dmp -m merged.dmp -s 1 -b 20 > bins.txt
	
	# Re-write fasta into separated files 
	./split_fna_bins.sh sequences.fna bins.txt output_folder/

Further options:
----------------

	# Add new sequences to existing bins (updating as few bins as possible)
	python3 TaxSBP.py add -f ADDED_acc_len_txid.txt -i bins.txt -n nodes.dmp -m merged.dmp > added_bins.txt
	
	# Add new sequences to existing bins distributin (distributing the sequences as much as possible)
	python3 TaxSBP.py add --distribute -f ADDED_acc_len_txid.txt -i bins.txt -n nodes.dmp -m merged.dmp > added_bins.txt
	
	# Remove sequences from existing bins
	python3 TaxSBP.py remove -f accessions_to_remove.txt -i bins.txt > all_bins_updated.txt

	# List sequences from existing bins
	python3 TaxSBP.py list -f accessions_to_list.txt -i bins.txt > bins_from_list.txt
	
Parameters:
-----------

$ python3 TaxSBP.py -h

	usage: TaxSBP [-h] [-v] {create,add,remove,list} ...

	positional arguments:
	  {create,add,remove,list}
	    create              Create new bins
	    add                 Add sequences to existing bins
	    remove              Remove sequences to existing bins
	    list                List bins based on sequence ids

	optional arguments:
	  -h, --help            show this help message and exit
	  -v                    show program's version number and exit

  
$ python3 TaxSBP.py create -h

	usage: TaxSBP create [-h] -f <input_file> -n <nodes_file> [-m <merged_file>]
	                     [-s <start_node>] [-b <bins>] [-l <bin_len>]

	optional arguments:
	  -h, --help        show this help message and exit
	  -f <input_file>   Tab-separated file with sequence id, sequence length and
	                    taxonomic id
	  -n <nodes_file>   nodes.dmp from NCBI Taxonomy
	  -m <merged_file>  merged.dmp from NCBI Taxonomy
	  -s <start_node>   Start node taxonomic id. Default: 1 (root node)
	  -b <bins>         Number of bins (estimated by sequence lenghts). Default:
	                    50 [Mutually exclusive -l]
	  -l <bin_len>      Maximum bin length. Use this parameter insted of -b to
	                    define the number of bins [Mutually exclusive -b]

$ python3 TaxSBP.py add -h

	usage: TaxSBP add [-h] -f <input_file> -i <bins_file> -n <nodes_file> -m
	                  <merged_file> [--distribute]

	optional arguments:
	  -h, --help        show this help message and exit
	  -f <input_file>   Tab-separated file with the NEW sequence ids, sequence
	                    length and taxonomic id
	  -i <bins_file>    Previously generated bins
	  -n <nodes_file>   nodes.dmp from NCBI Taxonomy (new sequences)
	  -m <merged_file>  merged.dmp from NCBI Taxonomy (new sequences)
	  --distribute      Distribute newly added sequences among more bins. Without
	                    this option, TaxSBP will try to update as few bins as
	                    possible.

$ python3 TaxSBP.py remove -h

	usage: TaxSBP remove [-h] -f <input_file> -i <bins_file>

	optional arguments:
	  -h, --help       show this help message and exit
	  -f <input_file>  List of sequence ids to be removed
	  -i <bins_file>   Previously generated bins	   

$ python3 TaxSBP.py list -h

	usage: TaxSBP list [-h] -f <input_file> -i <bins_file>

	optional arguments:
	  -h, --help       show this help message and exit
	  -f <input_file>  List of sequence ids
	  -i <bins_file>   Previously generated bins

References:
-----------

[1] Codenotti, B., De Marco, G., Leoncini, M., Montangero, M., & Santini, M. (2004). Approximation algorithms for a hierarchically structured bin packing problem. Information Processing Letters, 89(5), 215–221. http://doi.org/10.1016/j.ipl.2003.12.001

[2] Federhen, S. (2012). The NCBI Taxonomy database. Nucleic Acids Research, 40(D1), D136–D143. http://doi.org/10.1093/nar/gkr1178
