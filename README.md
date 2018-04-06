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

	# Create bins with exclusive species (sequences of the same species can go to several bins but bins are always species exclusive)
	python3 TaxSBP.py create -f acc_len_txid.txt -n nodes.dmp -m merged.dmp -s 1 -b 20 -r species > bins.txt

	# Create bins with pre-clustered species (one bin can have sequences from several species, but sequences of the same species are always on the same bin)
	python3 TaxSBP.py create -f acc_len_txid.txt -n nodes.dmp -m merged.dmp -s 1 -b 20 -p species > bins.txt

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
	                     [-p <pre_cluster>] [-r <bin_exclusive>] [--use-group]

	optional arguments:
	  -h, --help          show this help message and exit
	  -f <input_file>     Tab-separated with the fields: sequence id, sequence
	                      length, taxonomic id [, group]
	  -n <nodes_file>     nodes.dmp from NCBI Taxonomy
	  -m <merged_file>    merged.dmp from NCBI Taxonomy
	  -s <start_node>     Start node taxonomic id. Default: 1 (root node)
	  -b <bins>           Approximate number of bins (estimated by total
	                      length/bin number). Default: 50 [Mutually exclusive -l]
	  -l <bin_len>        Maximum bin length (in bp). Use this parameter insted of
	                      -b to define the number of bins [Mutually exclusive -b]
	  -p <pre_cluster>    Pre-cluster sequences into ranks/taxid/group, so they
	                      won't be splitted among bins
	                      [none,group,taxid,species,genus,...] Default: none
	  -r <bin_exclusive>  Make bins rank/taxid/group exclusive, so bins won't have
	                      mixed sequences. When the chosen rank is not present on
	                      a sequence lineage, this sequence will be taxid
	                      exclusive. [none,group,taxid,species,genus,...] Default:
	                      none
	  --use-group         If activated, TaxSBP will further cluster sequences on a
	                      specialized level after the taxonomic id (e.g. assembly
	                      accession, strain name, etc). The group should be
	                      provided as an extra collumns in the input_file

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
