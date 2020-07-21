# TaxSBP

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/taxsbp/README.html)

[![Build Status](https://travis-ci.org/pirovc/taxsbp.svg?branch=master)](https://travis-ci.org/pirovc/taxsbp) 

Implementation of the approximation algorithm for the hierarchically structured bin packing problem [1] based on the NCBI Taxonomy database [2] (uses LCA script from [3]).

## Installation

```shh
conda -c bioconda -c conda-forge taxsbp
```
or [manual installation](#manual-installation) without conda

## Usage

### Input

 * A tab-separated file:
	
	`sequence id <tab> sequence length <tab> taxonomic id [ <tab> specialization]`
 
 * nodes.dmp and merged.dmp from NCBI Taxonomy (ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)

### Output 

 * A tab-separated file:

 	`sequence id <tab> seq.start <tab> seq.end <tab> sequence length <tab> taxonomic id <tab> bin id [ <tab> specialization] `

## Running sample data:

13 Sequences (A-M), 8 Specializations (S1-S9) distributed in the following hierarchy with 5 levels (rank1-5). Root node = 1. Leaf nodes marked with *

    A    100   4.1    S1
    B    50    4.1    S1
    C    300   4.2    S2
    D    100   4.3    S3
    E    40    5.1    S4
    F    50    5.1    S4
    G    90    5.2    S5
    H    1000  5.3    S6
    I    10    4.5    S7
    J    17    4.6    S8
    K    5     4.6    S8
    L    300   4.6    S8
    M    733   1      S9

    rank-1                   1 ________________
                            / \           \    \
    rank-2                2.1 2.2 ______   \    \
                          / \    \      \   \    \
    rank-3             3.1  3.2   3.4    \   \    \
                        /   / \     \     \   \    \
    rank-4          *4.1 *4.2 *4.3   4.4  *4.5 *4.6 \
                     /    /   /     / | \   \   \    \
    rank-5          /    /   / *5.1*5.2*5.3  \   \    \
                   /    /   /    |   |   |    \   \    \
    spec.         S1   S2  S3   S4  S5   S6   S7  S8    S9
                 /  \   |   |   / \  |   |    |   /|\    |
    sequence    A    B  C   D  E  F  G   H    I  J K L   M
    length     100  50 300 100 40 50 90 1000 10 17 5 300 733

	taxsbp.py -i sample_data/seqinfo.tsv -n sample_data/nodes.dmp

## Parameters:

$ taxsbp -h

	usage: TaxSBP [-h] -f <input_file> -n <nodes_file> [-m <merged_file>]
              [-b <bins>] [-l <bin_len>] [-a <fragment_len>]
              [-o <overlap_len>] [-p <pre_cluster>] [-r <bin_exclusive>]
              [-z <specialization>] [-u <update_file>] [-v]

	optional arguments:
	  -h, --help           show this help message and exit
	  -f <input_file>      Tab-separated with the fields: sequence id <tab>
	                       sequence length <tab> taxonomic id [<tab>
	                       specialization]
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
	                       taxid/specialization exclusive. [none,specialization
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


## Manual Installation

### Dependencies:

- python>=3.4
- [binpacking](https://pypi.org/project/binpacking/)==1.4.1
- [pylca](https://github.com/pirovc/pylca)==1.0.0

### Pylca:

```shh
git clone https://github.com/pirovc/pylca
cd pylca
python setup.py install
```

### TaxSBP + binpacking:

```shh
git clone https://github.com/pirovc/taxsbp.git
cd taxsbp
python setup.py install
taxsbp -h
```


References:
-----------

[1] Codenotti, B., De Marco, G., Leoncini, M., Montangero, M., & Santini, M. (2004). Approximation algorithms for a hierarchically structured bin packing problem. Information Processing Letters, 89(5), 215–221. http://doi.org/10.1016/j.ipl.2003.12.001

[2] Federhen, S. (2012). The NCBI Taxonomy database. Nucleic Acids Research, 40(D1), D136–D143. http://doi.org/10.1093/nar/gkr1178

[3] https://www.ics.uci.edu/~eppstein/ in the package https://github.com/pirovc/pylca
