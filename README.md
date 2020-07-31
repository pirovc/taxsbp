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

## Examples with sample data:

The sample data is comprised of:

	- Root node (1)
	- hierarchy with 5 levels (rank1-rank5)
	- Leaf nodes (*)
	- 8 Specializations (S1-S9)
	- 13 targets (A-M) with equal length (100) 

	rank-1                   1 ________________
	                        / \           \    \
	rank-2                2.1 2.2 ______   \    \
	                      / \    \      \   \    \
	rank-3             3.1  3.2   3.4    \   \    \
	                    /   / \     \     \   \    \
	rank-4          *4.1 *4.2 *4.3   4.4  *4.5 *4.6 \
	                 /    /   /     / | \   \   \    \
	rank-5          /    /   / *5.1 *5.2 \   \   \    \
	               /    /   /     /   |   \   \   \    \
	spec.         S1   S2  S3   S4   S5   S6  S7  S8   S9
	             /  \   |   |   / \   |   |   |   /|\   |
	target      A    B  C   D  E   F  G   H   I  J K L  M

### Clusters of size 400

	taxsbp.py -i sample_data/seqinfo.tsv -n sample_data/nodes.dmp -l 400

Clusters are limited by size 400

<details>
  <summary>Results</summary>

	#id	st	end	len	tax	bin	
	F	1	100	100	5.1	0
	E	1	100	100	5.1	0
	H	1	100	100	4.4	0
	G	1	100	100	5.2	0
	D	1	100	100	4.3	1
	C	1	100	100	4.2	1
	B	1	100	100	4.1	1
	A	1	100	100	4.1	1
	L	1	100	100	4.6	2
	K	1	100	100	4.6	2
	J	1	100	100	4.6	2
	M	1	100	100	1	2
	I	1	100	100	4.5	3

</details>

### Pre-clustered by "rank-2"

	taxsbp.py -i sample_data/seqinfo.tsv -n sample_data/nodes.dmp -l 400 -p "rank-2"

ABCD (2.1) and EFGHI (2.2) are forced together (even if bigger than 400)

<details>
  <summary>Results</summary>

	#id	st	end	len	tax	bin	
	I	1	100	100	4.5	0
	G	1	100	100	5.2	0
	E	1	100	100	5.1	0
	F	1	100	100	5.1	0
	H	1	100	100	4.4	0
	C	1	100	100	4.2	1
	D	1	100	100	4.3	1
	A	1	100	100	4.1	1
	B	1	100	100	4.1	1
	J	1	100	100	4.6	2
	K	1	100	100	4.6	2
	L	1	100	100	4.6	2
	M	1	100	100	1	2

</details>

### Bin exclusive clusters by "rank-4"

	taxsbp.py -i sample_data/seqinfo.tsv -n sample_data/nodes.dmp -l 400 -e "rank-4"

Clusters are generated for each sub-tree of nodes from rank-4. The used rank is printed instead of the original.

<details>
  <summary>Results</summary>

	#id	st	end	len	tax	bin	
	F	1	100	100	4.4	0
	E	1	100	100	4.4	0
	H	1	100	100	4.4	0
	G	1	100	100	4.4	0
	L	1	100	100	4.6	1
	K	1	100	100	4.6	1
	J	1	100	100	4.6	1
	B	1	100	100	4.1	2
	A	1	100	100	4.1	2
	I	1	100	100	4.5	3
	D	1	100	100	4.3	4
	C	1	100	100	4.2	5
	M	1	100	100	1	6

</details>

### Cluster with specialization

	taxsbp.py -i sample_data/seqinfo.tsv -n sample_data/nodes.dmp -l 400 -s MySpec -e MySpec

Clusters are exclusive by specialization

<details>
  <summary>Results</summary>

	#id	st	end	len	tax	bin	spec	
	L	1	100	100	4.6	0	S8
	K	1	100	100	4.6	0	S8
	J	1	100	100	4.6	0	S8
	F	1	100	100	5.1	1	S4
	E	1	100	100	5.1	1	S4
	B	1	100	100	4.1	2	S1
	A	1	100	100	4.1	2	S1
	G	1	100	100	5.2	3	S5
	D	1	100	100	4.3	4	S3
	M	1	100	100	1	5	S9
	H	1	100	100	4.4	6	S6
	I	1	100	100	4.5	7	S7
	C	1	100	100	4.2	8	S2

</details>

### Cluster with fragmentation

	taxsbp.py -i sample_data/seqinfo.tsv -n sample_data/nodes.dmp -l 150 -f 50 -a 5 -e "rank-3"

Clusters of size 150. Fragment inputs in 50 with overlap of 5. Cluster exclusive of "rank-3"

<details>
  <summary>Results</summary>

	#id	st	end	len	tax	bin	
	F	1	55	55	3.4	0
	E	1	55	55	3.4	0
	B	1	55	55	3.1	1
	A	1	55	55	3.1	1
	L	1	55	55	4.6	2
	K	1	55	55	4.6	2
	G	1	55	55	3.4	3
	G	51	100	50	3.4	3
	H	1	55	55	3.4	4
	H	51	100	50	3.4	4
	D	1	55	55	3.2	5
	D	51	100	50	3.2	5
	C	1	55	55	3.2	6
	C	51	100	50	3.2	6
	I	1	55	55	4.5	7
	I	51	100	50	4.5	7
	J	1	55	55	4.6	8
	L	51	100	50	4.6	8
	M	1	55	55	1	9
	M	51	100	50	1	9
	F	51	100	50	3.4	10
	E	51	100	50	3.4	10
	B	51	100	50	3.1	11
	A	51	100	50	3.1	11
	K	51	100	50	4.6	12
	J	51	100	50	4.6	12

</details>

## Examples with real data:

### Prepare data:

	gzip -d sample_data/*.gz

### Clustering:

	taxsbp.py -i sample_data/20181219_abfv_refseq_cg.tsv -n sample_data/20181219_abfv_refseq_cg_nodes.dmp -l 10000000 -f 999500 -a 500

## Parameters:

$ taxsbp -h

	usage: TaxSBP [-h] [-i <input_file>] [-o <output_file>] [-n <nodes_file>] [-m <merged_file>] [-b <bins>] [-l <bin_len>]
	              [-f <fragment_len>] [-a <overlap_len>] [-p <pre_cluster>] [-e <bin_exclusive>] [-s <specialization>]
	              [-u <update_file>] [--output-unique-seqid] [-v]

	optional arguments:
	  -h, --help            show this help message and exit
	  -i <input_file>, --input-file <input_file>
	                        Tab-separated with the fields: sequence id <tab> sequence length <tab> taxonomic id [<tab>
	                        specialization]
	  -o <output_file>, --output-file <output_file>
	                        Path to the output tab-separated file with the fields. Default: STDOUT
	  -n <nodes_file>, --nodes-file <nodes_file>
	                        nodes.dmp from NCBI Taxonomy
	  -m <merged_file>, --merged-file <merged_file>
	                        merged.dmp from NCBI Taxonomy
	  -b <bins>, --bins <bins>
	                        Approximate number of bins (estimated by total length/bin number). [Mutually exclusive -l]
	  -l <bin_len>, --bin-len <bin_len>
	                        Maximum bin length (in bp). Use this parameter insted of -b to define the number of bins.
	                        Default: length of the biggest group [Mutually exclusive -b]
	  -f <fragment_len>, --fragment-len <fragment_len>
	                        Fragment sequences into pieces, output accession will be modified with positions:
	                        ACCESION/start:end
	  -a <overlap_len>, --overlap-len <overlap_len>
	                        Overlap length between fragments [Only valid with -a]
	  -p <pre_cluster>, --pre-cluster <pre_cluster>
	                        Pre-cluster sequences into rank/taxid/specialization, so they won't be splitted among bins
	                        [none,specialization name,taxid,species,genus,...] Default: none
	  -e <bin_exclusive>, --bin-exclusive <bin_exclusive>
	                        Make bins rank/taxid/specialization exclusive, so bins won't have mixed sequences. When the
	                        chosen rank is not present on a sequence lineage, this sequence will be taxid/specialization
	                        exclusive. [none,specialization name,taxid,species,genus,...] Default: none
	  -s <specialization>, --specialization <specialization>
	                        Specialization name (e.g. assembly, strain). If given, TaxSBP will cluster entries on a
	                        specialized level after the taxonomic id. The specialization identifier should be provided as an
	                        extra collumn in the input_file ans should respect the taxonomic hiercharchy (one taxid ->
	                        multiple specializations / one specialization -> one taxid). Default: ''
	  -u <update_file>, --update-file <update_file>
	                        Previously generated files to be updated. Default: ''
	  --output-unique-seqid
	                        Output unique sequence ids after fragmentation in the format: seq.id/seq.start:seq.end]
	  -v, --version         show program's version number and exit

## Manual Installation

### Dependencies:

- python>=3.4
- [binpacking](https://pypi.org/project/binpacking/)==1.4.3
- [pylca](https://github.com/pirovc/pylca)==1.0.0
- [pandas](https://pypi.org/project/pandas/)pandas>=0.22.0 (tests only)

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

### Testing:

```shh
pip install "pandas>=0.22.0"
cd taxsbp
python3 -m unittest discover -s tests/taxsbp/integration/
```

References:
-----------

[1] Codenotti, B., De Marco, G., Leoncini, M., Montangero, M., & Santini, M. (2004). Approximation algorithms for a hierarchically structured bin packing problem. Information Processing Letters, 89(5), 215–221. http://doi.org/10.1016/j.ipl.2003.12.001

[2] Federhen, S. (2012). The NCBI Taxonomy database. Nucleic Acids Research, 40(D1), D136–D143. http://doi.org/10.1093/nar/gkr1178

[3] https://www.ics.uci.edu/~eppstein/ in the package https://github.com/pirovc/pylca
