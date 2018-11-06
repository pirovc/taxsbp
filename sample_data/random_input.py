import random
from collections import defaultdict
import sys

def parse_nodes_ranks(nodes_file):
	# READ nodes -> fields (1:TAXID 2:PARENT_TAXID 3:RANK)
	nodes = {}
	ranks = {}
	children = defaultdict(set)
	with open(nodes_file,'r') as fnodes:
		for line in fnodes:
			taxid, parent_taxid, rank, _ = line.split('\t|\t',3)
			ranks[taxid] = rank
			nodes[taxid] = parent_taxid
			children[parent_taxid].add(taxid)
	nodes['1'] = '0' #Change parent taxid of the root node to 0 (it's usually 1 and causes infinite loop later)
	
	return nodes, children, ranks

nodes, children, ranks = parse_nodes_ranks(sys.argv[1])

n_of_sequences = 10
min_len = 100
max_len = 2000
specialization = True
selected_ranks = ['b','c', 'd', 'e', 'f']


# create list of taxids (for all nodes or only selected ranks)
if selected_ranks:
	taxids = [n for n in nodes.keys() if ranks[n] in selected_ranks]
else:
	taxids = list(nodes.keys())

# seqid:taxid structure, selecting random taxids for the each sequence
seqid_taxid = {i:taxids[random.randrange(len(taxids))] for i in range(n_of_sequences)}

# taxid:set(seqids) reverse structure to count number of entries of each taxid
taxid_seqids = defaultdict(set)
for s,t in seqid_taxid.items():
	taxid_seqids[t].add(s)

# print
for seqid in range(n_of_sequences):
	seq = "seq_" + str(seqid)
	seqlen = random.randint(min_len, max_len)
	taxid = seqid_taxid[seqid]
	
	if not specialization:
		print(seq, seqlen, taxid, sep="\t")
	else:
		spec = "spec_" + str(random.sample(taxid_seqids[taxid],1)[0]) 
		print(seq, seqlen, taxid, spec, sep="\t")

