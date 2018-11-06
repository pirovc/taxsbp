from pptree import *
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


def generate_tree(n):
	
	if n.name in children:
		for child in children[n.name]: 
			generate_tree(Node(child, n))
	
	return n

nodes, children, ranks = parse_nodes_ranks(sys.argv[1])

print_tree(generate_tree(Node('1')))
