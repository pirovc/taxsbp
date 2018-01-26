#!/usr/bin/python3
# The MIT License (MIT)
# 
# Copyright (c) 2017 - Vitor C. Piro - PiroV@rki.de - vitorpiro@gmail.com
# Robert Koch-Institut, Germany
# All rights reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import binpacking
import argparse
import random
from collections import defaultdict

def main():
	parser = argparse.ArgumentParser(prog='TaxSBP',conflict_handler="resolve")
	subparsers = parser.add_subparsers()

	# create
	create_parser = subparsers.add_parser('create', help='Create new bins')
	create_parser.set_defaults(which='create')
	create_parser.add_argument('-f', required=True, metavar='<input_file>', dest="input_file", help="Tab-separated file with sequence id, sequence length and taxonomic id")
	create_parser.add_argument('-n', required=True, metavar='<nodes_file>', dest="nodes_file", help="nodes.dmp from NCBI Taxonomy")
	create_parser.add_argument('-m', required=False, metavar='<merged_file>', dest="merged_file", help="merged.dmp from NCBI Taxonomy")
	create_parser.add_argument('-s', default=1, metavar='<start_node>', dest="start_node", type=int, help="Start node taxonomic id. Default: 1 (root node)")
	create_parser.add_argument('-b', default=50, metavar='<bins>', dest="bins", type=int, help="Number of bins (estimated by sequence lenghts). Default: 50 [Mutually exclusive -l]")
	create_parser.add_argument('-l', metavar='<bin_len>', dest="bin_len", type=int, help="Maximum bin length. Use this parameter insted of -b to define the number of bins [Mutually exclusive -b]")
	
	# add
	add_parser = subparsers.add_parser('add', help='Add sequences to existing bins')
	add_parser.set_defaults(which='add')
	add_parser.add_argument('-f', required=True, metavar='<input_file>', dest="input_file", help="Tab-separated file with the NEW sequence ids, sequence length and taxonomic id")
	add_parser.add_argument('-i', required=True, metavar='<bins_file>', dest="bins_file", help="Previously generated bins")
	add_parser.add_argument('-n', required=True, metavar='<nodes_file>', dest="nodes_file", help="nodes.dmp from NCBI Taxonomy (new sequences)")
	add_parser.add_argument('-m', required=True, metavar='<merged_file>', dest="merged_file", help="merged.dmp from NCBI Taxonomy (new sequences)")
	add_parser.add_argument('--distribute', dest="distribute", default=False, action='store_true', help="Distribute newly added sequences among more bins. Without this option, TaxSBP will try to update as few bins as possible.")

	# remove
	remove_parser = subparsers.add_parser('remove', help='Remove sequences to existing bins')
	remove_parser.set_defaults(which='remove')
	remove_parser.add_argument('-f', required=True, metavar='<input_file>', dest="input_file", help="List of sequence ids to be removed")
	remove_parser.add_argument('-i', required=True, metavar='<bins_file>', dest="bins_file", help="Previously generated bins")
	
	# list
	list_parser = subparsers.add_parser('list', help='List bins based on sequence ids')
	list_parser.set_defaults(which='list')
	list_parser.add_argument('-f', required=True, metavar='<input_file>', dest="input_file", help="List of sequence ids")
	list_parser.add_argument('-i', required=True, metavar='<bins_file>', dest="bins_file", help="Previously generated bins")
	
	parser.add_argument('-v', action='version', version='%(prog)s 0.04')
	args = parser.parse_args()

	global parents
	global leaves
	global bin_len
		
	if args.which=="create":

		parents, leaves, accessions, total_len = read_input(args.input_file, args.start_node, read_nodes(args.nodes_file), read_merged(args.merged_file))
	
		if not parents[args.start_node]:
			print("No children nodes found / invalid taxid -", str(args.start_node))
			return
		
		# Bin length (estimated from number of bins or directly)
		if args.bin_len:
			bin_len = args.bin_len
		else:
			# Estimate bin len based on number of requested bins
			bin_len = total_len/float(args.bins) 
		
		# Run taxonomic structured bin packing
		final_bins = ApproxSBP(args.start_node) 			## RECURSIVE
		#final_bins = ApproxSBP_stack(args.start_node)		## STACK

		# Print resuls (by sequence)
		for binid,bin in enumerate(final_bins):
			for id in bin[1:]:
				# Output: accession, seq len, taxid, bin
				print(id, accessions[id][0], accessions[id][1], binid, sep="\t")

	elif args.which=="add":
		nodes = read_nodes(args.nodes_file)
		merged = read_merged(args.merged_file)
		bins, lens = read_bins(args.bins_file, nodes, merged)
		parents, leaves, accessions, total_len = read_input(args.input_file, 1, nodes, merged)
		# Sort accessions for reproducible output (choice on multiple bins when adding sequence lens)
		for accession,(length,taxid) in sorted(accessions.items()):
			t = taxid
			# If taxid of the entry is not directly assigned, look for LCA with assignment
			while not bins[t] and t!=1: t = nodes[t]

			# If taxid is assigned among several bins, chooses smallest to make the assignment
			if len(bins[t])>1:
				d = {b:lens[b] for b in bins[t]}
				bin = min(d, key=d.get)
			else:
				bin = list(bins[t])[0]
			
			# Add length count to the bin to be accounted in the next sequences
			# (makes it distribute more evenly, but splits similar sequences and affects more bins)
			if args.distribute: lens[bin]+=length 

			print(accession, length, taxid, bin, sep="\t")			

	elif args.which=="remove":
		r = set(line.split('\n')[0] for line in open(args.input_file,'r'))

		with open(args.bins_file,'r') as file:
			for line in file:
				accession = line.split('\t')[0]
				if accession not in r:
					print(line, end='')

	elif args.which=="list":
		r = set(line.split('\n')[0] for line in open(args.input_file,'r'))
		with open(args.bins_file,'r') as file:
			for line in file:
				accession = line.split('\t')[0]
				if accession in r:
					print(line, end='')


# Input: list of tuples [(seqlen, seqid1 [, ..., seqidN])]
# Output: bin packed list of tuples [(seqlen, seqid1 [, ..., seqidN])]
# Returns multi-valued tuple: first [summed] length summed followed by the id[s]
def bpck(d): 
	# Only one bin, no need to pack
	if len(d)==1: return d
	else:
		ret = []
		for bin in binpacking.to_constant_volume(d, bin_len, weight_pos=0):
			if bin: #Check if the returned bin is not empty: it happens when the bin packing algorith cannot divide larger sequences
				# Convert the bin listed output to tuple format
				sum_length = 0
				ids = []
				for i in bin:
					sum_length+=i[0]
					ids.extend(i[1:])
				ret.append((sum_length,*ids))
		return ret

def ApproxSBP(v):
	children = parents[v]
	
	# If it doesn't have any children it's a leaf and should return the packed sequences
	if not children: return bpck(leaves[v])
		
	# Recursively bin pack children
	# Sort children to keep it more consistent with different versions of the taxonomy (new taxids)
	ret = []
	for child in sorted(children): ret.extend(ApproxSBP(child))

	# if current node has sequences assigned to it (but it's not a leaf), add it to the current bin packing (it will first pack with its own children nodes)
	## QUESTION: should I bin together those sequeneces or "distribute" along its children -- command: ret.update(leaves[v]) -- 
	if leaves[v]: ret.extend(bpck(leaves[v]))

	return bpck(ret)

def read_nodes(nodes_file):
	# READ nodes -> fields (1:TAXID 2:PARENT_TAXID)
	nodes = {}
	with open(nodes_file,'r') as fnodes:
		for line in fnodes:
			taxid, parent_taxid, _ = line.split('\t|\t',2)
			nodes[int(taxid)] = int(parent_taxid)
	nodes[1] = 0 #Change parent taxid of the root node to 0 (it's usually 1 and causes infinite loop later)
	return nodes
	
def read_merged(merged_file):
	# READ nodes -> fields (1:OLD TAXID 2:NEW TAXID)
	merged = {}
	if merged_file:
		with open(merged_file,'r') as fmerged:
			for line in fmerged:
				old_taxid, new_taxid, _ = line.rstrip().split('\t|',2)
				merged[int(old_taxid)] = int(new_taxid)
	return merged
	
def read_input(input_file, start_node, nodes, merged):
	# READ input file -> fields (0:ACCESSION 1:LENGTH 2:TAXID)
	parents = defaultdict(set) # set cause it has faster lookup and it does not accept duplicated values (no need to check for that)
	leaves = defaultdict(list)
	accessions = dict()
	total_len = 0
	with open(input_file,'r') as file:
		for line in file:
			fields = line.split('\t')
			accession = fields[0]
			if accession=="na": continue # SKIP ENTRY WITH NO ACCESSION - TODO log
			length = int(fields[1])
			taxid = int(fields[2])
			if taxid not in nodes: 
				if taxid not in merged: 
					continue # SKIP ENTRY WITH NO TAXONOMIC ASSIGNEMNT - TODO log
				else:
					taxid = merged[taxid] # Get new taxid from merged.dmp
			leaves[taxid].append((length,accession)) # Keep length and accession for each taxid (multiple entries)
			accessions[accession] = (length,taxid) # Keep length and taxid for each accession (input file)
			while True: #Check all taxids in the lineage
				if taxid==start_node: total_len+=length # Just account sequence to total when it's on the sub-tree
				parents[nodes[taxid]].add(taxid) # Create parent:children structure only for used taxids
				if taxid==1: break
				taxid = nodes[taxid]
				
	return parents, leaves, accessions, total_len

def read_bins(bins_file, nodes, merged):
	# READ bins -> fields (0:ACCESSION 1:BIN)
	bins = defaultdict(set) # set cause it has faster lookup and it does not accept duplicated values (no need to check for that)
	lens = defaultdict(int)
	with open(bins_file,'r') as file:
		for line in file:
			fields = line.split('\t')
			accession = fields[0]
			length = int(fields[1])
			taxid = int(fields[2])
			bin = int(fields[3])
			lens[bin]+=length
			while True: #Check all taxids in the lineage
				bins[taxid].add(bin) # Create parent:children structure only for used taxids
				if taxid==1: break
				# If taxid is not present on newer version of nodes.dmp, look for merged entry
				try:
					taxid = nodes[taxid]
				except KeyError:
					taxid = merged[taxid] # TODO log if not found on both
					
	return bins, lens
	
#def ApproxSBP_stack(v):
#	stack = [v]
#	parent_stack = [1]
#	bins = defaultdict(list)
#	while stack:
#		taxid = stack[-1]
#		parent_taxid = parent_stack[-1]
#
#		# return to parent
#		if taxid==parent_taxid:
#			# if current node has sequences assigned to it (but it's not a leaf), add it to the current bin packing (it will first pack with its own children nodes)
#			if leaves[taxid]: bins[parent_taxid].extend(bpck(leaves[taxid]))
#			stack.pop()
#			parent_stack.pop()
#			# Add children packs to parent
#			bins[parent_stack[-1]].extend(bpck(bins[parent_taxid]))
#			del bins[parent_taxid]
#			continue
#			
#		# Check if node has children
#		children = parents[taxid]
#		
#		# If it doesn't have any children it's a leaf and should return the packed sequences
#		if not children: 
#			bins[parent_taxid].extend(bpck(leaves[taxid]))
#			stack.pop()
#		else:
#			parent_stack.append(taxid)
#			stack.extend(sorted(children,reverse=True)) # Sort children (reversed because it's in a stack) to keep it more consistent with different versions of the taxonomy (new taxids)
#
#	return bins[1]
	
if __name__ == "__main__":
	main()
