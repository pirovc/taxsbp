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
import sys
from collections import defaultdict

def main():
	parser = argparse.ArgumentParser(prog='TaxSBP',conflict_handler="resolve")
	subparsers = parser.add_subparsers()

	# create
	create_parser = subparsers.add_parser('create', help='Create new bins')
	create_parser.set_defaults(which='create')
	create_parser.add_argument('-f', required=True, metavar='<input_file>', dest="input_file", help="Tab-separated with the fields: sequence id, sequence length, taxonomic id [, group]")
	create_parser.add_argument('-n', required=True, metavar='<nodes_file>', dest="nodes_file", help="nodes.dmp from NCBI Taxonomy")
	create_parser.add_argument('-m', required=False, metavar='<merged_file>', dest="merged_file", help="merged.dmp from NCBI Taxonomy")
	create_parser.add_argument('-s', default=1, metavar='<start_node>', dest="start_node", type=int, help="Start node taxonomic id. Default: 1 (root node)")
	create_parser.add_argument('-b', default=50, metavar='<bins>', dest="bins", type=int, help="Approximate number of bins (estimated by total length/bin number). Default: 50 [Mutually exclusive -l]")
	create_parser.add_argument('-l', metavar='<bin_len>', dest="bin_len", type=int, help="Maximum bin length. Use this parameter insted of -b to define the number of bins [Mutually exclusive -b]")
	create_parser.add_argument('-p', metavar='<pre_cluster>', dest="pre_cluster", type=str, default="none", help="Pre-cluster sequences into ranks/taxids, so they won't be splitted among bins [none,group,taxid,species,genus,...] Default: none")
	create_parser.add_argument('-r', metavar='<rank_exclusive>', dest="rank_exclusive", type=str, default="none", help="Make bins rank/taxid exclusive, so bins won't have mixed ranks. When the chosen rank is not available for a organism, this option will make this organism taxid exclusive . [none,group,taxid,species,genus,...] Default: none")
	create_parser.add_argument('--use-group', dest="use_group", default=False, action='store_true', help="If activated, TaxSBP will further classify sequences on a specialized level after the taxonomic id by the group (e.g. assembly accession, strain name, etc). Group should be provided in the input_file")
	
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
	
	parser.add_argument('-v', action='version', version='%(prog)s 0.06')
	args = parser.parse_args()

	if len(sys.argv[1:])==0: # Print help calling script without parameters
		parser.print_help() 
		return 0

	if args.which=="create":
		nodes, ranks = read_nodes(args.nodes_file, True if args.pre_cluster not in ["none", "taxid"] or args.rank_exclusive not in ["none", "taxid"] else False)
		leaves, accessions, total_len, nodes, ranks = read_input(args.input_file, args.start_node, nodes, read_merged(args.merged_file), ranks, args.use_group)
		parents = parents_tree(leaves, nodes)

		if not parents[args.start_node]:
			print("No children nodes found / invalid taxid - ", str(args.start_node))
			return
		
		# Bin length (estimated from number of bins or directly)
		if args.bin_len:
			bin_len = args.bin_len
		else:
			# Estimate bin len based on number of requested bins
			bin_len = total_len/float(args.bins) 

		# Check for valid rank
		if args.pre_cluster != "none":
			possible_ranks = set(ranks.values())
			if args.pre_cluster!="taxid" and args.pre_cluster not in possible_ranks:
				print("Rank not found: " + args.pre_cluster)
				print("Possible ranks: " + ', '.join(possible_ranks))
				return
		if args.rank_exclusive != "none":
			possible_ranks = set(ranks.values())
			if args.rank_exclusive!="taxid" and args.rank_exclusive not in possible_ranks:
				print("Rank not found: " + args.rank_exclusive)
				print("Possible ranks: " + ', '.join(possible_ranks))
				return

		# Pre-clustering
		if args.pre_cluster != "none":
			leaves = pre_cluster(args.pre_cluster, leaves, nodes, ranks)
			parents = parents_tree(leaves, nodes) #Re-calculate parents dict

		if args.rank_exclusive == "none":
			# Run taxonomic structured bin packing for the whole tree
			final_bins = ApproxSBP(args.start_node, leaves, parents, bin_len) ## RECURSIVE
		else:
			# Run taxonomic structured bin packing for each rank/taxid separetly
			final_bins = []
			unique_rank_taxids, single_taxids = get_unique_rank_taxids(args.rank_exclusive, leaves, nodes, ranks, args.start_node)
			
			if unique_rank_taxids:
				# clustering directly on the rank chosen, recursion required for children nodes
				for unique_rank_taxid in unique_rank_taxids:
					final_bins.extend(ApproxSBP(unique_rank_taxid, leaves, parents, bin_len))

			if single_taxids:
				# clustering directly on the taxid level, no recursion to children nodes necessary
				for single_taxid in single_taxids:
					final_bins.extend(bpck(leaves[single_taxid], bin_len))
			
		# Print resuls (by sequence)
		for binid,bin in enumerate(final_bins):
			for id in bin[1:]:
				if args.use_group:
					print(id, accessions[id][0], nodes[accessions[id][1]], accessions[id][1], binid, sep="\t")
				else:
					# Output: accession, seq len, taxid, bin
					print(id, accessions[id][0], accessions[id][1], binid, sep="\t")

	elif args.which=="add":
		nodes, _ = read_nodes(args.nodes_file, False)
		merged = read_merged(args.merged_file)
		bins, lens = read_bins(args.bins_file, nodes, merged)
		_ , accessions , _ , _ , _ = read_input(args.input_file, 1, nodes, merged, {}, args.use_group)
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

def sum_tuple_ids(bin):
	sum_length = 0
	ids = []
	for i in bin:
		sum_length+=i[0]
		ids.extend(i[1:])
	return sum_length, ids

# Input: list of tuples [(seqlen, seqid1 [, ..., seqidN])]
# Output: bin packed list of tuples [(seqlen, seqid1 [, ..., seqidN])]
# Returns multi-valued tuple: first [summed] length summed followed by the id[s]
def bpck(d, bin_len): 
	# Only one bin, no need to pack
	if len(d)==1: return d
	else:
		ret = []
		for bin in binpacking.to_constant_volume(d, bin_len, weight_pos=0):
			if bin: #Check if the returned bin is not empty: it happens when the bin packing algorith cannot divide larger sequences
				# Convert the bin listed output to tuple format
				sum_length, ids = sum_tuple_ids(bin)
				ret.append((sum_length,*ids))

		return ret

def ApproxSBP(v, leaves, parents, bin_len):
	children = parents[v]
	
	# If it doesn't have any children it's a leaf and should return the packed sequences
	if not children: return bpck(leaves[v], bin_len)
		
	# Recursively bin pack children
	# Sort children to keep it more consistent with different versions of the taxonomy (new taxids), str to for groups
	ret = []
	for child in sorted(children, key=str): ret.extend(ApproxSBP(child, leaves, parents, bin_len))

	# if current node has sequences assigned to it (but it's not a leaf), add it to the current bin packing (it will first pack with its own children nodes)
	## QUESTION: should I bin together those sequeneces or "distribute" along its children -- command: ret.update(leaves[v]) -- 
	if leaves[v]: ret.extend(bpck(leaves[v], bin_len))

	return bpck(ret, bin_len)

def read_nodes(nodes_file, pre_cluster):
	# READ nodes -> fields (1:TAXID 2:PARENT_TAXID)
	nodes = {}
	ranks = {} # Only collect ranks in case pre_cluster is required
	with open(nodes_file,'r') as fnodes:
		for line in fnodes:
			taxid, parent_taxid, rank, _ = line.split('\t|\t',3)
			if pre_cluster: ranks[int(taxid)] = rank
			nodes[int(taxid)] = int(parent_taxid)
	nodes[1] = 0 #Change parent taxid of the root node to 0 (it's usually 1 and causes infinite loop later)
	return nodes, ranks
	
def read_merged(merged_file):
	# READ nodes -> fields (1:OLD TAXID 2:NEW TAXID)
	merged = {}
	if merged_file:
		with open(merged_file,'r') as fmerged:
			for line in fmerged:
				old_taxid, new_taxid, _ = line.rstrip().split('\t|',2)
				merged[int(old_taxid)] = int(new_taxid)
	return merged
	
def read_input(input_file, start_node, nodes, merged, ranks, use_group):
	# READ input file -> fields (0:ACCESSION 1:LENGTH 2:TAXID [3:GROUP])
	leaves = defaultdict(list)
	accessions = dict()
	total_len = 0

	with open(input_file,'r') as file:
		for line in file:
			fields = line.rstrip().split('\t')
			accession = fields[0]
			if accession=="na": continue # SKIP ENTRY WITH NO ACCESSION - TODO log
			length = int(fields[1])
			taxid = int(fields[2])
		
			if taxid not in nodes: 
				if taxid not in merged: 
					continue # SKIP ENTRY WITH NO TAXONOMIC ASSIGNEMNT - TODO log
				else:
					taxid = merged[taxid] # Get new taxid from merged.dmp
			
			# Check if taxid is on the requested tree (under start node). 1=root, always on the tree
			# do it just one for each taxid (once it was found on leaves means it was already checked)
			if start_node==1 or taxid in leaves:
				ontree=True
			else:
				ontree = False
				lin_taxid = taxid
				while lin_taxid!=1: #Check all taxids in the lineage
					if lin_taxid==start_node: # If start node is found, is on the sub-tree
						ontree=True
						break
					lin_taxid = nodes[lin_taxid]
			
			# If yes, add to dicts
			if ontree:
				total_len+=length # Account for seq. length
				if use_group:
					group = fields[3]
					nodes[group] = taxid # Update nodes with new leaf
					ranks[group] = "group"
					leaves[group].append((length,accession)) # add group as leaf
					accessions[accession] = (length,group) 
				else:
					leaves[taxid].append((length,accession)) # Keep length and accession for each taxid (multiple entries)
					accessions[accession] = (length,taxid) # Keep length and taxid for each accession (input file)
			else:
				pass # TODO log
				
	return leaves, accessions, total_len, nodes, ranks

def parents_tree(leaves, nodes):
	# Define parent tree for faster lookup (set unique entries)
	parents = defaultdict(set) 
	for leaf_taxid in leaves.keys():
		while True: #Check all taxids in the lineage
			parents[nodes[leaf_taxid]].add(leaf_taxid) # Create parent:children structure only for used taxids
			if leaf_taxid==1: break
			leaf_taxid = nodes[leaf_taxid]
	return parents

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

def pre_cluster(rank_lock, leaves, nodes, ranks):
	if rank_lock!="taxid": # if pre-cluster is only by taxid, last step is needed
		# For every leaf, identifies the lineage until the chosen rank
		# If the lineage for the leaf taxid does not have the chosen rank, it will be skipped and pre-cluster by its taxid later
		lineage_merge = defaultdict(set)
		for leaf_taxid, leaf_seqs in leaves.items():
			t = leaf_taxid
			lin = []			
			while t!=1: # runs the tree until finds the chosen rank and save it on the target taxid (taxid from the rank chosen)
				if ranks[t]==rank_lock:
					# union is necessary because many taxids will have the same path and all of them should be united in one taxid
					lineage_merge[t] |= set(lin) # union set operator (same as .extend for list), using set to avoid duplicates
					break 
				lin.append(t)
				t = nodes[t]
		
		# Merge classification on leaves to its parents {chosen rank taxid: [children taxids]}
		for rank_lock_taxid, children_taxids in lineage_merge.items():
			for child_taxid in children_taxids:
				if child_taxid in leaves:
					leaves[rank_lock_taxid].extend(leaves[child_taxid])
					del leaves[child_taxid] # After moving it to the parent, remove it from leaves

	# Pre-cluster leaves by taxid (the ones without the chosen rank will be pre-cluster by its own taxid)
	for leaf_taxid, leaf_seqs in leaves.items():
		sum_length, ids = sum_tuple_ids(leaf_seqs)
		leaves[leaf_taxid] = [(sum_length,*ids)]

	return leaves

def get_unique_rank_taxids(rank_exclusive, leaves, nodes, ranks, start_node):
	unique_rank_taxids = set()
	single_taxids = set()
	if rank_exclusive!="taxid":
		for leaf_taxid in leaves.keys():
			t = leaf_taxid
			while ranks[t]!=rank_exclusive and t!=start_node: t = nodes[t]	
			if t==start_node: # If taxid was not found on the lineage, consider it as single
				single_taxids.add(leaf_taxid)
			else: # if it was found, add to unique rank list
				unique_rank_taxids.add(t)
	else:
		single_taxids = set(leaves.keys())
	return unique_rank_taxids, single_taxids

if __name__ == "__main__":
	main()
