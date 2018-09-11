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
	create_parser.add_argument('-f', metavar='<input_file>', required=True, dest="input_file", help="Tab-separated with the fields: sequence id <tab> sequence length <tab> taxonomic id [<tab> group]")
	create_parser.add_argument('-n', metavar='<nodes_file>', required=True, dest="nodes_file", help="nodes.dmp from NCBI Taxonomy")
	create_parser.add_argument('-m', metavar='<merged_file>', dest="merged_file", help="merged.dmp from NCBI Taxonomy")
	create_parser.add_argument('-b', metavar='<bins>', default=50, dest="bins", type=int, help="Approximate number of bins (estimated by total length/bin number). Default: 50 [Mutually exclusive -l]")
	create_parser.add_argument('-l', metavar='<bin_len>', dest="bin_len", type=int, help="Maximum bin length (in bp). Use this parameter insted of -b to define the number of bins [Mutually exclusive -b]")
	create_parser.add_argument('-p', metavar='<pre_cluster>', dest="pre_cluster", type=str, default="", help="Pre-cluster sequences into rank/taxid/specialization, so they won't be splitted among bins [none,specialization name,taxid,species,genus,...] Default: none")
	create_parser.add_argument('-r', metavar='<bin_exclusive>', dest="bin_exclusive", type=str, default="", help="Make bins rank/taxid/specialization exclusive, so bins won't have mixed sequences. When the chosen rank is not present on a sequence lineage, this sequence will be taxid/group exclusive. [none,specialization name,taxid,species,genus,...] Default: none")
	create_parser.add_argument('-z', metavar='<specialization>', dest="specialization", type=str, default="", help="Specialization name (e.g. assembly, strain). If given, TaxSBP will cluster entries on a specialized level after the taxonomic id. The specialization identifier should be provided as an extra collumn in the input_file ans should respect the taxonomic hiercharchy (one taxid -> multiple specializations / one specialization -> one taxid). Default: ''")
	create_parser.add_argument('--sorted-output', dest="sorted_output", default=False, action='store_true', help="Bins will be created based on the bin size in descending order")
	
	# add
	add_parser = subparsers.add_parser('add', help='add new bins')
	add_parser.set_defaults(which='add')
	add_parser.add_argument('-f', metavar='<input_file>', required=True, dest="input_file", help="Tab-separated with the fields: sequence id <tab> sequence length <tab> taxonomic id [<tab> group]")
	add_parser.add_argument('-n', metavar='<nodes_file>', required=True, dest="nodes_file", help="nodes.dmp from NCBI Taxonomy")
	add_parser.add_argument('-m', metavar='<merged_file>', dest="merged_file", help="merged.dmp from NCBI Taxonomy")
	add_parser.add_argument('-b', metavar='<bins>', default=50, dest="bins", type=int, help="Approximate number of bins (estimated by total length/bin number). Default: 50 [Mutually exclusive -l]")
	add_parser.add_argument('-l', metavar='<bin_len>', dest="bin_len", type=int, help="Maximum bin length (in bp). Use this parameter insted of -b to define the number of bins [Mutually exclusive -b]")
	add_parser.add_argument('-p', metavar='<pre_cluster>', dest="pre_cluster", type=str, default="", help="Pre-cluster sequences into rank/taxid/specialization, so they won't be splitted among bins [none,specialization name,taxid,species,genus,...] Default: none")
	add_parser.add_argument('-r', metavar='<bin_exclusive>', dest="bin_exclusive", type=str, default="", help="Make bins rank/taxid/specialization exclusive, so bins won't have mixed sequences. When the chosen rank is not present on a sequence lineage, this sequence will be taxid/group exclusive. [none,specialization name,taxid,species,genus,...] Default: none")
	add_parser.add_argument('-z', metavar='<specialization>', dest="specialization", type=str, default="", help="Specialization name (e.g. assembly, strain). If given, TaxSBP will cluster entries on a specialized level after the taxonomic id. The specialization identifier should be provided as an extra collumn in the input_file ans should respect the taxonomic hiercharchy (one taxid -> multiple specializations / one specialization -> one taxid). Default: ''")
	add_parser.add_argument('-u', metavar='<update>', dest="update", type=str, default="", help="")
	add_parser.add_argument('--sorted-output', dest="sorted_output", default=False, action='store_true', help="Bins will be created based on the bin size in descending order")

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
	
	parser.add_argument('-v', action='version', version='%(prog)s 0.1')
	args = parser.parse_args()

	if len(sys.argv[1:])==0: # Print help calling script without parameters
		parser.print_help() 
		return 0

	if args.which=="create" or args.which=="add":
		accessions_bins=None
		if args.which=="add": # read  bins
			# leaves_bins are pre-clustered by binid (so they cannot be splited again)
			leaves_bins, accessions_bins, binid_map, total_len_bins = read_bins(args.update, args.specialization, args.bin_exclusive)
			number_of_bins = len(leaves_bins)

		nodes, ranks = read_nodes(args.nodes_file)
		merged = read_merged(args.merged_file)
		leaves, accessions, total_len, nodes, ranks = read_input(args.input_file, nodes, merged, ranks, args.specialization, accessions_bins)

		possible_ranks = set(ranks.values())
		if args.pre_cluster and args.pre_cluster!="taxid" and args.pre_cluster not in possible_ranks:
			print("Rank for pre-clustering not found: " + args.pre_cluster)
			print("Possible ranks: " + ', '.join(possible_ranks))
			return
		if args.bin_exclusive and args.bin_exclusive!="taxid" and args.bin_exclusive not in possible_ranks:
			print("Rank for bin exclusive not found: " + args.bin_exclusive)
			print("Possible ranks: " + ', '.join(possible_ranks))
			return

		# Pre-clustering
		if args.pre_cluster: leaves = pre_cluster(args.pre_cluster, leaves, nodes, ranks, args.specialization)
			
		# build parent structure
		parents = parents_tree(leaves, nodes)

		# Bin length (estimated from number of bins or directly)
		if args.bin_len:
			bin_len = args.bin_len
		else:
			# Estimate bin len based on number of requested bins
			bin_len = total_len/float(args.bins) 

	if args.which=="create":
		rank_taxid_leaf = ""
		if not args.bin_exclusive: # Standard mode
			# Run taxonomic structured bin packing for the whole tree
			final_bins = ApproxSBP(1, leaves, parents, bin_len) # RECURSIVE
		else: # Bin exclusive mode
			final_bins, rank_taxid_leaf = cluster_bin_exclusive(args.bin_exclusive, leaves, nodes, parents, ranks, bin_len, args.specialization)

		if args.sorted_output: final_bins.sort(key=lambda tup: tup[0], reverse=True)

		# Print resuls (by sequence)
		for binid,bin in enumerate(final_bins):
			for id in bin[1:]:
				print_results(id, binid, accessions, rank_taxid_leaf, nodes, args.specialization, args.bin_exclusive)

	elif args.which=="add":


		rank_taxid_leaf = ""
		if not args.bin_exclusive:
			used_tree = {}
			for leaf_taxid in set.union(*binid_map.values()):
				t = leaf_taxid
				while t!=1:
					used_tree[t] = nodes[t]
					t = nodes[t]
					if t in used_tree: break # branch already in the dict

			from LCA import LCA
			L = LCA(used_tree)
			for binid in list(leaves_bins.keys()):
				if len(binid_map[binid])>1: # perform LCA
					t = L(binid_map[binid].pop(),binid_map[binid].pop())
					for i in range(len(binid_map[binid])): t = L(t, binid_map[binid].pop())
				else:
					t = binid_map[binid].pop() #only one taxid
				leaves_bins[t].extend(leaves_bins.pop(binid))

			# add bins to the new leaves to cluster together
			for lca_taxid in leaves_bins:
				leaves[lca_taxid].extend(leaves_bins[lca_taxid])

			# redine parents			
			parents = parents_tree(leaves, nodes)
			final_bins = ApproxSBP(1, leaves, parents, bin_len) # RECURSIVE
		else:
			# Convert leaves_bins[binid] to leaves_bins[taxid], merging bins with same taxid
			for binid in list(leaves_bins.keys()):
				t = list(binid_map[binid])[0] # get taxonomic assignment/specialization for the bin (unique for bin_exclusive)
				leaves_bins[t].extend(leaves_bins.pop(binid))

			# if not working on leafs, define taxid for the bin_exclusive rank
			if args.bin_exclusive!="taxid" and args.bin_exclusive!=args.specialization:
				for taxid in list(leaves.keys()):
					rank_taxid = get_rank_taxid(taxid, args.bin_exclusive, nodes, ranks)
					if rank_taxid!=1 and rank_taxid!=taxid:
						leaves[rank_taxid].extend(leaves.pop(taxid))

			# add bins to the new leaves when possible
			for taxid in list(leaves.keys()):
				if taxid in leaves_bins:
					leaves[taxid].extend(leaves_bins[taxid])
			
			# redine parents			
			parents = parents_tree(leaves, nodes)			
			final_bins, rank_taxid_leaf = cluster_bin_exclusive(args.bin_exclusive, leaves, nodes, parents, ranks, bin_len, args.specialization)
			
		if args.sorted_output: final_bins.sort(key=lambda tup: tup[0], reverse=True)

		# print new sequences
		for bin in final_bins:
			# Define if it's a new bin or sequences were incorporated to existing bins
			binid = ""
			for id in bin[1:]: # check if bin belongs to old idons
				# check if entry in cluster was in some bins (meaning it's a known bin)
				if id in accessions_bins: 
					binid = accessions_bins[id][2] # get old bin number
					break

			# no sequences in the bin were previously found, create new bin for the cluster
			if not binid:
				binid = number_of_bins
				number_of_bins += 1

			for id in bin[1:]:
				# Just print new sequences, remove to print all affected
				if id in accessions:
					print_results(id, binid, accessions, rank_taxid_leaf, nodes, args.specialization, args.bin_exclusive)
					

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
	if v in leaves: ret.extend(bpck(leaves[v], bin_len))

	return bpck(ret, bin_len)

def read_nodes(nodes_file):
	# READ nodes -> fields (1:TAXID 2:PARENT_TAXID)
	nodes = {}
	ranks = {}
	with open(nodes_file,'r') as fnodes:
		for line in fnodes:
			taxid, parent_taxid, rank, _ = line.split('\t|\t',3)
			ranks[int(taxid)] = rank
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

def read_bins(input_file, specialization, bin_exclusive):
	# READ bins file -> fields (0:ACCESSION 1:LENGTH 2:TAXID [SPECIALIZATION] 3:BINID)
	
	leaves_bins = defaultdict(list)
	accessions_bins = dict()
	binid_map = defaultdict(set)
	total_len_bins = 0

	with open(input_file,'r') as file:
		for line in file:
			fields = line.rstrip().split('\t')
			accession = fields[0]
			length = int(fields[1])
			taxid = int(fields[2])
			if specialization:
				spec = fields[3]
				binid = int(fields[4])
			else:
				binid = int(fields[3])
				
			if bin_exclusive==specialization and (bin_exclusive or specialization):
				accessions_bins[accession] = (length,spec,binid) 
				binid_map[binid].add(spec)
			else:
				accessions_bins[accession] = (length,taxid,binid) 
				binid_map[binid].add(taxid)

			leaves_bins[binid].append((length,accession)) # add as leaf
			total_len_bins+=length

	# pre-cluster bins by their binid, so they won't split
	leaves_bins = pre_cluster("binid", leaves_bins, None, None, "binid")

	return leaves_bins, accessions_bins, binid_map, total_len_bins

def read_input(input_file, nodes, merged, ranks, specialization, accessions_bins=None):
	# READ input file -> fields (0:ACCESSION 1:LENGTH 2:TAXID [3:SPECIALIZATION])
	
	leaves = defaultdict(list)
	accessions = dict()
	total_len = 0

	with open(input_file,'r') as file:
		for line in file:
			fields = line.rstrip().split('\t')
			accession = fields[0]
			if accession in accessions:
				print_log("[" + accession + "] skipped - duplicated sequence accession")
				continue
			if accessions_bins:
				if accession in accessions_bins:
					print_log("[" + accession + "] skipped - duplicated sequence accession already clustered")
					continue
			length = int(fields[1])
			taxid = int(fields[2])
			if taxid not in nodes: 
				if taxid not in merged: 
					print_log("[" + "\t".join(fields) + "] skipped - taxid not found in nodes/merged file")
					continue
				else:
					taxid = merged[taxid] # Get new taxid from merged.dmp
		
			total_len+=length # Account for seq. length
			if specialization:
				spec = fields[3]
				if spec in nodes and nodes[spec]!=taxid: # group specialization was found in more than one taxid (breaks the tree hiercharchy)
					print_log("[" + "\t".join(fields) + "] skipped - group assigned to multiple taxids, just first taxid-group linking will be considered (" + str(nodes[spec]) + ":" + spec + ")")
					continue
				# update taxonomy
				nodes[spec] = taxid # Update nodes with new leaf
				ranks[spec] = specialization
				leaves[spec].append((length,accession)) # add group as leaf
				accessions[accession] = (length,spec) 
			else:
				leaves[taxid].append((length,accession)) # Keep length and accession for each taxid (multiple entries)
				accessions[accession] = (length,taxid) # Keep length and taxid for each accession (input file)

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

def pre_cluster(rank_lock, leaves, nodes, ranks, specialization):
	
	if rank_lock!="taxid" and rank_lock!=specialization: # if pre-cluster is only by taxid or specialization, only last step is needed
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

def cluster_bin_exclusive(bin_exclusive, leaves, nodes, parents, ranks, bin_len, specialization):
	# Run taxonomic structured bin packing for each rank/taxid separetly
	final_bins = []
	rank_taxid_leaf = {} # keep taxid of the chosen rank
		
	unique_rank_taxids, single_taxids = get_unique_rank_taxids(bin_exclusive, leaves, nodes, ranks, specialization)

	if unique_rank_taxids:
		# clustering directly on the rank chosen, recursion required for children nodes
		for unique_rank_taxid in unique_rank_taxids:
			res = ApproxSBP(unique_rank_taxid, leaves, parents, bin_len)
			for bin in res: #Save taxid of the chosen rank for output
				for id in bin[1:]:
					rank_taxid_leaf[id] = unique_rank_taxid
			final_bins.extend(res)

	if single_taxids:
		# clustering directly on the taxid level, no recursion to children nodes necessary
		for single_taxid in single_taxids:
			res = bpck(leaves[single_taxid], bin_len)
			for bin in res: #Save taxid of the chosen rank for output
				for id in bin[1:]:
					if specialization: #If on group mode, get taxid not group leaf
						rank_taxid_leaf[id] = nodes[single_taxid]
					else:
						rank_taxid_leaf[id] = single_taxid
			final_bins.extend(res)

	return final_bins, rank_taxid_leaf

def get_unique_rank_taxids(bin_exclusive, leaves, nodes, ranks, specialization):
	unique_rank_taxids = set()
	single_taxids = set()
	if bin_exclusive!="taxid" and bin_exclusive!=specialization:
		for leaf_taxid in leaves.keys():
			t = get_rank_taxid(leaf_taxid, bin_exclusive, nodes, ranks)
			if t==1: # If taxid was not found on the lineage, consider it as single
				single_taxids.add(leaf_taxid)
			else: # if it was found, add to unique rank list
				unique_rank_taxids.add(t)
	else:
		single_taxids = set(leaves.keys())

	return unique_rank_taxids, single_taxids

def get_rank_taxid(taxid, rank, nodes, ranks):
	while ranks[taxid]!=rank and taxid!=1: taxid = nodes[taxid]
	return taxid

def print_results(id, binid, accessions, rank_taxid_leaf, nodes, specialization, bin_exclusive):
	if specialization:
		if bin_exclusive and bin_exclusive != specialization: # by some taxonomic rank
			# Output: accession, seq len, bin exclusive taxid, group, bin
			print(id, accessions[id][0], rank_taxid_leaf[id], accessions[id][1], binid, sep="\t")
		else:
			# Output: accession, seq len, sequence taxid, group, bin
			print(id, accessions[id][0], nodes[accessions[id][1]], accessions[id][1], binid, sep="\t")
	else:
		if bin_exclusive: # by some taxonomic rank
			# Output: accession, seq len, bin exclusive taxid, bin
			print(id, accessions[id][0], rank_taxid_leaf[id], binid, sep="\t")
		else:
			# Output: accession, seq len, sequence taxid, bin
			print(id, accessions[id][0], accessions[id][1], binid, sep="\t")


def print_log(text):
	sys.stderr.write(text+"\n")

if __name__ == "__main__":
	main()
