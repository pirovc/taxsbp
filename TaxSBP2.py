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
import sys
from collections import defaultdict
from collections import OrderedDict
from pprint import pprint

from taxsbp.Group import Group
from taxsbp.Cluster import Cluster
from taxsbp.TaxNodes import TaxNodes
from taxsbp.Sequence import Sequence
from scripts.LCA import LCA

def main():
	parser = argparse.ArgumentParser(prog='TaxSBP',conflict_handler="resolve")
	subparsers = parser.add_subparsers()

	# create
	cluster_parser = subparsers.add_parser('cluster', help='Cluster sequences')
	cluster_parser.set_defaults(which='cluster')
	cluster_parser.add_argument('-f', metavar='<input_file>', required=True, dest="input_file", help="Tab-separated with the fields: sequence id <tab> sequence length <tab> taxonomic id [<tab> group]")
	cluster_parser.add_argument('-n', metavar='<nodes_file>', required=True, dest="nodes_file", help="nodes.dmp from NCBI Taxonomy")
	cluster_parser.add_argument('-m', metavar='<merged_file>', dest="merged_file", help="merged.dmp from NCBI Taxonomy")
	cluster_parser.add_argument('-b', metavar='<bins>', default=50, dest="bins", type=int, help="Approximate number of bins (estimated by total length/bin number). Default: 50 [Mutually exclusive -l]")
	cluster_parser.add_argument('-l', metavar='<bin_len>', dest="bin_len", type=int, help="Maximum bin length (in bp). Use this parameter insted of -b to define the number of bins [Mutually exclusive -b]")
	cluster_parser.add_argument('-p', metavar='<pre_cluster>', dest="pre_cluster", type=str, default="", help="Pre-cluster sequences into rank/taxid/specialization, so they won't be splitted among bins [none,specialization name,taxid,species,genus,...] Default: none")
	cluster_parser.add_argument('-r', metavar='<bin_exclusive>', dest="bin_exclusive", type=str, default="", help="Make bins rank/taxid/specialization exclusive, so bins won't have mixed sequences. When the chosen rank is not present on a sequence lineage, this sequence will be taxid/group exclusive. [none,specialization name,taxid,species,genus,...] Default: none")
	cluster_parser.add_argument('-z', metavar='<specialization>', dest="specialization", type=str, default="", help="Specialization name (e.g. assembly, strain). If given, TaxSBP will cluster entries on a specialized level after the taxonomic id. The specialization identifier should be provided as an extra collumn in the input_file ans should respect the taxonomic hiercharchy (one taxid -> multiple specializations / one specialization -> one taxid). Default: ''")
	cluster_parser.add_argument('-u', metavar='<update_file>', dest="update_file", type=str, default="", help="Previously generated files to be updated. Default: ''")

	parser.add_argument('-v', action='version', version='%(prog)s 0.1')
	args = parser.parse_args()

	if len(sys.argv[1:])==0: # Print help calling script without parameters
		parser.print_help() 
		return 0

	special_ranks = ["taxid"] if not args.specialization else ["taxid", args.specialization]
	
	if args.which=="cluster":

		taxnodes = TaxNodes(args.nodes_file, args.merged_file)

		if args.pre_cluster and args.pre_cluster not in special_ranks and not taxnodes.has_rank(args.pre_cluster):
			print("Rank for pre-clustering not found: " + args.pre_cluster)
			print("Possible ranks: " + ', '.join(taxnodes.get_ranks()))
			return
		if args.bin_exclusive and args.bin_exclusive not in special_ranks and not taxnodes.has_rank(args.bin_exclusive):
			print("Rank for bin exclusive not found: " + args.bin_exclusive)
			print("Possible ranks: " + ', '.join(taxnodes.get_ranks()))
			return

		# keep track of sequences read
		sequences = {}

		if args.update_file:
			groups_bins = parse_input(args.update_file, taxnodes, args.specialization, sequences, True)
			# join sequences in the bins
			for group_bin in groups_bins.values(): group_bin.join()
			# join pre-clustered groups by their LCA
			lca_bins(groups_bins, taxnodes)
		groups = parse_input(args.input_file, taxnodes, args.specialization, sequences, False)


#		pprint(groups)
#		print()

		if args.update_file: 
			for g in groups_bins: groups[g].merge(groups_bins[g]) 

#		pprint(groups)
#		print()
		
		# Estimate bin len based on number of requested bins or direct by user
		bin_len = args.bin_len if args.bin_len else sum([g.get_length() for g in groups.values()])/float(args.bins)

		# Pre-clustering
		if args.pre_cluster: pre_cluster(args.pre_cluster, groups, taxnodes, args.specialization)

		final_bins = cluster(groups, taxnodes, bin_len, args.bin_exclusive, args.specialization)

		# sort by bin size
		final_bins.sort(key=lambda tup: tup[0], reverse=True)

		print_results(final_bins, taxnodes, sequences, args.bin_exclusive, args.specialization)

def cluster(groups, taxnodes, bin_len, bin_exclusive, specialization):
	# parent->children structure for fast loookup, only for used taxids 
	children = taxnodes.build_children(groups.keys())

	# bin_exclusive mode
	if bin_exclusive:
		final_bins = []			
		rank_taxids, orphan_taxids = get_rank_taxids(groups, taxnodes, bin_exclusive, specialization)

		if rank_taxids:
			# clustering directly on the rank chosen, recursion required for children nodes
			for rank_taxid in rank_taxids:
				final_bins.extend(ApproxSBP(rank_taxid, groups, children, bin_len))
		if orphan_taxids:
			# clustering directly on the taxid level, no recursion to children nodes necessary
			for orphan_taxid in orphan_taxids:
				final_bins.extend(bpck(groups[orphan_taxid].get_clusters_bpck(), bin_len))

	else: # default mode
		final_bins = ApproxSBP(1, groups, children, bin_len)

	return final_bins	
	
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

def ApproxSBP(v, groups, children, bin_len):
	ch = children[v]
	
	# If it doesn't have any children it's a leaf and should return the packed sequences
	if not ch: return bpck(groups[v].get_clusters_bpck(), bin_len)
		
	# Recursively bin pack children
	# Sort children to keep it more consistent with different versions of the taxonomy (new taxids), str to for groups
	ret = []
	for child in sorted(ch, key=str): ret.extend(ApproxSBP(child, groups, children, bin_len))

	# if current node has sequences assigned to it (but it's not a leaf), add it to the current bin packing (it will first pack with its own children nodes)
	## QUESTION: should I bin together those sequeneces or "distribute" along its children -- command: ret.update(leaves[v]) -- 
	if v in groups: ret.extend(bpck(groups[v].get_clusters_bpck(), bin_len))

	return bpck(ret, bin_len)

def lca_bins(groups_bins, taxnodes):

	#set.union()
	leaves = set()
	for g in groups_bins.values(): leaves.update(g.get_leaves())
	subtree = taxnodes.get_subtree(leaves)
	L = LCA(subtree)

	for binid in list(groups_bins.keys()):
		group_leaves = list(groups_bins[binid].get_leaves())
		if len(group_leaves)>1: # perform if more than one leaf
			lca_node = L(group_leaves[0], group_leaves[1])
			for i in range(len(group_leaves)-2): lca_node = L(lca_node, group_leaves[i])
		else:
			lca_node = group_leaves[0] #only one taxid
		
		groups_bins[lca_node].merge(groups_bins[binid])
		del groups_bins[binid]


def parse_input(input_file, taxnodes, specialization, sequences, bins=False):
	# input file -> fields (0:SEQID 1:LENGTH 2:TAXID [3:SPECIALIZATION] [3-4:BINID])
	groups = defaultdict(Group)
	with open(input_file,'r') as file:
		for line in file:
			try:

				fields = line.rstrip().split('\t')
				seqid= fields[0]
				seqlen = int(fields[1])
				taxid = int(fields[2])
				if specialization:
					spec = fields[3]
					binid = int(fields[4]) if bins else None
				else:
					spec = None
					binid = int(fields[3]) if bins else None

				# validations
				if seqid in sequences:
				 	print_log("[" + seqid + "] skipped - duplicated sequence identifier)")
				 	continue

				if not taxnodes.get_parent(taxid): 
					m = taxnodes.get_merged(taxid)
					if not m: 
						print_log("[" + seqid + "] skipped - taxid not found in nodes/merged file")
						continue
					else:
						taxid = m # Get updated version of the taxid from merged.dmp
				
				if specialization:
					s = taxnodes.get_parent(spec)
					if s!=None and s!=taxid: # group specialization was found in more than one taxid (breaks the tree hiercharchy)
						print_log("[" + seqid + "] skipped - group assigned to multiple taxids, just first taxid-group linking will be considered (" + str(s) + ":" + spec + ")")
						continue
					# update taxonomy
					taxnodes.add_node(taxid, spec, specialization) #add taxid as parent, specialization as rank
					

				# for internal check internal check
				sequences[seqid] = Sequence(seqlen, taxid, spec, binid)

				# Define leaf
				if bins:
					leaf = str(binid) # str to no conflict with taxids
				elif specialization:
					leaf = spec
				else:
					leaf = taxid

				# Add sequence as cluster and relative leaf nodes
				groups[leaf].add_cluster(spec if specialization else taxid,seqid,seqlen)

			except Exception as e:
				print_log(e)

	return groups


def pre_cluster(pre_cluster_rank, groups, taxnodes, specialization):	
	# if pre-cluster is only by taxid or specialization (leaves), only last step is needed
	if pre_cluster_rank!="taxid" and pre_cluster_rank!=specialization: 
		
		# For every leaf, identifies the lineage until the chosen rank
		# If the lineage for the leaf taxid does not have the chosen rank, it will be skipped and pre-cluster by its taxid later
		lineage_merge = defaultdict(set)
		for leaf in groups:
			t = leaf
			lin = []			
			while t!=1: # runs the tree until finds the chosen rank and save it on the target taxid (taxid from the rank chosen)
				if taxnodes.get_rank(t)==pre_cluster_rank:
					# union is necessary because many taxids will have the same path and all of them should be united in one taxid
					lineage_merge[t] |= set(lin) # union set operator (same as .extend for list), using set to avoid duplicates
					break 
				lin.append(t)
				t = taxnodes.get_parent(t)
		
		# Merge classification on leaves to its parents {chosen rank taxid: [children taxids]}
		for rank_taxid, children_taxids in lineage_merge.items():
			for child_taxid in children_taxids:
				if child_taxid in groups:
					# add the clusters from the children clusters to the parent rank taxid
					groups[rank_taxid].merge(groups[child_taxid])
					del groups[child_taxid] # After moving it to the parent, remove it from leaves

	# Pre-cluster leaves by taxid (the ones without the chosen rank will be pre-cluster by its own taxid)
	for group in groups.values():
		group.join()

def get_rank_taxids(groups, taxnodes, bin_exclusive, specialization):
	rank_taxids = set()
	orphan_taxids = set()
	# if not working on leaf level
	if bin_exclusive!="taxid" and bin_exclusive!=specialization:
		for leaf in groups:
			t = taxnodes.get_rank_node(leaf, bin_exclusive)
			if t==1: # If taxid was not found on the lineage, consider it as single
				orphan_taxids.add(leaf)
			else: # if it was found, add to unique rank list
				rank_taxids.add(t)
	else:
		orphan_taxids = set(groups.keys())

	return rank_taxids, orphan_taxids

def print_results(final_bins, taxnodes, sequences, bin_exclusive, specialization):

	for binid,bin in enumerate(final_bins):
		for id in bin[1:]:
			taxid=1
			
			if bin_exclusive and bin_exclusive!="taxid" and bin_exclusive!=specialization:
				taxid = taxnodes.get_rank_node(sequences[id].taxid, bin_exclusive)
			
			if taxid==1: 
				taxid = sequences[id].taxid

			if specialization:
				# Output: accession, seq len, sequence taxid, group, bin
				print(id, sequences[id].seqlen, taxid, sequences[id].specialization, binid, sep="\t")
			else:
				# Output: accession, seq len, sequence taxid, bin
				print(id, sequences[id].seqlen, taxid, binid, sep="\t")

def print_log(text):
	sys.stderr.write(text+"\n")

if __name__ == "__main__":
	main()