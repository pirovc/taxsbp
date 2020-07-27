#!/usr/bin/env python3

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
import math
import pandas as pd

from collections import defaultdict
from collections import OrderedDict
from pprint import pprint
from pylca.pylca import LCA
from taxsbp.Group import Group
from taxsbp.Cluster import Cluster
from taxsbp.TaxNodes import TaxNodes
from taxsbp.Sequence import Sequence

_taxsbp_silent = True

def main(arguments: str=None):

	if arguments is not None: sys.argv=arguments

	parser = argparse.ArgumentParser(prog='TaxSBP', conflict_handler="resolve", add_help=True)
	parser.add_argument('-i','--input-file', metavar='<input_file>', dest="input_file", help="Tab-separated with the fields: sequence id <tab> sequence length <tab> taxonomic id [<tab> specialization]")
	parser.add_argument('-o','--output-file', metavar='<output_file>', dest="output_file", help="Path to the output tab-separated file with the fields. Default: STDOUT")
	parser.add_argument('-n','--nodes-file', metavar='<nodes_file>', dest="nodes_file", help="nodes.dmp from NCBI Taxonomy")
	parser.add_argument('-m','--merged-file', metavar='<merged_file>', dest="merged_file", help="merged.dmp from NCBI Taxonomy")
	parser.add_argument('-b','--bins', metavar='<bins>', dest="bins", type=int, help="Approximate number of bins (estimated by total length/bin number). [Mutually exclusive -l]")
	parser.add_argument('-l','--bin-len', metavar='<bin_len>', dest="bin_len", type=int, help="Maximum bin length (in bp). Use this parameter insted of -b to define the number of bins. Default: length of the biggest group [Mutually exclusive -b]")
	parser.add_argument('-f','--fragment-len', metavar='<fragment_len>', dest="fragment_len", type=int, default=0, help="Fragment sequences into pieces, output accession will be modified with positions: ACCESION/start:end")
	parser.add_argument('-a','--overlap-len', metavar='<overlap_len>', dest="overlap_len", type=int, default=0, help="Overlap length between fragments [Only valid with -a]")
	parser.add_argument('-p','--pre-cluster', metavar='<pre_cluster>', dest="pre_cluster", type=str, default="", help="Pre-cluster sequences into rank/taxid/specialization, so they won't be splitted among bins [none,specialization name,taxid,species,genus,...] Default: none")
	parser.add_argument('-e','--bin-exclusive', metavar='<bin_exclusive>', dest="bin_exclusive", type=str, default="", help="Make bins rank/taxid/specialization exclusive, so bins won't have mixed sequences. When the chosen rank is not present on a sequence lineage, this sequence will be taxid/specialization exclusive. [none,specialization name,taxid,species,genus,...] Default: none")
	parser.add_argument('-s','--specialization', metavar='<specialization>', dest="specialization", type=str, default="", help="Specialization name (e.g. assembly, strain). If given, TaxSBP will cluster entries on a specialized level after the taxonomic id. The specialization identifier should be provided as an extra collumn in the input_file ans should respect the taxonomic hiercharchy (one taxid -> multiple specializations / one specialization -> one taxid). Default: ''")
	parser.add_argument('-t','--silent', dest="silent", default=False, action='store_true', help="Do not print warning to STDERR")
	parser.add_argument('-u','--update-file', metavar='<update_file>', dest="update_file", type=str, default="", help="Previously generated files to be updated. Default: ''")
	parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0.0')

	if len(sys.argv)<=1: # Print help calling script without parameters
		parser.print_help() 
		return False

	args = parser.parse_args() # read sys.argv[1:] by default

	return pack(**vars(args))

def pack(bin_exclusive: str=None, 
		bin_len: int=0, 
		bins: int=0, 
		fragment_len: int=0, 
		input_file: str=None,
		merged_file: str=None,
		nodes_file: str=None,
		overlap_len: int=0,
		output_file: str=None,
		pre_cluster: str=None,
		specialization: str=None,
		update_file: str=None,
		input_table: pd.DataFrame=None,
		update_table: pd.DataFrame=None,
		silent: bool=True):

	global _taxsbp_silent
	_taxsbp_silent = silent

	if not output_file: output_file = sys.stdout

	special_ranks = ["taxid"] if not specialization else ["taxid", specialization]
	
	taxnodes = TaxNodes(nodes_file, merged_file)
	
	if pre_cluster and pre_cluster not in special_ranks and not taxnodes.has_rank(pre_cluster):
		print_log("Rank for pre-clustering not found: " + pre_cluster)
		print_log("Possible ranks: " + ', '.join(taxnodes.get_ranks()))
		return False
	if bin_exclusive and bin_exclusive not in special_ranks and not taxnodes.has_rank(bin_exclusive):
		print_log("Rank for bin exclusive not found: " + bin_exclusive)
		print_log("Possible ranks: " + ', '.join(taxnodes.get_ranks()))
		return False

	# structure to keep sequence information
	sequences = {}
	

	# Parse from file
	if update_file:
		update_table = parse_update_file(update_file)
	if input_file:
		input_table = parse_input_file(input_file)

	# Check if update table exists, otherwise create an empty
	if update_table is None: update_table = pd.DataFrame()

	# Verify for duplicates and repetead entries
	check_tables(input_table, update_table)

	# Check tax entries
	set_tax_tables(input_table, update_table, taxnodes)

	# Check specialiation and add to taxnodes
	if specialization:
		set_specialization_tax(input_table, update_table, taxnodes, specialization)

	if input_table.empty or (update_file and update_table.empty):
		print_log("No entries to cluster")
		return False

	number_of_bins = 0
	if not update_table.empty:
		groups_bins = process_update_table(update_table, taxnodes, specialization, sequences)
		number_of_bins = len(groups_bins)
		
		# join sequences in the bins
		for binid, group_bin in groups_bins.items(): 
			group_bin.join()
			if bin_exclusive and len(group_bin.get_leaves())>1:
				print_log(binid + " bin with more than one assignment. Is the bin_exclusive rank used to update the same used to generate the bins?")

		# join pre-clustered groups by their LCA or unique leaf
		set_leaf_bins(groups_bins, taxnodes)

	groups = process_input_table(input_table, taxnodes, specialization, sequences)
	del input_table, update_table

	# fragment input
	fragment_groups(groups, sequences, fragment_len, overlap_len)

	# merge bins to current groups
	if update_file: 
		#  add new sequences to the current bins
		for g in groups_bins: groups[g].merge(groups_bins[g])

	# Define bin length
	if bin_len: # user defined
		blen = bin_len
	elif bins: # Estimate bin len based on number of requested bins or direct by user
		blen = sum([g.get_length() for g in groups.values()])/float(bins)
	else: # Default bin length on the max sequence
		blen = max([g.get_length() for g in groups.values()])

	# Pre-clustering
	if pre_cluster: pre_cluster_groups(pre_cluster, groups, taxnodes, specialization)

	# cluster
	final_bins = cluster(groups, taxnodes, blen, bin_exclusive, specialization)

	# sort by bin size
	final_bins.sort(key=lambda tup: tup[0], reverse=True)

	print_results(final_bins, taxnodes, sequences, bin_exclusive, specialization, number_of_bins, output_file)

	return True

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
		final_bins = ApproxSBP("1", groups, children, bin_len)
		
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
				ret.append((sum_length,) + tuple(ids)) # ret.append((sum_length,*ids)) 

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

def set_leaf_bins(groups_bins, taxnodes):
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


def process_input_table(input_table, taxnodes, specialization, sequences: set=None):
	groups = defaultdict(Group)
	for index, row in input_table.iterrows():
		# keep sequence information
		sequences[row['seqid']] = Sequence(row['length'], row['taxid'], row['specialization'], None)
		leaf = row['specialization'] if specialization else row['taxid']
		# Add sequence as cluster and relative leaf nodes
		groups[leaf].add_cluster(leaf, row['seqid'], row['length'])
	return groups

def process_update_table(update_table, taxnodes, specialization, sequences):
	groups = defaultdict(Group)
	for index, row in update_table.iterrows():
		row['seqid'] = make_unique_seqid(row['seqid'], row["seqstart"], row["seqend"])
		# keep sequence information
		sequences[row['seqid']] = Sequence(row['length'], row['taxid'], row['specialization'], row['binid'])
		leaf = row['specialization'] if specialization else row['taxid']
		# Add sequence as cluster and relative leaf nodes
		# binid integer, taxid always string
		groups[row['binid']].add_cluster(leaf, row['seqid'], row['length'])

	return groups

def check_tables(input_table, update_table):
	duplicated_input_idx = input_table.seqid.duplicated()
	if duplicated_input_idx.any():
		print_log("skipped duplicated entries: ")
		print_log(", ".join(input_table[duplicated_input_idx].seqid))	
		# drop duplicated input
		input_table.drop(input_table[duplicated_input_idx].index, inplace=True)

	if not update_table.empty:
		# Do not check for duplicates on update (slip sequences)
		# check for repeated entries on the input/update - remove from input
		repeated_entries_idx = input_table["seqid"].isin(update_table["seqid"])
		if repeated_entries_idx.any():
			#drop repeated
			print_log("skipped repeated entries: ")
			print_log(", ".join(input_table[repeated_entries_idx].seqid))	
			# drop duplicated input
			input_table.drop(input_table[repeated_entries_idx].index, inplace=True)

def set_tax_tables(input_table, update_table, taxnodes):
	# Verify taxonomic entries from the input and taxonomy provided

	# Get unique_taxids from all inputs
	if not update_table.empty:
		unique_taxids = pd.unique(pd.concat([input_table.taxid,update_table.taxid]))
	else:
		unique_taxids = pd.unique(input_table.taxid)

	# Check taxid in nodes and merged
	merged_nodes = {}
	taxids_not_found = set()
	for txid in unique_taxids:
		if not taxnodes.get_parent(txid): 
			m = taxnodes.get_merged(txid)
			if not m: 
				print_log("taxid not found in nodes/merged file (" + txid + ")")
				taxids_not_found.add(txid)
			else:
				print_log("outdated taxid (old: "+txid+" -> new: "+m+")")
				merged_nodes[txid] = m

	# Print invalidated entries
	if taxids_not_found:
		print_log("skipped entries without valid taxid: ")
		print_log(", ".join(input_table[input_table["taxid"].isin(taxids_not_found)].seqid))
		if not update_table.empty:
			print_log(", ".join(update_table[update_table["taxid"].isin(taxids_not_found)].seqid))
		
	# Drop entries without tax
	input_table.drop(input_table[input_table["taxid"].isin(taxids_not_found)].index, inplace=True)
	# Apply merged nodes
	for old_txid, new_txid in merged_nodes.items():
		input_table.loc[input_table["taxid"]==old_txid,'taxid'] = new_txid
	
	# For update table
	if not update_table.empty:
		# Drop entries without tax
		update_table.drop(update_table[update_table["taxid"].isin(taxids_not_found)].index, inplace=True)
		# Apply merged nodes
		for old_txid, new_txid in merged_nodes.items():
			update_table.loc[update_table["taxid"]==old_txid,'taxid'] = new_txid
		
def set_specialization_tax(input_table, update_table, taxnodes, specialization):
	# Set specialiation to taxnodes, dropping inconsistent entries

	# Get unique specialization (taxid + spec) from all inputs
	if not update_table.empty:
		unique_spec = pd.concat([input_table[["taxid","specialization"]],update_table[["taxid","specialization"]]]).drop_duplicates()
	else:
		unique_spec = input_table[["taxid","specialization"]].drop_duplicates()

	#Check if specialization has duplicates
	idx_duplicated = unique_spec.specialization.duplicated()
	if idx_duplicated.any():
		for index, row in unique_spec[idx_duplicated].iterrows():
			print_log("specialization assigned to multiple taxids ("+row["taxid"]+","+row["specialization"]+"), just first taxid-specialization linking will be kept")
		
		invalid_input_table_idx = input_table[["taxid","specialization"]].isin(unique_spec[idx_duplicated]).all(axis=1)
		# Print invalidated entries
		print_log("skipped entries without valid specialization: ")
		print_log(", ".join(input_table[invalid_input_table_idx].seqid))
		# Drop duplicated entries from table
		input_table.drop(input_table[invalid_input_table_idx].index, inplace=True)

		if not update_table.empty:
			invalid_update_table_idx = update_table[["taxid","specialization"]].isin(unique_spec[idx_duplicated]).all(axis=1)
			# Print invalidated entries
			print_log(", ".join(update_table[invalid_update_table_idx].seqid))
			# Drop duplicated entries from table
			update_table.drop(update_table[invalid_update_table_idx].index, inplace=True)

	# Add to taxnodes valid entries (not duplicated)
	for index, row in unique_spec[~idx_duplicated].iterrows():
		taxnodes.add_node(row['taxid'], row['specialization'], specialization) #add taxid as parent, specialization as rank

def fragment_groups(groups, sequences, fragment_len, overlap_len):
	# it will separate into clusters sequences already together

	for leaf,group in groups.items():
		frag_clusters = []
		frag_group = Group()
		for cluster in group.get_clusters():
			for seqid, seqlen in cluster.get_seqlen().items():
				nfrags = 0
				if fragment_len>0: # If fragramentation is required
					nfrags = math.ceil(seqlen / fragment_len) # number of fragments

				if nfrags:
					fragid = ""
					fraglen = 0
					for i in range(nfrags): #range i=0..nfrags
						frag_start = (fragment_len*i) + 1 # +1 to start counting sequence at 1
						
						# if current fragment + overlap is smaller than the total seqlen
						if (fragment_len*(i+1))+overlap_len <= seqlen:
							fraglen = fragment_len+overlap_len # full fragment
						else:
							fraglen = seqlen - (fragment_len*i) # shorted last fragment
							if fraglen<=overlap_len and i>=1: continue # if last fragment is smaller than overlap and its not the first (fraglen<overlap_len)already covered in the last fragment), skip
						
						frag_end = frag_start + fraglen - 1 # -1 offset to count sequences
						fragid=make_unique_seqid(seqid,frag_start,frag_end)
						

						frag_clusters.append(Cluster(fragid,fraglen)) # create cluster for fragment
						# update sequence entry
						sequences[fragid] = Sequence(fraglen, sequences[seqid].taxid, sequences[seqid].specialization, sequences[seqid].binid)
				else:
					fraglen=seqlen
					fragid=seqid+"/1:"+str(seqlen)
					frag_clusters.append(Cluster(fragid,fraglen))
					sequences[fragid] = Sequence(fraglen, sequences[seqid].taxid, sequences[seqid].specialization, sequences[seqid].binid)
			del sequences[seqid] # delete from sequences
		frag_group.add_clusters(group.get_leaves(), frag_clusters)
		groups[leaf] = frag_group #overwrite group

def pre_cluster_groups(pre_cluster_rank, groups, taxnodes, specialization):	
	# if pre-cluster is only by taxid or specialization (leaves), only last step is needed
	if pre_cluster_rank!="taxid" and pre_cluster_rank!=specialization: 
		
		# For every leaf, identifies the lineage until the chosen rank
		# If the lineage for the leaf taxid does not have the chosen rank, it will be skipped and pre-cluster by its taxid later
		lineage_merge = defaultdict(set)
		for leaf in groups:
			t = leaf
			lin = []
			while t!="1": # runs the tree until finds the chosen rank and save it on the target taxid (taxid from the rank chosen)
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
			if t=="1": # If taxid was not found on the lineage, consider it as single
				orphan_taxids.add(leaf)
			else: # if it was found, add to unique rank list
				rank_taxids.add(t)
	else:
		orphan_taxids = set(groups.keys())

	return rank_taxids, orphan_taxids

def print_results(final_bins, taxnodes, sequences, bin_exclusive, specialization, number_of_bins, output_file):

	if output_file!=sys.stdout: output_file=open(output_file,"w")
	new_bins_count=number_of_bins # start bin count in the number of bins (0 if new cluster)
	for bins in final_bins:
		if number_of_bins: #if is updating
			# check the sequences in the bin and check if any already has a binid
			binid = set([sequences[seqid].binid for seqid in bins[1:] if sequences[seqid].binid is not None])
			if not binid: #new bin
				binid = new_bins_count 
				new_bins_count+=1
			elif len(binid)>1: # bins were merge, log
				print_log(str(binid) + " bins were merged. Update parameters differ from ones used to cluster.")
				binid=binid.pop()
			else: # new sequence joined a existing bin
				binid=binid.pop()
		else:
			# create new binid for new bins
			binid = new_bins_count 
			new_bins_count+=1
		
		for seqid in bins[1:]:
			
			# if updating, skip sequences from bins
			if number_of_bins and sequences[seqid].binid is not None: continue

			# if bin_exclusive, recover bin_exclusive rank taxid
			if bin_exclusive and bin_exclusive!="taxid" and bin_exclusive!=specialization:
				taxid = taxnodes.get_rank_node(sequences[seqid].taxid, bin_exclusive)
				if taxid=="1": taxid = sequences[seqid].taxid
			else: # otherwise, just use the input taxid
				taxid = sequences[seqid].taxid

			parsed_seqid, st, en = split_unique_seqid(seqid)
			print(parsed_seqid, st, en, sequences[seqid].seqlen, taxid, str(binid) + ("\t" + sequences[seqid].specialization if specialization else ""), sep="\t", file=output_file)

	if output_file!=sys.stdout: output_file.close()

def print_log(text):
	if not _taxsbp_silent: sys.stderr.write(text+"\n")

def split_unique_seqid(unique_seqid):
	i, pos = unique_seqid.split("/")
	st, en = pos.split(":")
	return i, st, en

def make_unique_seqid(seqid, st, en):
	return seqid+"/"+str(st)+":"+str(en)

def parse_input_file(file):
    colums=['seqid', 'length', 'taxid', 'specialization']
    types={'seqid': 'str', 'length': 'uint64', 'taxid': 'str', 'specialization': 'str'}
    return pd.read_table(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_update_file(file):
    colums=['seqid', 'seqstart', 'seqend', 'length', 'taxid', 'binid', 'specialization']
    types={'seqid': 'str', 'start': 'uint64', 'end': 'uint64', 'length': 'uint64', 'taxid': 'str', 'binid': 'uint64', 'specialization': 'str'}
    return pd.read_table(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

if __name__ == "__main__":
	main()
