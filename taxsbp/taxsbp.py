#!/usr/bin/env python3

# The MIT License (MIT)
# 
# Copyright (c) 2020 - Vitor C. Pir  - pirovc@posteo.net
#
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
from math import ceil
from io import StringIO
from collections import defaultdict
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
	parser.add_argument('-u','--update-file', metavar='<update_file>', dest="update_file", type=str, default="", help="Previously generated files to be updated. Output only new sequences. Default: ''")
	parser.add_argument('-w','--allow-merge', dest="allow_merge", default=False, action='store_true', help="When updating, allow merging of existing bins. Will output the whole set, not only new bins")
	parser.add_argument('-t','--silent', dest="silent", default=False, action='store_true', help="Do not print warning to STDERR")
	parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0.0')

	if len(sys.argv)<=1: # Print help calling script without parameters
		parser.print_help() 
		return False

	args = parser.parse_args() # read sys.argv[1:] by default

	# if calling from cli and no output, set to stdout
	if not args.output_file: args.output_file=sys.stdout
	return pack(**vars(args))

def pack(bin_exclusive: str=None, 
		bin_len: int=0, 
		bins: int=0, 
		fragment_len: int=0, 
		input_file: str=None,
		input_table: str=None,
		merged_file: str=None,
		nodes_file: str=None,
		overlap_len: int=0,
		output_file: str=None,
		pre_cluster: str=None,
		specialization: str=None,
		update_file: str=None,
		update_table: str=None,
		silent: bool=True,
		allow_merge: bool=False):

	global _taxsbp_silent
	_taxsbp_silent = silent

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
	# keep track of sequences read
	control_seqid = set()
	# main structure to organize groups and clusters
	groups = defaultdict(Group)

	used_bins = set()
	if update_file or update_table:
		# Parse update file to groups/sequences 
		# return grouping by their binid: groups[binid] = Group(...)
		parse_input(update_file, update_table, groups, taxnodes, specialization, sequences, control_seqid, fragment_len, overlap_len, bins=True)
		# get used bins
		used_bins = set(groups.keys())

		# Join clusters of for each group (they should not be splitted)
		for binid, group in groups.items(): 
			group.join_clusters()
			if bin_exclusive and len(group.get_leaves())>1:
				print_log(str(binid) + " bin with more than one leaf assignment. Clusters are not bin exclusive.")
		# replace binid of groups by their LCA or unique leaf
		set_leaf_bins(groups, taxnodes)

	# Parse input files and add to groups/sequences
	parse_input(input_file, input_table, groups, taxnodes, specialization, sequences, control_seqid, fragment_len, overlap_len, bins=False)

	if not groups:
		print_log("No entries to cluster")
		return False

	# Define bin length
	if bin_len: # user defined
		blen = bin_len
	elif bins: # Estimate bin len based on number of requested bins or direct by user
		blen = sum([g.get_length() for g in groups.values()])/float(bins)
	else: # Default bin length on the max sequence
		blen = max([g.get_length() for g in groups.values()])

	# Pre-clustering
	if pre_cluster: pre_cluster_groups(pre_cluster, groups, taxnodes, specialization)

	# Remove binid from clusters to allow them to merge
	# Don't do it before so they won't be joined together (Group func. join_clusters())
	if (update_file or update_table) and allow_merge: clear_binids(groups)

	# cluster
	cluster(groups, taxnodes, blen, bin_exclusive, specialization)

	# Set tax entries to their bin_exclusive 
	if bin_exclusive and bin_exclusive!="taxid" and bin_exclusive!=specialization:
		set_tax(sequences, taxnodes, bin_exclusive)

	# Set binids for groups
	set_bins(groups, sequences, used_bins, allow_merge)
	
	# yield a generator for results
	res = generate_results(groups, sequences, specialization, allow_merge)

	if not output_file:
		# Call from python without output_file set
		return [r for r in res]
	elif output_file==sys.stdout:
		# call from cli without output_file set
		for r in res:
			print(*r, sep="\t")	
	else:
		# call with output_file set
		with open(output_file,"w") as file:
			for r in res:
				print(*r, sep="\t", file=file)	

	return True

def cluster(groups, taxnodes, bin_len, bin_exclusive, specialization):
	# parent->children structure for fast loookup, only for used taxids 
	children = taxnodes.build_children(groups.keys())

	# bin_exclusive mode
	if bin_exclusive:	
		rank_taxids, orphan_taxids = get_rank_taxids(groups, taxnodes, bin_exclusive, specialization)
		if rank_taxids:
			# clustering directly on the rank chosen, recursion required for children nodes
			for rank_taxid in rank_taxids:
				ApproxSBP(rank_taxid, None, groups, children, bin_len)
		if orphan_taxids:
			# clustering directly on the taxid level, no recursion to children nodes necessary
			for orphan_taxid in orphan_taxids:
				bpck(groups, orphan_taxid, orphan_taxid, bin_len)

	else: # default mode
		ApproxSBP("1", None, groups, children, bin_len)

def bpck(groups, node, parent, bin_len): 
	# Perform bin packing on a single node
	# it packs the clusters on groups[node] and add to groups[parent]
	# if node and parent are equal, root was reached
	at_root = True if node==parent else False

	# If there is only one cluster, do not need to pack
	if groups[node].get_cluster_count()==1:
		if not at_root: # transfer cluster to parent if not root
			groups[parent].merge(groups[node])
			del groups[node]
	else:
		# Perform bin packing
		clusters = binpacking.to_constant_volume(groups[node].get_clusters_to_bpck(), bin_len, weight_pos=1)
		if clusters:
			if not at_root:
				# Parse clustered results into parent node and remove actual node
				groups[parent].add_clusters_from_bpck(clusters, leaves=groups[node].get_leaves())
				del groups[node]
			else: #if root
				# Parse clustered results into same node (clear it before)
				groups[parent].clear_clusters()
				groups[parent].add_clusters_from_bpck(clusters)

def ApproxSBP(node, parent, groups, children, bin_len):
	# Function to perform hiearchical bin packing recursively
	# Recursively call to pack sorted list of children (to get always same results)
	for child in sorted(children[node], key=str): 
		ApproxSBP(child, node, groups, children, bin_len)
	else:
		# If node is a leaf - no child in children[node]
		# or
		# After all children of a node were packed in the for loop, pack node itself into parent
		bpck(groups, node, parent if parent is not None else node, bin_len)
	
def set_leaf_bins(groups, taxnodes):
	leaves = set()
	for g in groups.values(): leaves.update(g.get_leaves())
	subtree = taxnodes.get_subtree(leaves)
	L = LCA(subtree)

	for binid in list(groups.keys()):
		group_leaves = list(groups[binid].get_leaves())
		if len(group_leaves)>1: # perform if more than one leaf
			lca_node = L(group_leaves[0], group_leaves[1])
			for i in range(len(group_leaves)-2): lca_node = L(lca_node, group_leaves[i])
		else:
			lca_node = group_leaves[0] #only one taxid
		
		groups[lca_node].merge(groups[binid])
		del groups[binid]


def parse_input(file, table, groups, taxnodes, specialization, sequences, control_seqid, fragment_len, overlap_len, bins: bool=False):
	# Parser for input_file or update_file
	#input_file: 0:SEQID 1:LENGTH 2:TAXID [3:SPECIALIZATION] 
	#update_file: 0:SEQID 1:SEQSTART 2:SEQEND 3:LENGTH 4:TAXID 5:BINID [6:SPECIALIZATION]

	fields_pos = {"seqid": 0}
	if bins:
		fields_pos["seqstart"] = 1
		fields_pos["seqend"] = 2
		fields_pos["seqlen"] = 3
		fields_pos["taxid"] = 4
		fields_pos["binid"] = 5
		fields_pos["specialization"] = 6
	else:
		fields_pos["seqlen"] = 1
		fields_pos["taxid"] = 2
		fields_pos["specialization"] = 3

	if table is not None:
		inf = StringIO(table)
	else:
		inf = open(file,'r')

	for line in inf:
		try:
			fields = line.rstrip().split('\t')
			seqid = fields[fields_pos["seqid"]]
			seqlen = int(fields[fields_pos["seqlen"]])
			taxid = fields[fields_pos["taxid"]]
			binid = int(fields[fields_pos["binid"]]) if bins else None
			spec = fields[fields_pos["specialization"]] if specialization else None

			# if reading main input (after loaded bins), do not add duplicated sequences
			if not bins and seqid in control_seqid:
			 	print_log("[" + seqid + "] skipped - duplicated sequence identifier")
			 	continue

			# add entry before other checks
			# when parsing bins, do not ignore the ones with failing tax/spec
			# since they have to be kept
			control_seqid.add(seqid)

			if not taxnodes.get_parent(taxid): 
				m = taxnodes.get_merged(taxid)
				if not m: 
					print_log("[" + seqid + "] skipped - taxid not found in nodes/merged file")
					continue
				else:
					print_log("[" + seqid + "] outdated taxid (old: "+str(taxid)+" -> new:"+str(m)+")")
					taxid = m # Get updated version of the taxid from merged.dmp
			
			if specialization:
				s = taxnodes.get_parent(spec)
				if s!=None and s!=taxid: # group specialization was found in more than one taxid (breaks the tree hiercharchy)
					print_log("[" + seqid + "] skipped - specialization assigned to multiple leaves, just first leaf-group linking will be considered (" + str(s) + ":" + spec + ")")
					continue
				# update taxonomy
				taxnodes.add_node(taxid, spec, specialization) #add taxid as parent, specialization as rank

			leaf = spec if specialization else taxid
			# Define leaf
			if bins:
				seqid = make_unique_seqid(seqid, fields[fields_pos["seqstart"]], fields[fields_pos["seqend"]])
				sequences[seqid] = Sequence(seqlen, taxid, spec, binid)
				# Use binid as groupid
				# Do not save binid information if bins can be merged
				groups[binid].add_cluster(leaf,seqid,seqlen,binid)
			else:
				# Fragment input, add to sequences and clusters
				groups[leaf].add_clusters([leaf], fragment_input(seqid, seqlen, taxid, spec, fragment_len, overlap_len, sequences))
		except Exception as e:
			print_log(str(e))
	inf.close()

def fragment_input(seqid, seqlen, taxid, specialization, fragment_len, overlap_len, sequences):
	frag_clusters = []
	
	nfrags = 0
	if fragment_len>0: # If fragramentation is required
		nfrags = ceil(seqlen / fragment_len) # number of fragments

	if nfrags:
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
			frag_clusters.append(Cluster([fragid],fraglen)) # create cluster for fragment
			sequences[fragid] = Sequence(fraglen, taxid, specialization, None)
	else:
		fraglen=seqlen
		fragid=seqid+"/1:"+str(seqlen)
		frag_clusters.append(Cluster([fragid],fraglen))
		sequences[fragid] = Sequence(fraglen, taxid, specialization, None)

	return frag_clusters

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
		group.join_clusters()

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

def clear_binids(groups):
	for group in groups.values():
		for c in group.get_clusters():
			c.set_binid(None)

def set_tax(sequences, taxnodes, bin_exclusive):
	for seqid in sequences:
		taxid = taxnodes.get_rank_node(sequences[seqid].taxid, bin_exclusive)
		if taxid!="1": sequences[seqid].taxid = taxid

def set_bins(groups, sequences, used_bins, allow_merge):
	# Initialize bin count
	update=False
	free_binids = set()
	if used_bins: # there are already previously generated bins
		#place the binid count on last bin
		binid_count=max(used_bins)
		# Create a set of free binids in between the last binid
		free_binids.update(set(range(max(used_bins)+1)).difference(used_bins))
		# define update mode
		update=True
	else:
		binid_count=-1

	for v, group in groups.items():
		for cluster in group.get_clusters():
			binid=None
			#if is updating existing bins
			if update: 
				# if merging bins is not allowed
				if not allow_merge:
					# If binid is assigned, use it or create new bin
					bin_cluster=cluster.get_binid()
					if bin_cluster is not None:
						binid = bin_cluster
					else:
						# If some binid is free, use it
						if free_binids:
							binid = free_binids.pop()
						else:
							# create new bin
							binid_count+=1
							binid=binid_count
				else:
					# check the sequences in the bin and check if any already has a binid
					existing_binids = set([sequences[seqid].binid for seqid in cluster.get_ids() if sequences[seqid].binid is not None])
					if not existing_binids: 
						# If some binid is free, use it
						if free_binids:
							binid = free_binids.pop()
						else:
							binid_count+=1
							binid=binid_count
					elif len(existing_binids)==1: # cluster with only one assigned bin
						binid=existing_binids.pop()
					else: # more than one bin per cluster
						# get the smallest
						binid=min(existing_binids) 
						# report merging
						print_log("bins " + ", ".join([str(b) for b in existing_binids]) + " were merged into bin " + str(binid))
						# remove it from list
						existing_binids.remove(binid)
						# set remaining as free
						free_binids.update(existing_binids) 
			else: # cluster new file
				binid_count+=1
				binid=binid_count

			# set bin for cluster
			cluster.set_binid(binid)

	if len(free_binids)>0:
		print_log(", ".join([str(b) for b in free_binids])  + " bins were left empty")
					
def generate_results(groups, sequences, specialization, allow_merge):
	for v, group in groups.items():
		for cluster in group.get_clusters():
			for seqid in cluster.get_ids():
				# if updating and not merging, do not print older entries
				if sequences[seqid].binid is not None and not allow_merge: continue
				parsed_seqid, st, en = split_unique_seqid(seqid)
				spec = sequences[seqid].specialization if specialization else ""
				yield [parsed_seqid, st, en, sequences[seqid].seqlen, sequences[seqid].taxid, str(cluster.get_binid()), spec]

def print_log(text):
	if not _taxsbp_silent: sys.stderr.write(text+"\n")

def split_unique_seqid(unique_seqid):
	i, pos = unique_seqid.split("/")
	st, en = pos.split(":")
	return i, st, en

def make_unique_seqid(seqid, st, en):
	return seqid+"/"+str(st)+":"+str(en)

if __name__ == "__main__":
	main()
