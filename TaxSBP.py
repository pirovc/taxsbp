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
from collections import defaultdict

def main():
	parser = argparse.ArgumentParser(prog='TaxSBP')
	parser.add_argument('-a', required=True, metavar='<assembly_reports>', dest="assembly_reports")
	parser.add_argument('-n', required=True, metavar='<nodes>', dest="nodes")
	parser.add_argument('-s', default=2, metavar='<start_node>', dest="start_node", type=int)
	parser.add_argument('-b', default=50, metavar='<bins>', dest="bins", type=int)
	parser.add_argument('-l', metavar='<bin_len>', dest="bin_len", type=int)
	parser.add_argument('-v', action='version', version='%(prog)s 0.01')
	args = parser.parse_args()
	
	global parents
	global leaves
	global bin_len
	
	# READ nodes -> fields (1:TAXID 2:PARENT_TAXID)
	nodes = {}
	with open(args.nodes,'r') as fnodes:
		for line in fnodes:
			taxid, parent_taxid, _ = line.split('\t|\t',2)
			nodes[int(taxid)] = int(parent_taxid)
			
	# READ assembly reports -> fields (6:ACCESSION 8:LENGTH 10:TAXID)
	parents = defaultdict(set)
	leaves = defaultdict(list)
	total_len = 0
	with open(args.assembly_reports,'r') as file:
		for line in file:
			fields = line.split('\t')
			accession = fields[6]
			if accession=="na": continue
			taxid = int(fields[10])
			length = int(fields[8])
			leaves[taxid].append((length,accession))
			while taxid!=1: #Check all taxids in the lineage
				if taxid==args.start_node: total_len+=length # Just account sequence to total when it's on the sub-tree
				parents[nodes[taxid]].add(taxid) # Create parent:children structure only for used taxids
				taxid = nodes[taxid]
	del nodes

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
			print(id,binid,sep="\t")

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

def ApproxSBP_stack(v):
	stack = [v]
	parent_stack = [1]
	bins = defaultdict(list)
	while stack:
		taxid = stack[-1]
		parent_taxid = parent_stack[-1]

		# return to parent
		if taxid==parent_taxid:
			# if current node has sequences assigned to it (but it's not a leaf), add it to the current bin packing (it will first pack with its own children nodes)
			if leaves[taxid]: bins[parent_taxid].extend(bpck(leaves[taxid]))
			stack.pop()
			parent_stack.pop()
			# Add children packs to parent
			bins[parent_stack[-1]].extend(bpck(bins[parent_taxid]))
			del bins[parent_taxid]
			continue
			
		# Check if node has children
		children = parents[taxid]
		
		# If it doesn't have any children it's a leaf and should return the packed sequences
		if not children: 
			bins[parent_taxid].extend(bpck(leaves[taxid]))
			stack.pop()
		else:
			parent_stack.append(taxid)
			stack.extend(sorted(children,reverse=True)) # Sort children (reversed because it's in a stack) to keep it more consistent with different versions of the taxonomy (new taxids)

	return bins[1]
	
if __name__ == "__main__":
	main()