from collections import defaultdict

class TaxNodes:

	def __init__(self, nodes_file, merged_file=None):
		self.nodes_file = nodes_file
		self.merged_file = merged_file
		self.nodes, self.ranks = self.parse_nodes_ranks(nodes_file)
		self.merged = self.parse_merged(merged_file) if merged_file else {}

	def add_node(self, parent_taxid, node, rank):
		self.nodes[node] = parent_taxid
		self.ranks[node] = rank

	def parse_nodes_ranks(self, nodes_file):
		# READ nodes -> fields (1:TAXID 2:PARENT_TAXID 3:RANK)
		nodes = {}
		ranks = {}
		with open(nodes_file,'r') as fnodes:
			for line in fnodes:
				taxid, parent_taxid, rank, _ = line.split('\t|\t',3)
				taxid = int(taxid)
				ranks[taxid] = rank
				nodes[taxid] = int(parent_taxid)
		nodes[1] = 0 #Change parent taxid of the root node to 0 (it's usually 1 and causes infinite loop later)
		
		return nodes, ranks

	def parse_merged(self,merged_file):
		# READ nodes -> fields (1:OLD TAXID 2:NEW TAXID)
		merged = {}
		with open(merged_file,'r') as fmerged:
			for line in fmerged:
				old_taxid, new_taxid, _ = line.rstrip().split('\t|',2)
				merged[int(old_taxid)] = int(new_taxid)
		return merged

	def build_children(self, nodes):
		# Define parent->children tree for faster lookup (set unique entries) only for used taxids	
		children = defaultdict(set)
		for node in nodes:
			while True:
				children[self.nodes[node]].add(node) # Create parent:children structure only for used taxids
				if node==1: break # root
				node = self.nodes[node]
		return children

	def get_subtree(self, nodes):
		subtree={}
		for node in nodes:
			t=node
			# while root or branch already on the tree
			while t!=1 and t not in subtree:
				subtree[t] = self.get_parent(t)
				t = self.get_parent(t)
		return subtree

	def get_rank_node(self, node, rank):
		while self.ranks[node]!=rank and node!=1: node = self.nodes[node]
		return node
	
	def get_parent(self, node):
		return self.nodes[node] if node in self.nodes else None

	def get_merged(self, node):
		return self.merged[node] if node in self.merged else None

	def get_ranks(self):
		return set(self.ranks.values())

	def get_rank(self, node):
		return self.ranks[node] if node in self.ranks else None

	def has_rank(self, rank):
		return True if rank in self.get_ranks() else False

