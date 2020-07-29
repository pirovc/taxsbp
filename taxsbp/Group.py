from taxsbp.Cluster import Cluster

class Group:

	def __init__(self):
		self.leaves = set()
		self.elements = []

	def add_cluster(self,leaf_node,seqid,seqlen,binid: int=None):
		self.leaves.add(leaf_node)
		self.elements.append(Cluster(seqid,seqlen,binid))

	def add_clusters(self, leaf_nodes, clusters):
		self.leaves.update(leaf_nodes)
		self.elements.extend(clusters)

	def get_leaves(self):
		return self.leaves

	def get_clusters(self):
		return self.elements

	def clear_elements(self):
		self.elements = []

	def get_clusters_to_bpck(self):
		ret = []
		for c in self.elements:
			ret.extend(c.get_tuples())
		return ret

	def add_clusters_from_bpck(self, clusters, leaves: set=None):
		for cluster in clusters:
			c = Cluster()
			for seqlen,seqid in cluster:
				c.add(seqid,seqlen)
			self.elements.append(c)
			if leaves: self.leaves.update(leaves)

	def join_clusters(self):
		# Join clusters in the group by binid
		final_clusters = {}
		for c in self.elements:
			if c.binid not in final_clusters: final_clusters[c.binid] = Cluster(binid=c.binid)
			final_clusters[c.binid].update(c)
		self.elements = list(final_clusters.values())

	def merge(self, group):
		self.add_clusters(group.get_leaves(), group.get_clusters())
		
	def get_length(self):
		return sum([c.get_cluster_length() for c in self.elements])

	def get_cluster_count(self):
		return len(self.elements)

	def print_bins(self):
		for c in self.elements:
			print(c.seqlen)

	def __repr__(self):
		args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
		return 'Group({})'.format(', '.join(args))