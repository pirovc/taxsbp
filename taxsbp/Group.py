from taxsbp.Cluster import Cluster

class Group:

	def __init__(self):
		self.leaves = set()
		self.elements = []

	def add_cluster(self,leaf_node,seqid,seqlen):
		self.leaves.add(leaf_node)
		self.elements.append(Cluster(seqid,seqlen))

	def add_clusters(self, leaf_node, clusters):
		self.leaves.update(leaf_node)
		self.elements.extend(clusters)

	def get_leaves(self):
		return self.leaves

	def pop_leaf(self):
		return self.leaves.pop()
		
	def get_clusters(self):
		return self.elements

	def get_clusters_bpck(self):
		bins = []
		for c in self.elements:
			bins.append(tuple([c.get_length(),*c.get_ids()]))
		return bins

	def join(self):
		final_cluster = Cluster()
		for c in self.elements:
			final_cluster.update_cluster(c)
		self.elements = [final_cluster]

	def merge(self, group):
		self.add_clusters(group.get_leaves(), group.get_clusters())
		
	def get_length(self):
		return sum([c.get_length() for c in self.elements])

	def get_cluster_count(self):
		return len(self.elements)

	def __repr__(self):
		args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
		return 'Group({})'.format(', '.join(args))