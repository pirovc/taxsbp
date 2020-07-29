from taxsbp.Cluster import Cluster

class Group:

	def __init__(self):
		self.leaves = set()
		self.clusters = []

	def add_cluster(self,leaf_node,seqid,seqlen,binid: int=None):
		self.leaves.add(leaf_node)
		self.clusters.append(Cluster([seqid],seqlen,binid))

	def add_clusters(self, leaf_nodes, new_clusters):
		self.leaves.update(leaf_nodes)
		self.clusters.extend(new_clusters)

	def get_leaves(self):
		return self.leaves

	def get_clusters(self):
		return self.clusters

	def clear_clusters(self):
		self.clusters = []

	def get_clusters_to_bpck(self):
		ret = []
		for c in self.clusters:
			ret.append(c.get_tuples())
		return ret

	def add_clusters_from_bpck(self, bpck_clusters, leaves: set=None):
		for cluster in bpck_clusters:
			ids=[]
			l=0
			for e in cluster:
				ids.extend([ids for ids in e[1:]])
				l+=e[0]
			self.clusters.append(Cluster(ids=ids, length=l))
		if leaves: self.leaves.update(leaves)

	def join_clusters(self):
		# Join clusters in the group by binid
		final_clusters = {}
		for c in self.clusters:
			if c.binid not in final_clusters: final_clusters[c.binid] = Cluster(binid=c.binid)
			final_clusters[c.binid].update(c)
		self.clusters = list(final_clusters.values())

	def merge(self, group):
		self.add_clusters(group.get_leaves(), group.get_clusters())
		
	def get_length(self):
		return sum([c.get_length() for c in self.clusters])

	def get_cluster_count(self):
		return len(self.clusters)

	def __repr__(self):
		args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
		return 'Group({})'.format(', '.join(args))