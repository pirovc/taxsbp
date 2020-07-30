from taxsbp.Cluster import Cluster

class Group:

	def __init__(self):
		self.leaves = set()
		self.clusters = []

	def add_cluster(self,leaf,uid,length,binid: int=None):
		self.leaves.add(leaf)
		self.clusters.append(Cluster([uid],length,binid))

	def add_clusters(self, leaves, clusters):
		self.leaves.update(leaves)
		self.clusters.extend(clusters)

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
		# For each cluster returned by binpaking
		for cluster in bpck_clusters:
			# split clusters in their respective binid assigned (or None)
			cluster_binids = {}

			for e in cluster:
				if e[0] not in cluster_binids: cluster_binids[e[0]] = [0]
				cluster_binids[e[0]][0]+=e[1]
				cluster_binids[e[0]].extend(e[2:])

			# in case of multiple bins and new sequences added to the same cluster, re-assign to just one
			if len(cluster_binids)>=2 and None in cluster_binids:
				none_cluster = cluster_binids.pop(None)
				smallest_binid = min(cluster_binids)
				cluster_binids[smallest_binid][0]+=none_cluster[0]
				cluster_binids[smallest_binid].extend(none_cluster[1:])

			# Add clusters to group
			for bid, cl in cluster_binids.items():
				self.clusters.append(Cluster(ids=cl[1:], length=cl[0], binid=bid))
		if leaves: self.leaves.update(leaves)

	def join_clusters(self):
		# Join all clusters inside the group, do not join clusters wiht different binids
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