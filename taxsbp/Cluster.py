class Cluster:
	def __init__(self, seqid=None, seqlen=None):
		self.seqlen = {seqid: seqlen} if seqid and seqlen else {}

	def update_cluster(self,cluster):
		self.seqlen.update(cluster.get_seqlen())

	def get_seqlen(self):
		return self.seqlen

	def get_length(self):
		return sum(self.seqlen.values())

	def get_ids(self):
		return self.seqlen.keys()

	def __repr__(self):
		args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
		return 'Cluster({})'.format(', '.join(args))