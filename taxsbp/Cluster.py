class Cluster:
	def __init__(self, seqid=None, seqlen=None, binid: int=None):
		self.seqlen = [(seqlen,seqid)] if seqid and seqlen else []
		self.binid = binid

	def set_binid(self, binid):
		self.binid = binid

	def get_tuples(self):
		return self.seqlen

	def get_cluster_length(self):
		return sum(i[0] for i in self.seqlen)

	def get_ids(self):
		return [i[1] for i in self.seqlen]

	def add(self, seqid, seqlen):
		self.seqlen.append((seqlen,seqid))

	def update(self,cluster):
		self.seqlen.extend(cluster.seqlen)

	def __repr__(self):
		args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
		return 'Cluster({})'.format(', '.join(args))