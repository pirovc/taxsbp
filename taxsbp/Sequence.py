from taxsbp.Cluster import Cluster

class Sequence:

	def __init__(self,seqlen,taxid,specialization=None,binid=None):
		self.seqlen = seqlen
		self.taxid = taxid
		self.specialization = specialization
		self.binid = binid

	def __repr__(self):
		args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
		return 'Sequence({})'.format(', '.join(args))