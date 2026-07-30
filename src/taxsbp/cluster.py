class cluster:
    def __init__(self, ids: list | None = None, length: int = 0):
        self.length = length
        self.ids = set(ids) if ids is not None else set()
        self.binid = None

    def set_binid(self, binid):
        self.binid = binid

    def get_binid(self):
        return self.binid

    def get_tuples(self):
        return (self.get_length(),) + tuple(self.get_ids())

    def get_length(self):
        return self.length

    def get_ids(self):
        return self.ids

    def add(self, uid, length):
        self.ids.add(uid)
        self.length += length

    def update(self, cluster):
        self.ids.update(cluster.get_ids())
        self.length += cluster.get_length()

    def __repr__(self):
        args = [f"{k}={v!r}" for (k, v) in vars(self).items()]
        return "Cluster({})".format(", ".join(args))
