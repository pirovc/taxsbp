class cluster:
    def __init__(self, ids: list | None = None, length: int = 0):
        self.length = length
        self.ids = set(ids) if ids is not None else set()

    def get_tuples(self):
        return (self.get_length(),) + tuple(self.get_ids())

    def get_length(self):
        return self.length

    def get_ids(self):
        return self.ids

    def __repr__(self):  # pragma: no cover
        args = [f"{k}={v!r}" for (k, v) in vars(self).items()]
        return "Cluster({})".format(", ".join(args))
