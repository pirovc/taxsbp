from taxsbp.cluster import cluster


class group:
    def __init__(self):
        self.clusters = []

    def add_clusters(self, clusters):
        self.clusters.extend(clusters)

    def get_clusters(self, with_id=""):
        if with_id:
            for c in self.clusters:
                if with_id in c.ids:
                    return c
        else:
            return self.clusters

    def clear_clusters(self, with_id=""):
        if with_id:
            for i, c in enumerate(self.clusters):
                if with_id in c.ids:
                    del self.clusters[i]
                    break
        else:
            self.clusters = []

    def get_clusters_to_bpck(self):
        # Return list of tuples with the clusters in the format necessary for the binpacking
        # Format [(length,seq1,seq2,...,seqN),...]
        # Example: [(500,A,B,C),(300,D),(200,E)]
        return [c.get_tuples() for c in self.clusters]

    def add_clusters_from_bpck(self, bpck_clusters):
        # Parse binpacking output - list of lists with tuples generated with get_clusters_to_bpck
        # Example: [[(500,A,B,C)],[(300,D),(200,E)]]

        # For each cluster returned by binpaking
        for c in bpck_clusters:
            if c:
                slen = 0
                ids = []
                for e in c:
                    slen += e[0]
                    ids.extend(e[1:])
                self.clusters.append(cluster(ids=ids, length=slen))

    def merge(self, group):
        self.add_clusters(group.get_clusters())

    def get_length(self):
        return sum([c.get_length() for c in self.clusters])

    def get_cluster_count(self):
        return len(self.clusters)

    def __repr__(self):  # pragma: no cover
        args = [f"{k}={v!r}" for (k, v) in vars(self).items()]
        return "Group({})".format(", ".join(args))
