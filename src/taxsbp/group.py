from taxsbp.cluster import cluster


class group:
    def __init__(self):
        self.leaves = set()
        self.clusters = []

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
        # Return list of tuples with the clusters in the format necessary for the binpacking
        # Format [(length,seq1,seq2,...,seqN),...]
        # Example: [(500,A,B,C),(300,D),(200,E)]
        return [c.get_tuples() for c in self.clusters]

    def add_clusters_from_bpck(self, bpck_clusters, leaves: set = ()):
        # Parse binpacking output - list of lists with tuples generated with get_clusters_to_bpck
        # Example: [[(500,A,B,C)],[(300,D),(200,E)]]

        # For each cluster returned by binpaking
        for c in bpck_clusters:
            if not c:
                continue
            # split clusters in their respective binid assigned (or None)
            slen = 0
            ids = []
            for e in c:
                slen += e[0]
                ids.extend(e[1:])
            self.clusters.append(cluster(ids=ids, length=slen))
        if leaves:
            self.leaves.update(leaves)

    def join_clusters(self, ids: list | None = None):
        # Join all clusters inside the group, do not join clusters wiht different binids
        joined_clusters = cluster()
        single_clusters = []
        for c in self.clusters:
            if ids:
                if c.ids in ids:
                    joined_clusters.update(c)
                else:
                    single_clusters.append(c)
            else:
                joined_clusters.update(c)

        if joined_clusters.ids:
            self.clusters = single_clusters + [joined_clusters]
        else:
            self.clusters = single_clusters

    def merge(self, group):
        self.add_clusters(group.get_leaves(), group.get_clusters())

    def get_length(self):
        return sum([c.get_length() for c in self.clusters])

    def get_cluster_count(self):
        return len(self.clusters)

    def __repr__(self):
        args = [f"{k}={v!r}" for (k, v) in vars(self).items()]
        return "Group({})".format(", ".join(args))
