import argparse
import sys

from binpacking.numpy import to_constant_volume
from multitax import CustomTx

from taxsbp import __version__
from taxsbp.cluster import cluster
from taxsbp.group import group


def main(arguments: str | None = None):

    if arguments is not None:
        sys.argv = arguments

    parser = argparse.ArgumentParser(
        prog="taxsbp", conflict_handler="resolve", add_help=True
    )
    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        help="Input file: unique_id <tab> weigth <tab> node",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        help="Output file: id <tab> weigth <tab> node <tab> bin. Default: STDOUT",
    )
    parser.add_argument(
        "-t", "--taxonomy-file", help="Taxonomy file: node <tab> parent <tab> rank"
    )
    parser.add_argument(
        "-l",
        "--bin-len",
        type=int,
        help="Maximum bin length (in bp). Use this parameter insted of -b to define the number of bins. Default: length of the biggest group [Mutually exclusive -b]",
    )
    parser.add_argument(
        "-b",
        "--n-bins",
        type=int,
        help="Approximate number of bins (estimated by total length/bin number). [Mutually exclusive -l]",
    )
    parser.add_argument(
        "-p",
        "--pre-cluster",
        type=str,
        default="",
        help="",
    )
    parser.add_argument(
        "--pre-cluster-file",
        type=str,
        default="has precedence over --pre-cluster",
        help="",
    )
    parser.add_argument(
        "-e",
        "--bin-exclusive",
        type=str,
        default="",
        help="",
    )
    parser.add_argument(
        "--bin-exclusive-file",
        type=str,
        default="",
        help="",
    )
    parser.add_argument(
        "-s",
        "--silent",
        default=False,
        action="store_true",
        help="Ignore warning",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="version: %(prog)s " + __version__,
        help="Show program's version number and exit.",
    )

    if len(sys.argv) == 1:  # Print help calling script without parameters
        parser.print_help()
        return False

    args = parser.parse_args()  # read sys.argv[1:] by default

    tax = CustomTx(files=args.taxonomy_file)
    groups, lens = parse_input(
        args.input_file, tax, args.pre_cluster, args.pre_cluster_file
    )

    # Keep only used nodes on tax
    tax.filter(groups.keys())

    # Define bin length
    if args.bin_len:  # user defined
        blen = args.bin_len
    elif args.n_bins:
        blen = sum([g.get_length() for g in groups.values()]) / float(args.n_bins)
    else:  # Default bin length on the max group length
        blen = max([g.get_length() for g in groups.values()])

    clusterx(groups, tax, blen, args.bin_exclusive)

    set_bins(groups)
    print_stats(groups)
    res = generate_results(groups, lens)
    if args.output_file:
        with open(args.output_file, "w") as file:
            for r in res:
                print(*r, sep="\t", file=file)
    else:
        for r in res:
            print(*r, sep="\t", file=sys.stdout)


def clusterx(groups, tax, blen, bin_exclusive_rank):
    if bin_exclusive_rank:
        rank_taxids, orphan_taxids = get_rank_taxids(bin_exclusive_rank, groups, tax)
        if rank_taxids:
            # clustering directly on the rank chosen, recursion required for children nodes
            for rank_taxid in rank_taxids:
                ApproxSBP(rank_taxid, None, groups, tax, blen)
        if orphan_taxids:
            # clustering directly on the taxid level, no recursion to children nodes necessary
            for orphan_taxid in orphan_taxids:
                bpck(groups, orphan_taxid, orphan_taxid, blen)
    else:  # default mode
        ApproxSBP(tax.root_node, None, groups, tax, blen)


def get_rank_taxids(bin_exclusive_rank, groups, tax):
    rank_taxids = set()
    orphan_taxids = set()
    # if not working on leaf level
    if bin_exclusive_rank != "leaves":
        for leaf in groups:
            t = tax.parent_rank(leaf, rank=bin_exclusive_rank)
            if t == tax.undefined_node:
                orphan_taxids.add(leaf)
            else:
                rank_taxids.add(t)
    else:
        orphan_taxids = set(groups.keys())
    return rank_taxids, orphan_taxids


def bpck(groups, node, parent, blen):
    # Perform bin packing on a single node
    # it packs the clusters on groups[node] and add to groups[parent]
    # if node and parent are equal, root was reached
    if node not in groups:
        return

    at_root = node == parent

    # If there is only one cluster, do not need to pack
    if groups[node].get_cluster_count() == 1:
        if not at_root:  # transfer cluster to parent if not root
            if parent not in groups:
                groups[parent] = group()
            groups[parent].merge(groups[node])
            del groups[node]
    else:
        # Perform bin packing
        clusters = to_constant_volume(
            groups[node].get_clusters_to_bpck(), blen, weight_pos=0
        )

        if clusters:
            if parent not in groups:
                groups[parent] = group()
            if not at_root:
                # Parse clustered results into parent node and remove actual node
                groups[parent].add_clusters_from_bpck(
                    clusters, leaves=groups[node].get_leaves()
                )
                del groups[node]
            else:  # if root
                # Parse clustered results into same node (clear it before)
                groups[parent].clear_clusters()
                groups[parent].add_clusters_from_bpck(clusters)


def ApproxSBP(node, parent, groups, tax, blen):
    # Function to perform hiearchical bin packing recursively
    # Recursively call to pack sorted list of children (to get always same results)
    for child in sorted(tax.children(node), key=str):
        ApproxSBP(child, node, groups, tax, blen)

    # If node is a leaf - no child in children[node]
    # or
    # After all children of a node were packed in the for loop, pack node itself into parent
    bpck(groups, node, parent if parent is not None else node, blen)


def set_bins(groups):
    binid_count = -1
    for g in groups.values():
        for c in g.get_clusters():
            binid_count += 1
            c.set_binid(binid_count)


def generate_results(groups, lens):
    for g in groups.values():
        for c in g.get_clusters():
            for seqid in c.get_ids():
                yield [seqid, lens[seqid], str(c.get_binid())]


def print_stats(groups):

    c_lens = []
    c_ids = []
    cnt = 0
    for g in groups.values():
        for c in g.get_clusters():
            c_ids.append(len(c.ids))
            c_lens.append(c.length)
            cnt += 1

    print(f"Min. cluster w: {min(c_lens)}", file=sys.stderr)
    print(f"Avg. cluster w: {sum(c_lens) / cnt}", file=sys.stderr)
    print(f"Max. cluster w: {max(c_lens)}", file=sys.stderr)

    print(f"Min. ids/cluster: {min(c_ids)}", file=sys.stderr)
    print(f"Avg. ids/cluster: {sum(c_ids) / cnt}", file=sys.stderr)
    print(f"Max. ids/cluster: {max(c_ids)}", file=sys.stderr)


def parse_input(filename, tax, pre_cluster_rank, pre_cluster_file):
    groups = {}
    lens = {}
    pre_clustered_nodes = set()
    pre_clustered_ids = {}
    if pre_cluster_file:
        uid_node = {}
        with open(filename, "r") as infile:
            for line in infile:
                uid, _, node = line.rstrip().split("\t")
                unode = tax.latest(node)
                if unode and unode != node:
                    print(f"{node} -> {unode}", file=sys.stderr)
                uid_node[uid] = unode

        with open(pre_cluster_file, "r") as pinfile:
            for line in pinfile:
                clusters = line.rstrip().split("\t")
                lca = tax.lca([uid_node[uid] for uid in clusters])
                pre_clustered_ids.update({uid: lca for uid in clusters})
                pre_clustered_nodes.add(lca)

    with open(filename, "r") as infile:
        for line in infile:
            uid, w, node = line.rstrip().split("\t")

            if pre_cluster_file and uid in pre_clustered_ids:
                unode = pre_clustered_ids[uid]
            else:
                unode = tax.latest(node)
                if unode and unode != node:
                    print(f"{node} -> {unode}", file=sys.stderr)

            # Pre-cluster by rank
            if pre_cluster_rank:
                if pre_cluster_rank != "leaves":
                    pnode = tax.parent_rank(unode, rank=pre_cluster_rank)
                    if pnode != tax.undefined_node:
                        unode = pnode
                pre_clustered_nodes.add(unode)

            if unode != tax.undefined_node:
                if unode not in groups:
                    groups[unode] = group()

                groups[unode].add_clusters([unode], [cluster([uid], int(w))])
                lens[uid] = int(w)
            else:
                print(node + " not found", file=sys.stderr)

    print(groups)
    # Join pre-clustered entries
    # Only join selected unique ids if pre_cluster_file
    if pre_cluster_rank or pre_cluster_file:
        for node, g in groups.items():
            if node in pre_clustered_nodes:
                g.join_clusters(
                    ids=[{i} for i in pre_clustered_ids] if pre_clustered_ids else []
                )
    print(groups)
    return groups, lens


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
