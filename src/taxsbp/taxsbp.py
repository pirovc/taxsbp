import argparse
import sys
from statistics import mean, median, stdev

from binpacking.numpy import to_constant_volume
from multitax import CustomTx

from taxsbp import __version__
from taxsbp.cluster import cluster
from taxsbp.group import group


def main():
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
        "-t",
        "--taxonomy-file",
        required=True,
        help="Taxonomy file: node <tab> parent <tab> rank",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        help="Output file: id <tab> weigth <tab> target node <tab> binno. Default: STDOUT",
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
        default="",
        help="has precedence over --pre-cluster",
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
        "--stats",
        default=False,
        action="store_true",
        help="Output stats to stderr",
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

    args = parser.parse_args()
    opt = vars(args)

    bins, sts = taxsbp(**opt)

    # Print bins to STDOUT if called cli without --output-file
    if not args.output_file:
        for b in bins:
            print(*b, sep="\t")

    # Print stats to STDOUT if called cli with --stats
    if args.stats:
        import pprint

        pprint.pprint(sts, compact=True, sort_dicts=False, stream=sys.stderr)

    sys.exit(0)


def taxsbp(
    input_file: str,
    taxonomy_file: str,
    output_file: str | None = None,
    bin_len: int | None = None,
    n_bins: int | None = None,
    pre_cluster: str | None = None,
    pre_cluster_file: str | None = None,
    bin_exclusive: str | None = None,
    bin_exclusive_file: str | None = None,
    stats: bool = False,
):
    bin_exclusive_prefix = "@@be@@-"
    tax = CustomTx(files=taxonomy_file)
    groups, info = parse_input(
        input_file,
        tax,
        pre_cluster,
        pre_cluster_file,
        bin_exclusive,
        bin_exclusive_file,
        bin_exclusive_prefix,
    )

    # Keep only used nodes on taxonomy for faster
    tax.filter([n for _, n in info.values()])

    # Define bin length
    if bin_len:  # user defined
        blen = bin_len
    elif n_bins:
        blen = int(sum([g.get_length() for g in groups.values()]) / float(n_bins))
    else:  # Default bin length on the max group length
        blen = max([g.get_length() for g in groups.values()])

    # Cluster bin exclusive groups indivudually
    need_recursive = False
    if bin_exclusive or bin_exclusive_file:
        for node in groups:
            if node.startswith(bin_exclusive_prefix):
                bpck(groups, node, node, blen)
            else:
                need_recursive = True
    else:
        need_recursive = True

    if need_recursive:
        ApproxSBP(tax.root_node, None, groups, tax, blen)

    if output_file:
        with open(output_file, "w") as outf:
            for b in generate_bins(groups, info):
                print(*b, sep="\t", file=outf)

    sts = generate_stats(generate_bins(groups, info), info, blen) if stats else None

    return list(generate_bins(groups, info)), sts


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
                groups[parent].add_clusters_from_bpck(clusters)
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


def generate_bins(groups, info):
    binid = 0
    for g in groups.values():
        for c in g.get_clusters():
            for uid in c.get_ids():
                yield [uid, info[uid][0], info[uid][1], binid]
            binid += 1


def generate_stats(bins, info, blen):
    cnts = {}
    for uid, weigth, _, binid in bins:
        if binid not in cnts:
            cnts[binid] = {"nodes": set(), "ids": set(), "weight": 0}
        cnts[binid]["nodes"].add(info[uid][1])
        cnts[binid]["ids"].add(uid)
        cnts[binid]["weight"] += weigth

    sts = {}

    def general_stats(vals):
        return {
            "min": min(vals),
            "max": max(vals),
            "avg": mean(vals),
            "med": median(vals),
            "sd": stdev(vals) if len(vals) > 1 else None,
        }

    sts["target_weigth"] = blen
    sts["total_bins"] = binid + 1
    sts["weigths"] = general_stats([c["weight"] for c in cnts.values()])
    sts["ids"] = general_stats([len(c["ids"]) for c in cnts.values()])
    sts["nodes"] = general_stats([len(c["nodes"]) for c in cnts.values()])

    return sts


def parse_input(
    filename,
    tax,
    pre_cluster_rank,
    pre_cluster_file,
    bin_exclusive_rank,
    bin_exclusive_file,
    bin_exclusive_prefix,
):
    groups = {}
    info = {}
    pre_clustered_nodes = set()
    pre_clustered_ids = {}

    if pre_cluster_file:
        uid_node = {}

        # Map unique id to node from input file
        with open(filename, "r") as infile:
            for line in infile:
                uid, _, node = line.rstrip().split("\t")
                unode = tax.latest(node)
                if unode and unode != node:
                    print(f"{node} -> {unode}", file=sys.stderr)
                uid_node[uid] = unode

        # get LCA of the ids on the pre cluster file
        with open(pre_cluster_file, "r") as pinfile:
            for line in pinfile:
                clusters = line.rstrip().split("\t")
                lca = tax.lca([uid_node[uid] for uid in clusters])
                pre_clustered_ids.update({uid: lca for uid in clusters})
                pre_clustered_nodes.add(lca)

    # Mark bin exclusive entries with prefix to be separated from the taxonomy tree
    bin_exclusive_ids = {}
    if bin_exclusive_file:
        # get LCA of the ids on the pre cluster file
        with open(bin_exclusive_file, "r") as binfile:
            for c, line in enumerate(binfile):
                bins = line.rstrip().split("\t")
                bin_exclusive_ids.update(
                    {uid: f"{bin_exclusive_prefix}{c}" for uid in bins}
                )

    with open(filename, "r") as infile:
        for line in infile:
            uid, w, node = line.rstrip().split("\t")

            if bin_exclusive_file and uid in bin_exclusive_ids:
                unode = bin_exclusive_ids[
                    uid
                ]  # Add node as custom node with bin exclusive prefix
            elif pre_cluster_file and uid in pre_clustered_ids:
                unode = pre_clustered_ids[
                    uid
                ]  # get LCA node if pre_cluster_file was given
            else:
                unode = tax.latest(node)
                if unode and unode != node:
                    print(f"{node} -> {unode}", file=sys.stderr)

            # Get parent rank node
            if pre_cluster_rank:
                if pre_cluster_rank == "leaves":
                    pre_clustered_nodes.add(unode)
                else:
                    # Only add to pre-cluster if given rank is available
                    pnode = tax.parent_rank(unode, rank=pre_cluster_rank)
                    if pnode != tax.undefined_node:
                        unode = pnode
                        pre_clustered_nodes.add(unode)

            if bin_exclusive_rank:
                # Use the taxid of the rank for the entry
                # If entry does not have the rank keep as a loose node (to be clustered)
                if bin_exclusive_rank == "leaves":
                    unode = f"{bin_exclusive_prefix}{unode}"
                else:
                    pnode = tax.parent_rank(unode, rank=bin_exclusive_rank)
                    if pnode != tax.undefined_node:
                        unode = f"{bin_exclusive_prefix}{pnode}"

            if unode != tax.undefined_node:
                if unode not in groups:
                    groups[unode] = group()
                groups[unode].add_clusters([cluster([uid], int(w))])
                info[uid] = (int(w), unode.replace(bin_exclusive_prefix, "", 1))
            else:
                print(node + " not found", file=sys.stderr)

    # Join pre-clustered entries
    # Only join selected unique ids if pre_cluster_file
    if pre_cluster_rank or pre_cluster_file:
        for node, g in groups.items():
            if node in pre_clustered_nodes:
                g.join_clusters(
                    ids=[{i} for i in pre_clustered_ids] if pre_clustered_ids else []
                )

    return groups, info


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
