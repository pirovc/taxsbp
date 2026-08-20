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
    # groups, info = parse_input(
    #     input_file,
    #     tax,
    #     pre_cluster,
    #     pre_cluster_file,
    #     bin_exclusive,
    #     bin_exclusive_file,
    #     bin_exclusive_prefix,
    # )

    groups, info = parse_input3(
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


def parse_input3(
    filename,
    tax,
    pre_cluster_rank,
    pre_cluster_file,
    bin_exclusive_rank,
    bin_exclusive_file,
    bin_exclusive_prefix,
):
    # Parse input file into info dict
    info = {}
    with open(filename, "r") as infile:
        for line in infile:
            uid, weigth, node = line.rstrip().split("\t")
            latest_node = tax.latest(node)
            if latest_node and latest_node != node:
                print(f"{node} -> {latest_node}", file=sys.stderr)
            info[uid] = [int(weigth), latest_node]

    # Create pre-clusters
    # Store pre-clusters into tuple since there can be more than
    # one pre-cluster based on the same node due to the lca
    pre_clusters = []
    pre_cluster_node = {}
    if pre_cluster_file:
        # get LCA of the ids on the pre cluster file
        with open(pre_cluster_file, "r") as pinfile:
            for line in pinfile:
                clusters = set(line.rstrip().split("\t"))
                lca_pc = tax.lca([info[uid][1] for uid in clusters])
                pre_clusters.append((lca_pc, clusters))
                for uid in clusters:
                    pre_cluster_node[uid] = lca_pc

    # Map bin exclusive to unique groups (incremental int)
    bin_exclusive_groups = []
    bin_exclusive_groups_idx = {}
    if bin_exclusive_file:
        with open(bin_exclusive_file, "r") as binfile:
            for line in binfile:
                be_groups = set(line.rstrip().split("\t"))
                bin_exclusive_groups.append(be_groups)
                for uid in be_groups:
                    bin_exclusive_groups_idx[uid] = len(bin_exclusive_groups)

    if bin_exclusive_rank or pre_cluster_rank:
        pre_clusters_aux = {}
        bin_exclusive_groups_aux = {}
        for uid, (_, latest_node) in info.items():
            # bin exclusive node
            be_node = None
            if bin_exclusive_rank == "leaves":
                be_node = latest_node
            elif bin_exclusive_rank:
                parent_node = tax.parent_rank(latest_node, rank=bin_exclusive_rank)
                if parent_node != tax.undefined_node:
                    be_node = parent_node
            if be_node:
                if be_node not in bin_exclusive_groups_aux:
                    bin_exclusive_groups_aux[be_node] = set()
                bin_exclusive_groups_aux[be_node].add(uid)
                bin_exclusive_groups_idx[uid] = len(bin_exclusive_groups_aux)

            # pre cluster node
            pc_node = None
            if pre_cluster_rank == "leaves":
                pc_node = latest_node
            elif pre_cluster_rank:
                parent_node = tax.parent_rank(latest_node, rank=pre_cluster_rank)
                if parent_node != tax.undefined_node:
                    pc_node = parent_node
            if pc_node:
                if pc_node not in pre_clusters_aux:
                    pre_clusters_aux[pc_node] = set()
                pre_clusters_aux[pc_node].add(uid)
                pre_cluster_node[uid] = pc_node
        pre_clusters = list(pre_clusters_aux.items())
        bin_exclusive_groups = list(bin_exclusive_groups_aux.values())

    # Check bin exclusive overlap pre-cluster
    if pre_clusters and bin_exclusive_groups:
        for i, be_group in enumerate(bin_exclusive_groups):
            for _, pre_cluster in pre_clusters:
                # Reject if bin exclusive group partially overlaps with a pre-cluster
                intersect = pre_cluster.intersection(be_group)
                if intersect and len(intersect) < len(pre_cluster):
                    print(
                        f"{','.join(be_group)} cannot be bin exclusive due to a partial overlap with a pre-cluster",
                        file=sys.stderr,
                    )
                    for uid in be_group:
                        del bin_exclusive_groups_idx[uid]
                    del bin_exclusive_groups[i]

    groups = {}
    # Add pre clusters before
    for node, uids in pre_clusters:
        if node not in groups:
            groups[node] = group()
        groups[node].add_clusters([cluster(uids, sum(info[uid][0] for uid in uids))])

    print(pre_clusters)
    print(bin_exclusive_groups)
    print(bin_exclusive_groups_idx)

    for uid, (weigth, latest_node) in info.items():
        pc_node = pre_cluster_node.get(uid, None)
        be_idx = bin_exclusive_groups_idx.get(uid, None)

        print(uid, pc_node, be_idx)

        if be_idx:
            # Create bin exclusive node with prefix
            node = f"{bin_exclusive_prefix}{be_idx}"

            if node not in groups:
                groups[node] = group()

            # If node was previously pre-clustered, move specific cluster to bin exclusive group
            if pc_node:
                # May not be in groups since it was already moved from a previous entry
                if pc_node in groups:
                    # May not be in clusters, already moved from previous pre-clustered uid
                    c = groups[pc_node].get_clusters(with_id=uid)
                    if c:
                        groups[node].add_clusters([c])
                        groups[pc_node].clear_clusters(with_id=uid)
                        # Delete entry if was only clusters
                        if groups[pc_node].get_cluster_count() == 0:
                            del groups[pc_node]
            else:
                groups[node].add_clusters([cluster([uid], weigth)])

        elif not pc_node:
            if latest_node not in groups:
                groups[latest_node] = group()
            groups[latest_node].add_clusters([cluster([uid], weigth)])

    print(groups)
    return groups, info


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
