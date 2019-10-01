import networkx as nx
import argparse
import os
from main import get_bn_from_file, reconstruct_bn, draw_subplot_graphs


def plot_graphs(name, graphs, graphs_subtitles,
                plot_diff=True, max_plot_rec=-1):
    graphs = [g for g in graphs]
    graphs_subtitles = [g for g in graphs_subtitles]

    main_title = name + ", reconstructions: " + str(len(graphs) - 1)

    if plot_diff:
        diff_countings = []  # [(#nodes, #edges)]
        for rec_bn in [g for g in graphs[1:]]:
            counting = (len(rec_bn.dag.nodes()), len(rec_bn.dag.edges()))
            if counting not in diff_countings:
                diff_countings.append(counting)
            else:
                idx = graphs.index(rec_bn)
                graphs.pop(idx)
                graphs_subtitles.pop(idx)
        main_title += ", Diff Countings: " + str(len(diff_countings))

    if max_plot_rec > 0:
        graphs = graphs[:(max_plot_rec + 1)]
        graphs_subtitles = graphs_subtitles[:(max_plot_rec + 1)]
        main_title += ", Max Plot: " + str(max_plot_rec)

    if plot_diff and len(graphs) < 3:
        print("No different reconstructions")
    else:
        draw_subplot_graphs(graphs, graphs_subtitles, main_title)


def all_top_ord(bn_file):

    bn = get_bn_from_file(bn_file)
    all_elim_ord = list(nx.algorithms.dag.all_topological_sorts(bn.dag))
    total_elim_ord = len(all_elim_ord)
    print("Total elimination oderings: " + str(total_elim_ord))

    graphs = []
    graphs_subtitles = []

    graphs.append(bn)
    graphs_subtitles.append("Original BN")
    for idx, elim_ord in enumerate(all_elim_ord):
        elim_ord = list(reversed(elim_ord))
        rec_bn = reconstruct_bn(bn, elim_ord, plot=False)
        graphs.append(rec_bn)
        graphs_subtitles.append(
            "#" + str(idx) + " - " + ",".join(elim_ord))

    return (graphs, graphs_subtitles)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compile a BN into an AC")
    parser.add_argument("--bn", type=str, help="Path to '.bn' file.")
    parser.add_argument(
        "--max-plot-rec", type=int,
        help="Max amount of plotting reconstructions.",
        default=-1)
    parser.add_argument(
        "--plot-diff",
        help="Force to print all reconstructions.",
        action="store_true")
    args = parser.parse_args()

    bn_files = []
    is_file = args.bn.endswith(".bn")
    if is_file:
        bn_files = [args.bn]
    else:
        bn_files = os.listdir(args.bn)

    for bn_file_name in bn_files:
        if bn_file_name.endswith(".bn"):
            bn_file = os.path.join("bns", bn_file_name) if not is_file \
                else bn_file_name
            graphs, graphs_subtitles = all_top_ord(bn_file)
            name = bn_file[:-3]
            plot_graphs(
                name,
                graphs,
                graphs_subtitles,
                plot_diff=args.plot_diff,
                max_plot_rec=args.max_plot_rec)
