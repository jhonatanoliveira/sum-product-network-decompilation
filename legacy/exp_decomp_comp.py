import argparse
import os
from main import get_bn_from_file, reconstruct_bn,\
    compile_variable_elimination, draw_subplot_graphs

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compile a BN into an AC,\
        decompile the AC into a BN, then compile this BN into an AC")
    parser.add_argument("--bn", type=str, help="Path to '.bn' file.")
    parser.add_argument("--elim-ord", type=str,
                        help="Elimination ordering to be used.", default="rev")
    parser.add_argument(
        "--plot-diff",
        help="Plot only the graphs with differences in #nodes, edges.",
        action="store_true")
    args = parser.parse_args()

    bn_files = []
    is_file = args.bn.endswith(".bn")
    if is_file:
        bn_files = [args.bn]
    else:
        bn_files = [os.path.join("bns", name)
                    for name in os.listdir(args.bn)]

    for bn_file_name in bn_files:
        if bn_file_name.endswith(".bn"):

            bn_file = bn_file_name if not is_file else bn_file_name
            graphs = []
            graphs_subtitles = []

            ori_bn = get_bn_from_file(bn_file)
            graphs.append(ori_bn)
            graphs_subtitles.append("Original BN - N,E: (" + str(
                len(ori_bn.dag.nodes())) + "," + str(
                    len(ori_bn.dag.edges())) + ")")
            elim_ord = ori_bn.get_elimination_ordering(args.elim_ord)
            print("Elimination Ordering: " + str(elim_ord))
            ori_ac = compile_variable_elimination(ori_bn, elim_ord)
            graphs.append(ori_ac)
            graphs_subtitles.append("Original AC - N,E: (" + str(
                len(ori_ac.dag.nodes())) + "," + str(
                    len(ori_ac.dag.edges())) + ")")
            bn1 = reconstruct_bn(ori_bn, elim_ord, plot=False)
            graphs.append(bn1)
            graphs_subtitles.append("BN 1 - N,E: (" + str(
                len(bn1.dag.nodes())) + "," + str(
                    len(bn1.dag.edges())) + ")")

            elim_ord = bn1.get_elimination_ordering(args.elim_ord)
            print("Elimination Ordering: " + str(elim_ord))
            ac1 = compile_variable_elimination(bn1, elim_ord)
            graphs.append(ac1)
            graphs_subtitles.append("AC 1 - N,E: (" + str(
                len(ac1.dag.nodes())) + "," + str(
                    len(ac1.dag.edges())) + ")")
            bn2 = reconstruct_bn(bn1, elim_ord, plot=False)
            graphs.append(bn2)
            graphs_subtitles.append("BN 2 - N,E: (" + str(
                len(bn2.dag.nodes())) + "," + str(
                    len(bn2.dag.edges())) + ")")

            if args.plot_diff:
                size_bn1 = (len(bn1.dag.nodes()), len(bn1.dag.edges()))
                size_bn2 = (len(bn2.dag.nodes()), len(bn2.dag.edges()))
                if size_bn1 != size_bn2:
                    main_title = bn_file + " - Decompiling a decompiled BN"
                    draw_subplot_graphs(graphs, graphs_subtitles, main_title)
                else:
                    print("No difference found between BNs.")
            else:
                main_title = bn_file + " - Decompiling a decompiled BN"
                draw_subplot_graphs(graphs, graphs_subtitles, main_title)
