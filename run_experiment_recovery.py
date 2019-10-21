from models import BayesianNetwork, SubplotDrawer
from compilation import compile_variable_elimination, EliminationOrdering
from convert_ac2spn import convert_ac2spn
from comp_marg_spn import compile_marginalized_spn
from decompilation import decompile
import operator as op
from functools import reduce
import os
import argparse


def ncr(n, r):
    # Code from:
    # https://stackoverflow.com/questions/
    # 4941753/is-there-a-math-ncr-function-in-python
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom


def are_bns_same(ori_bn, decomp_bn):
    amt_parents_v_struct = [len(list(ori_bn.dag.predecessors(n)))
                            for n in ori_bn.dag.nodes()
                            if len(list(ori_bn.dag.predecessors(n))) > 1]
    extra_edges = sum(
        list(
            map(lambda n: ncr(n, 2), amt_parents_v_struct)
        )
    )
    ori_amt_edges = len(ori_bn.dag.edges)
    decomp_amt_edges = len(decomp_bn.dag.edges)
    max_amt_edges = len(ori_bn.dag.edges) + extra_edges

    has_expec_amt_edges = decomp_amt_edges >= ori_amt_edges and\
        decomp_amt_edges <= max_amt_edges
    has_eq_amt_nodes = len(
        ori_bn.dag.nodes()) == len(decomp_bn.dag.nodes())

    return has_eq_amt_nodes and has_expec_amt_edges


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Experiment for recoverying BN")
    parser.add_argument(
        "--bn", type=str, help="Path to '.bn' file or folder name.",
        required=True)
    parser.add_argument(
        "--elim-ord", type=str, help="Elimination ordering type.",
        choices=["rev", "top", "mn", "mw", "all_rev"], required=True)
    parser.add_argument(
        "--plot-all",
        help="Plot all graphs.",
        action="store_true")
    args = parser.parse_args()

    # Parameters
    bn_path = args.bn
    is_collapse_sums = False
    elim_ord_type = args.elim_ord
    comp_marg_spn = "internal"
    plot_diff = True
    plot_all = args.plot_all

    bn_files = []
    is_file = bn_path.endswith(".bn")
    if is_file:
        bn_files = [bn_path]
    else:
        bn_files = [os.path.join("bns", name)
                    for name in os.listdir(bn_path)]

    diff_decomp = []
    for bn_file_name in bn_files:
        # plotting
        bn_name = bn_file_name.replace(".bn", "")

        bn = BayesianNetwork.get_bn_from_file(bn_file_name)

        all_elim_ord = EliminationOrdering.get_elimination_ordering(
            bn, elim_ord_type)
        all_elim_ord = all_elim_ord if type(all_elim_ord[0]) is list\
            else [all_elim_ord]
        for i, elim_ord in enumerate(all_elim_ord):
            plotter = SubplotDrawer("Compiling and Decompiling BN - " +
                                    bn_file_name)
            plotter.add(bn, "BN")

            # --- Compile and Decompile
            # 1) Compilation
            ac = compile_variable_elimination(bn, elim_ord)
            plotter.add(ac, "Compiled AC - " + str(elim_ord))

            # 2) Convert to SPN
            spn = convert_ac2spn(ac, subplot_drawer=plotter)

            # 3) Marginalized SPN
            marg_nodes = []
            if comp_marg_spn == "internal":
                marg_nodes = [var for var in bn.dag.nodes()
                              if len(list(bn.dag.successors(var))) > 0]
            spn = compile_marginalized_spn(spn, marg_nodes,
                                           collapse_sums=is_collapse_sums,
                                           subplot_drawer=None)
            comp_marg_subplot_title = "Compiled Marginalized SPN"
            if is_collapse_sums:
                comp_marg_subplot_title += " - Collpase Sum"
            if comp_marg_spn == "none":
                comp_marg_subplot_title += " - None Marg"
            elif comp_marg_spn == "internal":
                comp_marg_subplot_title += " - Marg Internal"
            plotter.add(spn, comp_marg_subplot_title)

            # 4) Decompilation
            decomp_bn = decompile(spn, comp_assumption="ve")
            plotter.add(decomp_bn, "Decompiled BN")

            # --- Compare
            print("------ " + bn_name + " Report - Elim Ord: "
                  + str(i + 1) + "/" + str(len(all_elim_ord)) + " ------")
            print("-> Ori BN #nodes: " + str(len(bn.dag.nodes())))
            print("-> Decomp BN #nodes: " + str(len(decomp_bn.dag.nodes())))
            print("-> Ori BN #edges: " + str(len(bn.dag.edges())))
            print("-> Decomp BN #edges: " + str(len(decomp_bn.dag.edges())))
            bns_same = are_bns_same(bn, decomp_bn)
            if not bns_same:
                diff_decomp.append(bn_name)
            if (plot_diff and not bns_same) or plot_all:
                plotter.plot()

    print("--- Report ---")
    print("-> Different decompilations: " + str(len(diff_decomp))
          + "/" + str(len(bn_files)))
    print("-> Diff decomp BNs:" + str(diff_decomp))
