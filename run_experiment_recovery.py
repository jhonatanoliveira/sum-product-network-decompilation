from models import BayesianNetwork, Factor, SubplotDrawer
from compilation import compile_variable_elimination, EliminationOrdering
from convert_ac2spn import convert_ac2spn
from comp_marg_spn import compile_marginalized_spn
from decompilation import decompile
import os
import argparse
import networkx as nx
import numpy as np
import itertools


def connected_dags(N):
    assert N > 1
    prev_BNs = [np.zeros((1, 1), dtype=np.int64)] if N == 2 else\
        connected_dags(N-1)
    for B in prev_BNs:
        for k in range(1, 2**(N-1)):
            new_BN = np.zeros((N, N), dtype=np.int64)
            new_BN[0, 1:] = np.array(
                [int(b) for b in np.binary_repr(k, width=N-1)])
            new_BN[1:, 1:] = B
            yield new_BN


def generate_all_bns(num_vars, var_cardinalities=None):
    if var_cardinalities is None:
        var_cardinalities = {var: 2 for var in range(num_vars)}
    bns = []
    for adj_matrix in connected_dags(num_vars):
        dag = nx.convert_matrix.from_numpy_matrix(
            adj_matrix, create_using=nx.DiGraph)
        factors = Factor.construct_factors(dag, var_cardinalities)
        bn = BayesianNetwork(dag, var_cardinalities, factors)
        bns.append(bn)

    return bns


def count_ancestor_moral_graph_edges(dag):
    vstructures = [n for n in dag.nodes()
                   if len(list(dag.predecessors(n))) > 1]
    H = dag.to_undirected()
    new_edges = 0
    for vstructure in vstructures:
        ancestors = nx.algorithms.dag.ancestors(dag, vstructure)
        poss_new_edges = itertools.combinations(ancestors, r=2)
        for (v1, v2) in poss_new_edges:
            if (v1, v2) not in H.edges():
                H.add_edge(v1, v2)
                new_edges += 1

    return new_edges


def triangulate(dag, elim_ord):
    moral_graph = EliminationOrdering.moral_graph(dag)
    graph_cp = moral_graph.copy()
    for var in elim_ord:
        poss_edges = itertools.combinations(
            list(graph_cp.neighbors(var)), r=2)
        for edge in poss_edges:
            if edge not in graph_cp.edges():
                moral_graph.add_edge(edge[0], edge[1])
                graph_cp.add_edge(edge[0], edge[1])
        graph_cp.remove_node(var)

    return moral_graph


# def is_exaustive_isomorphic(graph1, graph2):
#     is_iso = True or len(graph1.nodes()) == len(graph2.nodes())
#     # Check leaf equality
#     leaves_g1 = [n for n in graph1.nodes()
#                  if graph1.successors(n) == 0]
#     leaves_g2 = [n for n in graph2.nodes()
#                  if graph2.successors(n) == 0]
#     is_iso = is_iso or len(leaves_g1) == len(leaves_g2)
#     is_iso = is_iso or all([l1 in leaves_g2 for l1 in leaves_g1])
#     # Check internal nodes equality
#     comp_graph = graph2.copy()
#     node_labels = [n for n in graph1.nodes()]
#     old_labels = [n for n in comp_graph.nodes()]
#     for new_labels in itertools.permutations(node_labels):
        


def is_decomp_bn_same(ori_bn, decomp_bn, elim_ord):
    triang_ori = triangulate(ori_bn.dag, elim_ord)
    undirected_decomp = decomp_bn.dag.to_undirected()

    return nx.is_isomorphic(triang_ori, undirected_decomp)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Experiment for recoverying BN")
    parser.add_argument(
        "--bn", type=str, help="Path to '.bn' file or folder name.",
        required=True)
    parser.add_argument(
        "--elim-ord", type=str, help="Elimination ordering type.",
        choices=["rev", "top", "mn", "mw", "all-rev"], required=True)
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
            plotter.add(bn, "BN" +
                        "- #N:" + str(len(bn.dag.nodes())) +
                        ",#E:" + str(len(bn.dag.edges()))
                        )

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
            plotter.add(decomp_bn, "Decompiled BN" +
                        "- #N:" + str(len(decomp_bn.dag.nodes())) +
                        ",#E:" + str(len(decomp_bn.dag.edges()))
                        )

            # --- Compare
            print("------ " + bn_name + " Report - Elim Ord: "
                  + str(i + 1) + "/" + str(len(all_elim_ord)) + " ------")
            print("-> Ori BN #nodes: " + str(len(bn.dag.nodes())) +
                  ", #edges: " + str(len(bn.dag.edges())))
            print("-> Decomp BN #nodes: " + str(len(decomp_bn.dag.nodes()))
                  + ", #edges: " + str(len(decomp_bn.dag.edges())))
            bns_same = is_decomp_bn_same(bn, decomp_bn, elim_ord)
            if not bns_same:
                diff_decomp.append(bn_name)
            if (plot_diff and not bns_same) or plot_all:
                plotter.plot()

    print("******************")
    print("***** Report *****")
    print("******************")
    print("-> Different decompilations: " + str(len(diff_decomp))
          + "/" + str(len(bn_files)))
    print("-> Diff decomp BNs:" + str(diff_decomp))
