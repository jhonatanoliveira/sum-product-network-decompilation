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
import time
from joblib import Parallel, delayed


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
        var_cardinalities = {"X" + str(var): 2
                             for var in range(num_vars)}
    bns = []
    i = 0
    for adj_matrix in connected_dags(num_vars):
        dag = nx.convert_matrix.from_numpy_matrix(
            adj_matrix, create_using=nx.DiGraph)
        dag.graph["name"] = "BN-" + str(i)
        str_node_labels = {old_lbl: "X" + str(old_lbl)
                           for old_lbl in dag.nodes()}
        dag = nx.relabel_nodes(dag, str_node_labels)
        i += 1
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


def permute(A, permutation):
    assert A.shape[0] == len(permutation)
    assert A.shape[1] == len(permutation)
    return A[permutation, :][:, permutation]


def is_exaustive_isomorphic(
        A, B, matched_idx_A=None, matched_idx_B=None):

    if matched_idx_A is None:
        matched_idx_A = []
    if matched_idx_B is None:
        matched_idx_B = []

    assert A.shape[0] == A.shape[1]
    assert B.shape[0] == B.shape[1]
    assert A.shape[0] == B.shape[0]

    N = A.shape[0]

    assert len(matched_idx_A) == len(matched_idx_B)
    assert len(matched_idx_A) == len(set(matched_idx_A))
    assert len(matched_idx_B) == len(set(matched_idx_B))
    if len(matched_idx_A):
        assert min(matched_idx_A + matched_idx_B) >= 0
        assert max(matched_idx_A + matched_idx_B) < N

    other_idx_A = [i for i in range(N) if i not in matched_idx_A]
    other_idx_B = [i for i in range(N) if i not in matched_idx_B]

    A = permute(A, other_idx_A + matched_idx_A)
    B = permute(B, other_idx_B + matched_idx_B)

    for p in itertools.permutations(range(len(other_idx_B))):
        if np.all(
            A == permute(
                B, list(p) + list(range(len(other_idx_B), N)))):
            return True

    return False


def directed_moralize(dag, elim_ord):
    for preds in dag.pred.values():
        for pred_comb in itertools.combinations(preds, r=2):
            right_dir = pred_comb \
                if elim_ord.index(pred_comb[0]) > elim_ord.index(pred_comb[1])\
                else (pred_comb[1], pred_comb[0])
            dag.add_edge(right_dir[0], right_dir[1])


def higher_order_moralization(dag, elim_ord):
    while True:
        mor_dag = dag.copy()
        directed_moralize(mor_dag, elim_ord)
        # OBS: the order of the difference matters: mor_dag - dag is different
        # from dag - mor_dag
        diff_graph = nx.algorithms.operators.binary.difference(mor_dag, dag)
        if len(list(diff_graph.edges())) == 0:
            break
        else:
            dag = mor_dag

    return dag


def is_decomp_bn_same(ori_bn, decomp_bn, elim_ord):
    triang_ori = triangulate(ori_bn.dag, elim_ord)
    ori_nds_ord = triang_ori.nodes()
    ori_leaf_idxs = [i for i, n in enumerate(ori_nds_ord)
                     if ori_bn.dag.successors(n) == 0]
    undirected_decomp = decomp_bn.dag.to_undirected()
    decomp_nds_ord = undirected_decomp.nodes()
    decomp_leaf_idxs = [i for i, n in enumerate(decomp_nds_ord)
                        if decomp_bn.dag.successors(n) == 0]
    adj_ori = nx.convert_matrix.to_numpy_matrix(triang_ori)
    adj_decomp = nx.convert_matrix.to_numpy_matrix(undirected_decomp)

    return is_exaustive_isomorphic(adj_ori, adj_decomp,
                                   matched_idx_A=ori_leaf_idxs,
                                   matched_idx_B=decomp_leaf_idxs)


def printProgressBar(iteration, total, prefix='', suffix='',
                     decimals=1, length=100, fill='█',
                     printEnd="\r"):
    """
    https://stackoverflow.com/questions/3173320/
    text-progress-bar-in-the-console
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals
                      in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(
        100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (
        prefix, bar, percent, suffix), end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


def run_bn_exp(bn, bn_idx, elim_ord_type):
    bn_name = bn.dag.graph['name']

    all_elim_ord = EliminationOrdering.get_elimination_ordering(
        bn, elim_ord_type)
    all_elim_ord = all_elim_ord if type(all_elim_ord[0]) is list\
        else [all_elim_ord]
    for i, elim_ord in enumerate(all_elim_ord):
        plotter = SubplotDrawer("Compiling and Decompiling BN - " +
                                bn_name)
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

        # --- Compare and print/plot
        bns_same = is_decomp_bn_same(bn, decomp_bn, elim_ord)
        if print_compare:
            print("------ " + bn_name + " Report - Elim Ord: "
                  + str(i + 1) + "/" + str(len(all_elim_ord)) + " ------")
            print("-> Ori BN #nodes: " + str(len(bn.dag.nodes())) +
                  ", #edges: " + str(len(bn.dag.edges())))
            print("-> Decomp BN #nodes: " + str(len(decomp_bn.dag.nodes()))
                  + ", #edges: " + str(len(decomp_bn.dag.edges())))
            print(">>> Decompilation same as original: " + str(bns_same))
        printProgressBar(bn_idx + 1, total_bns,
                         prefix='Progress:',
                         suffix='Complete', length=50)
        if not bns_same:
            diff_decomp.append(bn_name)
        if (plot_diff and not bns_same) or plot_all:
            plotter.plot()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Experiment for recoverying BN")
    parser.add_argument(
        "--bn", type=str, help="Path to '.bn' file or folder name.\
            For generated BNs, use 'all-bns'.",
        required=True)
    parser.add_argument(
        "--elim-ord", type=str, help="Elimination ordering type.",
        choices=["rev", "top", "mn", "mw", "all-rev"], required=True)
    parser.add_argument(
        "--plot-all",
        help="Plot all graphs.",
        action="store_true")
    parser.add_argument(
        "--all-bns-nvars",
        help="Number of variables for all BN generation.",
        type=int)
    args = parser.parse_args()

    # Parameters
    bn_path = args.bn
    is_collapse_sums = False
    elim_ord_type = args.elim_ord
    comp_marg_spn = "internal"
    plot_diff = True
    plot_all = args.plot_all
    print_compare = False
    run_parallel = False

    bns = []
    if bn_path == "all-bns":
        bns = generate_all_bns(args.all_bns_nvars)
    else:
        bn_files = []
        is_file = bn_path.endswith(".bn")
        if is_file:
            bn_files = [bn_path]
        else:
            bn_files = [os.path.join("bns", name)
                        for name in os.listdir(bn_path)]
        bns = [
            BayesianNetwork.get_bn_from_file(bn_file_name)
            for bn_file_name in bn_files
        ]

    diff_decomp = []
    total_bns = len(bns)
    start_time = time.time()

    if run_parallel:
        Parallel(n_jobs=10)(delayed(run_bn_exp)(bn, bn_idx, elim_ord_type)
                            for bn_idx, bn in enumerate(bns))
    else:
        for bn_idx, bn in enumerate(bns):
            run_bn_exp(bn, bn_idx, elim_ord_type)

    end_time = time.time()
    print("******************")
    print("***** Report *****")
    print("******************")
    print("-> Different decompilations: " + str(len(diff_decomp)))
    print("-> Diff decomp BNs:" + str(diff_decomp))
    print("Execution: " + str(end_time - start_time))
