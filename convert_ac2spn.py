from models import SumProductNetwork
import itertools
import networkx as nx


def redistribute_parameters(ac):
    for node in [n for n in ac.dag.nodes()]:
        if "P" in node:
            ac.dag.remove_node(node)
    return ac


def remove_barren_prod(ac):
    found = True
    while found:
        found = False
        for node in [n for n in ac.dag.nodes()]:
            children = list(ac.dag.successors(node))
            if "*" in node and len(children) == 1:
                parents = list(ac.dag.predecessors(node))
                ac.dag.add_edges_from(
                    list(itertools.product(parents, children)))
                ac.dag.remove_node(node)
                found = True
    return ac


def remove_doubled_prod(ac):
    found = True
    while found:
        found = False
        for node in [n for n in ac.dag.nodes()]:
            if "*" in node:
                parents = list(ac.dag.predecessors(node))
                for parent in parents:
                    if "*" in parent:
                        children = list(ac.dag.successors(node))
                        ac.dag.add_edges_from(
                            list(itertools.product([parent], children)))
                        ac.dag.remove_edge(parent, node)
                        found = True
                if len(list(ac.dag.predecessors(node))) == 0:
                    ac.dag.remove_node(node)
    return ac


def is_different(dag1, dag2):
    diff_dag = None
    is_diff = False
    try:
        diff_dag = nx.difference(dag1, dag2)
        is_diff = len(diff_dag.edges()) != 0
    except nx.exception.NetworkXError:
        is_diff = True

    return is_diff


def _apply_transformation_check_diff(spn, trans_func,
                                     drawer, subtitle):
    is_diff = False
    ori_spn = spn.copy()
    trans_spn = trans_func(ori_spn)
    if is_different(ori_spn.dag, trans_spn.dag):
        is_diff = True
        drawer.add(trans_spn, subtitle)
    return trans_spn, is_diff


def convert_ac2spn(ac, subplot_drawer=None):

    spn_dag = ac.dag.copy()
    spn = SumProductNetwork(spn_dag)

    spn = redistribute_parameters(spn)
    subplot_drawer.add(spn, "Redistribute Parameters")

    has_changed = True
    while has_changed:
        has_changed = False

        spn, is_diff = _apply_transformation_check_diff(
            spn, remove_barren_prod, subplot_drawer, "Remove Barren Products")
        has_changed = (has_changed or is_diff)

    return spn
