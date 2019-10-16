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


def add_terminal_node(ac):
    terminal_counter = 0
    for node in [n for n in ac.dag.nodes()]:
        children = list(ac.dag.successors(node))
        if "+" in node and len(children) > 0\
                and all(["I" in c for c in children]):
            a_child = children[0]
            term_node = a_child.replace("I(", "").replace(")", "")
            term_node = term_node[:term_node.index("-")]
            term_node += "-" + str(terminal_counter)
            term_node = "T(" + term_node + ")"
            terminal_counter += 1
            ac.dag.add_node(term_node)
            parents = list(ac.dag.predecessors(node))
            ac.dag.add_edges_from(
                list(itertools.product(parents, [term_node])))
            ac.dag.remove_node(node)
    for node in [n for n in ac.dag.nodes()]:
        children = list(ac.dag.successors(node))
        parents = list(ac.dag.predecessors(node))
        if len(children) == 0 and len(parents) == 0:
            ac.dag.remove_node(node)


def remove_prod_of_prod(ac):
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


def lump_prod_over_same_children(ac):
    seen_prods, seen_children = [], []  # Key, Value
    for node in [n for n in ac.dag.nodes()]:
        if "*" in node:
            ch = set(list(ac.dag.successors(node)))
            if ch in seen_children:
                parents = list(ac.dag.predecessors(node))
                prod_to_connec = seen_prods[seen_children.index(ch)]
                ac.dag.add_edges_from(
                    itertools.product(parents, [prod_to_connec]))
                ac.dag.remove_node(node)
            else:
                seen_children.append(ch)
                seen_prods.append(node)


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
    trans_spn = spn.copy()
    trans_func(trans_spn)
    if is_different(spn.dag, trans_spn.dag):
        is_diff = True
        if drawer:
            drawer.add(trans_spn, subtitle)
    return trans_spn, is_diff


def simplify(spn, subplot_drawer=None):

    has_changed = True
    while has_changed:
        has_changed = False

        spn, is_diff = _apply_transformation_check_diff(
            spn, add_terminal_node, subplot_drawer, "Add Terminal Nodes")
        has_changed = (has_changed or is_diff)

        spn, is_diff = _apply_transformation_check_diff(
            spn, remove_barren_prod, subplot_drawer, "Remove Barren Products")
        has_changed = (has_changed or is_diff)

        spn, is_diff = _apply_transformation_check_diff(
            spn, remove_prod_of_prod, subplot_drawer,
            "Remove Products of Products")
        has_changed = (has_changed or is_diff)

        spn, is_diff = _apply_transformation_check_diff(
            spn, lump_prod_over_same_children, subplot_drawer,
            "Lump Products Over Same Children")
        has_changed = (has_changed or is_diff)

    return spn


def convert_ac2spn(ac, subplot_drawer=None):

    spn_dag = ac.dag.copy()
    spn = SumProductNetwork(spn_dag, ac.var_cardinalities)

    spn = redistribute_parameters(spn)
    subplot_drawer.add(spn, "Redistribute Parameters")

    spn = simplify(spn, subplot_drawer)

    return spn
