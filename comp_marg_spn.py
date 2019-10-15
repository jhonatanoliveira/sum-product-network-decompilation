import itertools
from models import SumProductNetwork
from convert_ac2spn import simplify


def compile_marginalized_spn(spn, variables,
                             collapse_sums=False, subplot_drawer=None):

    spn_dag = spn.dag.copy()
    spn = SumProductNetwork(spn_dag)

    for node in [n for n in spn.dag.nodes()]:
        if (("I" in node) or ("T" in node))\
                and any([var in node for var in variables]):
            spn.dag.remove_node(node)

    spn = simplify(spn, subplot_drawer)

    if collapse_sums:
        for node in [n for n in spn.dag.nodes()]:
            if "+" in node and all(
                    ["+" in child for child in spn.dag.successors(node)]):
                new_children = []
                for child in spn.dag.successors(node):
                    new_children += list(spn.dag.successors(child))
                # OBS: remove nodes first then add edges later,
                # otherwise wrong nodes will be removed because of new edges
                spn.dag.remove_nodes_from(list(spn.dag.successors(node)))
                spn.dag.add_edges_from(
                    list(itertools.product([node], new_children))
                )

    return spn
