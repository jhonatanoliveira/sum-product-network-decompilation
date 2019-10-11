from models import ArithmeticCircuit, SumProductNetwork
import itertools


def redistribute_parameters(mod_dag):
    for node in [n for n in mod_dag.nodes()]:
        if "P" in node:
            mod_dag.remove_node(node)
    return mod_dag


def remove_barren_prod(ac):
    mod_dag = ac.dag.copy()
    found = True
    while found:
        found = False
        for node in [n for n in mod_dag.nodes()]:
            children = list(mod_dag.successors(node))
            if "*" in node and len(children) == 1:
                parents = list(mod_dag.predecessors(node))
                mod_dag.add_edges_from(
                    list(itertools.product(parents, children)))
                mod_dag.remove_node(node)
                found = True
    return ArithmeticCircuit(mod_dag)


def remove_doubled_prod(ac):
    mod_dag = ac.dag.copy()
    found = True
    while found:
        found = False
        for node in [n for n in mod_dag.nodes()]:
            if "*" in node:
                parents = list(mod_dag.predecessors(node))
                for parent in parents:
                    if "*" in parent:
                        children = list(mod_dag.successors(node))
                        mod_dag.add_edges_from(
                            list(itertools.product([parent], children)))
                        mod_dag.remove_edge(parent, node)
                        found = True
                if len(list(mod_dag.predecessors(node))) == 0:
                    mod_dag.remove_node(node)
    return ArithmeticCircuit(mod_dag)


def convert_ac2spn(ac):

    spn_dag = ac.dag.copy()
    spn_dag = redistribute_parameters(spn_dag)

    return SumProductNetwork(spn_dag)
