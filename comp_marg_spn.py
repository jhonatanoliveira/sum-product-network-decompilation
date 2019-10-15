from models import SumProductNetwork
from convert_ac2spn import simplify


def compile_marginalized_spn(spn, variables, subplot_drawer):

    spn_dag = spn.dag.copy()
    spn = SumProductNetwork(spn_dag)

    for node in [n for n in spn.dag.nodes()]:
        if (("I" in node) or ("T" in node))\
                and any([var in node for var in variables]):
            spn.dag.remove_node(node)

    return simplify(spn, subplot_drawer)
