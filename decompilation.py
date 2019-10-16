import networkx as nx
from models import BayesianNetwork, Factor


def spn_assign_lv_contract(spn_dag, spn_scopes, contract_scope=True):
    latent_vars = {}
    lv_scopes = []
    lv_assigned = []
    sum_lv_idx = 0
    for node in reversed(list(nx.topological_sort(spn_dag))):
        if "T" in node:
            latent_vars[node] = list(spn_scopes[node])[0]
        elif "+" in node:
            if contract_scope:
                if spn_scopes[node] in lv_scopes:
                    idx_lv = lv_scopes.index(spn_scopes[node])
                    latent_vars[node] = lv_assigned[idx_lv]
                else:
                    new_lv_assignment = "Zs-" + str(sum_lv_idx)
                    sum_lv_idx += 1
                    lv_scopes.append(spn_scopes[node])
                    lv_assigned.append(new_lv_assignment)
                    latent_vars[node] = new_lv_assignment
            else:
                latent_vars[node] = "Zs-" + str(sum_lv_idx)
                sum_lv_idx += 1
    return latent_vars


def spn_assign_lv_sum_depth(spn_dag, spn_scopes):
    sum_depth_layers = make_sum_depth_layers(spn_dag)
    latent_vars = {}
    sum_lv_idx = 0
    for layer in sum_depth_layers:
        lv_scopes = []
        lv_assigned = []
        for node in layer:
            if "T" in node:
                latent_vars[node] = list(spn_scopes[node])[0]
            elif "+" in node:
                if spn_scopes[node] in lv_scopes:
                    idx_lv = lv_scopes.index(spn_scopes[node])
                    latent_vars[node] = lv_assigned[idx_lv]
                else:
                    new_lv_assignment = "Zs-" + str(sum_lv_idx)
                    sum_lv_idx += 1
                    lv_scopes.append(spn_scopes[node])
                    lv_assigned.append(new_lv_assignment)
                    latent_vars[node] = new_lv_assignment
    return latent_vars


def make_sum_depth_layers(spn_dag):
    layers = []
    root_layer = [n for n in spn_dag.nodes()
                  if len(list(spn_dag.predecessors(n))) == 0]
    layers.append(root_layer)
    to_search_sum_layers = [root_layer[:]]
    marked = []
    while len(to_search_sum_layers) > 0:
        curr_layer = to_search_sum_layers.pop(0)
        to_search_in_layer = curr_layer[:]
        marked += curr_layer[:]
        next_sum_layer = []
        while len(to_search_in_layer) > 0:
            curr_node = to_search_in_layer.pop(0)
            next_nodes = []
            for child in spn_dag.successors(curr_node):
                if child not in next_sum_layer and\
                        all([p in marked
                             for p in spn_dag.predecessors(child)]):
                    if "T" in child or "+" in child:
                        next_sum_layer.append(child)
                    elif "*" in child:
                        next_nodes.append(child)
                    else:
                        raise ValueError
            to_search_in_layer += next_nodes
            marked += next_nodes
        if len(next_sum_layer) > 0:
            layers.append(next_sum_layer)
            to_search_sum_layers.append(next_sum_layer)
    return layers


# algorithm for detecting conditioning
def is_conditioning(anc, sum_node, dag):
    not_in_desc = False
    for child in dag.successors(anc):
        if sum_node not in nx.algorithms.dag.descendants(dag, child):
            not_in_desc = True
            break
    return not_in_desc


def decompile(spn, comp_assumption=None, lv_cardinalities=None):
    if comp_assumption == "ve":
        assign_lv = spn_assign_lv_sum_depth
    else:
        assign_lv = spn_assign_lv_contract
    scopes = spn.scopes()
    # assign latent variable
    latent_vars = assign_lv(spn.dag, scopes)
    # construct bn DAG
    bn_dag = nx.DiGraph()
    for node in reversed(list(nx.topological_sort(spn.dag))):
        if "+" in node or "T" in node:
            for anc in nx.algorithms.dag.ancestors(spn.dag, node):
                if "+" in anc and is_conditioning(anc, node, spn.dag):
                    parent = latent_vars[anc]
                    child = latent_vars[node]
                    bn_dag.add_edge(parent, child)
    # Assume cardinalities 2 for latent variables cardinalities
    if lv_cardinalities is None:
        lv_cardinalities = {
            var: 2 for var in bn_dag.nodes()
            if var not in spn.var_cardinalities}
    cardinalities = {**spn.var_cardinalities, **lv_cardinalities}
    # Construct the factors
    factors = Factor.construct_factors(bn_dag, cardinalities)

    return BayesianNetwork(bn_dag, cardinalities, factors)
