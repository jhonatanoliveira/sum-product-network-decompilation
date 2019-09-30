import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from functools import reduce
import argparse
import itertools
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


class BayesianNetwork:

    def __init__(self, dag, var_cardinalities, factors):
        # nx Graph
        self.dag = dag
        # Dict {'A': 2, 'B': 3}
        self.var_cardinalities = var_cardinalities
        self.factors = factors

    def get_elimination_ordering(self, ord_type, print_ordering=True):
        ordering = None
        if ord_type == "rev":
            ordering = list(reversed(list(nx.topological_sort(self.dag))))
        elif ord_type == "top":
            ordering = list(nx.topological_sort(self.dag))
        elif ord_type == "mn":
            ordering = greedy_ordering(
                self.moral_graph(), heuristic_min_neighbour,
                self.var_cardinalities)
        elif ord_type == "mw":
            ordering = greedy_ordering(
                self.moral_graph(), heuristic_min_weight,
                self.var_cardinalities)
        elif "," in ord_type:
            ordering = [v.strip() for v in ord_type.split(",")]
            print_ordering = False
        else:
            raise NotImplementedError
        if print_ordering:
            print("Elimination Ordering: " + str(ordering))

        return ordering

    def moral_graph(self):
        """
        Source: https://networkx.github.io/documentation/stable/
        _modules/networkx/algorithms/moral.html#moral_graph
        """
        H = self.dag.to_undirected()
        for preds in self.dag.pred.values():
            predecessors_combinations = itertools.combinations(preds, r=2)
            H.add_edges_from(predecessors_combinations)
        return H

    def draw(self, show=True):
        NODE_COLOR = "#DA667B"
        color_map = []
        labels = {}
        for node in self.dag.nodes():
            color_map.append(NODE_COLOR)
            labels[node] = node
        draw_graph(self.dag, show, color_map, labels)


class Factor:

    def __init__(self, variables, cardinalities, values):
        self.variables = variables
        self.total_variables = len(self.variables)
        self.values = values
        self.total_values = len(self.values)
        self.cardinalities = cardinalities
        # Compute strides
        self._strides = {}
        for idx, variable in enumerate(self.variables):
            if idx == 0:
                self._strides[variable] = 1
            else:
                last_variable = self.variables[idx - 1]
                self._strides[variable] = self._strides[last_variable] *\
                    self.cardinalities[last_variable]

    def stride(self, variable):
        return self._strides[variable] if variable in self._strides else 0


class ArithmeticCircuit:

    def __init__(self, dag):
        self.dag = dag

    def draw(self, show=True):
        SUM_COLOR, PROD_COLOR, IND_COLOR, PROB_COLOR, TERM_COLOR =\
            "#ECE4B7", "#D9DD92", "#EABE7C", "#DD6031", "#6B717E"
        SUM_LABEL, PROD_LABEL = "+", "x"
        color_map = []
        labels = {}
        for node in self.dag.nodes():
            if "+" in node:
                color_map.append(SUM_COLOR)
                labels[node] = SUM_LABEL
            elif "*" in node:
                color_map.append(PROD_COLOR)
                labels[node] = PROD_LABEL
            elif "I" in node:
                color_map.append(IND_COLOR)
                labels[node] = node.replace("I(", "").replace(")", "")
            elif "P" in node:
                color_map.append(PROB_COLOR)
                labels[node] = ""
            elif "T" in node:
                color_map.append(TERM_COLOR)
                tmp = node.replace("T(", "").replace(")", "")
                labels[node] = tmp[:tmp.index("-")]
        draw_graph(self.dag, show, color_map, labels)


def get_bn_from_file(bn_file_name):
    # Construct the DAG
    bn_dag = nx.DiGraph()
    var_cardinalities = {}
    with open(bn_file_name) as bn_file:
        vars_line = bn_file.readline().replace("\n", "")
        for var_line in vars_line.split(","):
            var, cardinality = var_line.split(":")
            cardinality = int(cardinality)
            bn_dag.add_node(var)
            var_cardinalities[var] = cardinality
        edge_line = bn_file.readline().replace("\n", "")
        while edge_line:
            parent, child = edge_line.split(",")
            bn_dag.add_edge(parent, child)
            edge_line = bn_file.readline().replace("\n", "")
    # Construct the factors
    factors = {}
    for node in bn_dag.nodes():
        variables = [node] + list(bn_dag.predecessors(node))
        cardinalities = {var: var_cardinalities[var] for var in variables}
        total_values = reduce(
            (lambda x, y: x * y), [cardinalities[k] for k in cardinalities])
        values = [0 for _ in range(0, total_values)]
        factor = Factor(variables, cardinalities, values)
        factors[node] = factor
    return BayesianNetwork(bn_dag, var_cardinalities, factors)


def ac_factor_product(ac_dag, factor1, factor2, node_counters):
    prod_variables = factor1.variables + [
        v for v in factor2.variables if v not in factor1.variables]
    prod_cardinalities = {
        **factor1.cardinalities,
        **{v: factor2.cardinalities[v] for v in factor2.cardinalities
           if v not in factor1.variables}
    }
    prod_total_values = reduce(
            (lambda x, y: x * y),
            [prod_cardinalities[k] for k in prod_cardinalities])
    prod_values = [0 for _ in range(0, prod_total_values)]
    prod_factor = Factor(prod_variables, prod_cardinalities, prod_values)

    j, k = 0, 0
    assignment = [0 for _ in range(prod_factor.total_variables)]
    for i in range(0, prod_factor.total_values):
        n_prod = "*" + str(node_counters[1])
        node_counters[1] += 1  # increment product counter
        prod_factor.values[i] = n_prod
        ac_dag.add_node(n_prod)
        ac_dag.add_edge(n_prod, factor1.values[j])
        ac_dag.add_edge(n_prod, factor2.values[k])
        for l, variable in enumerate(prod_factor.variables):
            assignment[l] += 1
            if assignment[l] == prod_factor.cardinalities[variable]:
                assignment[l] = 0
                j -= (prod_factor.cardinalities[variable] - 1) *\
                    factor1.stride(variable)
                k -= (prod_factor.cardinalities[variable] - 1) *\
                    factor2.stride(variable)
            else:
                j += factor1.stride(variable)
                k += factor2.stride(variable)
                break
    return prod_factor


def ac_factor_sum(ac_dag, factor1, variable, node_counters):
    sum_variables = [v for v in factor1.variables if v != variable]
    sum_cardinalities = {v: factor1.cardinalities[v]
                         for v in factor1.variables if v != variable}
    sum_total_values = reduce(
            (lambda x, y: x * y),
            [sum_cardinalities[k] for k in sum_cardinalities])\
        if len(sum_cardinalities) > 0 else 1
    sum_values = [0 for _ in range(0, sum_total_values)]
    sum_factor = Factor(sum_variables, sum_cardinalities, sum_values)

    j = 0
    assignment = [0 for _ in range(factor1.total_variables)]
    for i in range(0, factor1.total_values):
        if sum_factor.values[j] == 0:
            n_sum = "+" + str(node_counters[0])
            node_counters[0] += 1  # increment sum counter
            sum_factor.values[j] = n_sum
            ac_dag.add_node(n_sum)
            ac_dag.add_edge(n_sum, factor1.values[i])
        else:
            ac_dag.add_edge(sum_factor.values[j], factor1.values[i])
        for l, variable in enumerate(factor1.variables):
            assignment[l] += 1
            if assignment[l] == factor1.cardinalities[variable]:
                assignment[l] = 0
                j -= (factor1.cardinalities[variable] - 1) *\
                    sum_factor.stride(variable)
            else:
                j += sum_factor.stride(variable)
                break

    return sum_factor


def ac_leaves(ac_dag, bn, node_counters):
    for variable in bn.factors:
        factor = bn.factors[variable]
        assignment = [0 for _ in range(0, factor.total_variables)]
        for i in range(factor.total_values):  # per row
            indc_var = "I(" +\
                factor.variables[0] + "-" + str(assignment[0]) + ")"
            param_var = "P(" +\
                factor.variables[0] + str(assignment[0]) + "|" +\
                ",".join([factor.variables[k] + str(assignment[k])
                          for k in range(1, factor.total_variables)]) + ")"\
                if factor.total_variables > 1\
                else "P(" + factor.variables[0] + str(assignment[0]) + ")"
            ac_dag.add_node(indc_var, node_color='b')
            ac_dag.add_node(param_var, node_color='g')
            n_prod = "*" + str(node_counters[1])
            node_counters[1] += 1  # increment product counter
            factor.values[i] = n_prod
            ac_dag.add_node(n_prod)
            ac_dag.add_edge(n_prod, indc_var)
            ac_dag.add_edge(n_prod, param_var)
            for l, var in enumerate(factor.variables):  # per column
                assignment[l] += 1
                if assignment[l] == factor.cardinalities[var]:
                    assignment[l] = 0
                else:
                    break


def compile_variable_elimination(bn, elim_ord):
    ac_dag = nx.DiGraph()
    node_counters = [0, 0]  # [number of sum, number of products]
    ac_leaves(ac_dag, bn, node_counters)
    curr_bn_factors = [bn.factors[k] for k in bn.factors]
    for variable in elim_ord:
        var_factors = list(
            filter(lambda factor: variable in factor.variables,
                   curr_bn_factors))
        curr_bn_factors = list(
            filter(lambda factor: variable not in factor.variables,
                   curr_bn_factors))
        prod_factor = None
        for factor in var_factors:
            if prod_factor is None:
                prod_factor = factor
            else:
                pass
                prod_factor = ac_factor_product(
                    ac_dag, prod_factor, factor, node_counters)
        sum_factor = ac_factor_sum(
            ac_dag, prod_factor, variable, node_counters)
        curr_bn_factors.append(sum_factor)

    return ArithmeticCircuit(ac_dag)


def draw_graph(graph, show, color_map=None, labels=None):
    pos = graphviz_layout(graph, prog='dot')
    nx.draw(graph, pos, with_labels=True, arrows=True,
            node_color=color_map, node_size=400,
            font_size=8, labels=labels)
    if show:
        plt.show()


def greedy_ordering(moral_graph, heuristic, cardinalities):
    mgraph = moral_graph.copy()
    to_visit = [n for n in mgraph.nodes()]
    ordering = []
    while len(to_visit) > 0:
        min_n = to_visit[0]
        min_score = heuristic(min_n, mgraph, cardinalities)
        for n in to_visit[1:]:
            curr_score = heuristic(n, mgraph, cardinalities)
            if curr_score < min_score:
                min_n = n
                min_score = curr_score
        ordering.append(min_n)
        mgraph.add_edges_from(
            list(itertools.combinations(list(mgraph.neighbors(min_n)), 2))
            )
        mgraph.remove_node(min_n)
        to_visit.remove(min_n)
    return ordering


def heuristic_min_neighbour(node, moral_graph, cardinalities):
    return len(list(moral_graph.neighbors(node)))


def heuristic_min_weight(node, moral_graph, cardinalities):
    neighbours = list(moral_graph.neighbors(node))
    return reduce(
            (lambda x, y: x * y),
            [cardinalities[n] for n in neighbours])\
        if len(neighbours) > 0 else 1


def remove_parameters_ac(ac):
    mod_dag = ac.dag.copy()
    for node in [n for n in mod_dag.nodes()]:
        if "P" in node:
            mod_dag.remove_node(node)
    return ArithmeticCircuit(mod_dag)


def remove_barren_prod_ac(ac):
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


def add_terminal_node_ac(ac):
    mod_dag = ac.dag.copy()
    terminal_counter = 0
    for node in [n for n in mod_dag.nodes()]:
        children = list(mod_dag.successors(node))
        if "+" in node and len(children) > 0\
                and all(["I" in c for c in children]):
            a_child = children[0]
            term_node = a_child.replace("I(", "").replace(")", "")
            term_node = term_node[:term_node.index("-")]
            term_node += "-" + str(terminal_counter)
            term_node = "T(" + term_node + ")"
            terminal_counter += 1
            mod_dag.add_node(term_node)
            parents = list(mod_dag.predecessors(node))
            mod_dag.add_edges_from(
                list(itertools.product(parents, [term_node])))
            mod_dag.remove_node(node)
    for node in [n for n in mod_dag.nodes()]:
        children = list(mod_dag.successors(node))
        parents = list(mod_dag.predecessors(node))
        if len(children) == 0 and len(parents) == 0:
            mod_dag.remove_node(node)

    return ArithmeticCircuit(mod_dag)


def remove_indicators_ac(ac):
    mod_dag = ac.dag.copy()
    for node in [n for n in mod_dag.nodes()]:
        if "I" in node:
            all_sum_prod = True
            for parent in mod_dag.predecessors(node):
                for sibling in mod_dag.successors(parent):
                    if sibling != node and not (
                            "+" in sibling or
                            "*" in sibling or "T" in sibling):
                        all_sum_prod = False
            if all_sum_prod:
                mod_dag.remove_node(node)
    return ArithmeticCircuit(mod_dag)


def draw_subplot_graphs(graphs, subtitles=None, main_title=None):
    total_graphs = len(graphs)
    subplot_amt = math.ceil(math.sqrt(total_graphs))
    subplot_rows = (subplot_amt - 1)\
        if (subplot_amt * subplot_amt - total_graphs) >= subplot_amt\
        else subplot_amt
    plt.figure()
    for i, graph in enumerate(graphs):
        plt.subplot(subplot_rows, subplot_amt, i + 1)
        graph.draw(show=False)
        if subtitles:
            plt.title(subtitles[i])
        plt.grid(True)
    if main_title:
        plt.suptitle(main_title)
    plt.show()


def remove_doubled_prod_ac(ac):
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


def get_bn(ac, contract_scope=True):
    # algorithm for detecting conditioning
    def _is_conditioning(anc, sum_node, dag):
        not_in_desc = False
        for child in dag.successors(anc):
            if node not in nx.algorithms.dag.descendants(dag, child):
                not_in_desc = True
                break
        return not_in_desc
    # compute node scopes
    scopes = {}
    for node in reversed(list(nx.topological_sort(ac.dag))):
        if "+" in node or "*" in node:
            scope = set()
            for child in ac.dag.successors(node):
                scope = scope.union(scopes[child])
            scopes[node] = scope
        elif "T" in node:
            tmp = node.replace("T(", "").replace(")", "")
            scope = set([tmp[:tmp.index("-")]])
            scopes[node] = scope
        else:
            raise NotImplementedError("Node scope not defined")
    # assign latent variable
    latent_vars = {}
    lv_scopes = []
    lv_assigned = []
    sum_lv_idx = 0
    for node in reversed(list(nx.topological_sort(ac.dag))):
        if "T" in node:
            latent_vars[node] = list(scopes[node])[0]
        elif "+" in node:
            if contract_scope:
                if scopes[node] in lv_scopes:
                    idx_lv = lv_scopes.index(scopes[node])
                    latent_vars[node] = lv_assigned[idx_lv]
                else:
                    new_lv_assignment = "Zs-" + str(sum_lv_idx)
                    lv_scopes.append(scopes[node])
                    lv_assigned.append(new_lv_assignment)
                    latent_vars[node] = new_lv_assignment
            else:
                latent_vars[node] = "Zs-" + str(sum_lv_idx)
            sum_lv_idx += 1
    # construct bn DAG
    bn_dag = nx.DiGraph()
    for node in reversed(list(nx.topological_sort(ac.dag))):
        if "+" in node or "T" in node:
            for anc in nx.algorithms.dag.ancestors(ac.dag, node):
                if "+" in anc and _is_conditioning(anc, node, ac.dag):
                    parent = latent_vars[anc]
                    child = latent_vars[node]
                    bn_dag.add_edge(parent, child)

    return BayesianNetwork(bn_dag, None, None)


def remove_duplicated_prod_ac(ac):
    mod_dag = ac.dag.copy()
    seen_prods, seen_children = [], []  # Key, Value
    for node in [n for n in mod_dag.nodes()]:
        if "*" in node:
            ch = set(list(mod_dag.successors(node)))
            if ch in seen_children:
                parents = list(mod_dag.predecessors(node))
                prod_to_connec = seen_prods[seen_children.index(ch)]
                mod_dag.add_edges_from(
                    itertools.product(parents, [prod_to_connec]))
                mod_dag.remove_node(node)
            else:
                seen_children.append(ch)
                seen_prods.append(node)
    return ArithmeticCircuit(mod_dag)


def reconstruct_bn(bn, elim_ord, plot=False):
    graphs = []
    graphs_subtitles = []

    graphs.append(bn)
    graphs_subtitles.append("BN")

    ac = compile_variable_elimination(bn, elim_ord)
    graphs.append(ac)
    graphs_subtitles.append("AC")

    ac = remove_parameters_ac(ac)
    graphs.append(ac)
    graphs_subtitles.append("Remove Parameters")

    ac = remove_barren_prod_ac(ac)
    graphs.append(ac)
    graphs_subtitles.append("Remove Barren Products")

    ac = add_terminal_node_ac(ac)
    graphs.append(ac)
    graphs_subtitles.append("Add Terminal Nodes")

    ac = remove_indicators_ac(ac)
    graphs.append(ac)
    graphs_subtitles.append("Remove Indicators")

    ac = remove_barren_prod_ac(ac)
    graphs.append(ac)
    graphs_subtitles.append("Remove Barren Products")

    ac = remove_doubled_prod_ac(ac)
    graphs.append(ac)
    graphs_subtitles.append("Remove Doubled Products")

    ac = remove_duplicated_prod_ac(ac)
    graphs.append(ac)
    graphs_subtitles.append("Remove Duplicated Products")

    bn = get_bn(ac, contract_scope=True)
    graphs.append(bn)
    graphs_subtitles.append("Reconstructed BN")

    main_title = "Elimination Ordering: " + ",".join(elim_ord)

    if plot:
        draw_subplot_graphs(graphs, graphs_subtitles, main_title)

    return bn


def remove_observables(bn):
    mod_dag = bn.dag.copy()
    mod_dag.remove_nodes_from([n for n in mod_dag.nodes()
                               if len(list(mod_dag.successors(n))) == 0])

    return BayesianNetwork(mod_dag, None, None)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile a BN into an AC")
    parser.add_argument("--bn", type=str, help="Path to '.bn' file.")
    parser.add_argument("--elim-ord", type=str,
                        help="Elimination ordering to be used.", default="rev")
    args = parser.parse_args()

    bn = get_bn_from_file(args.bn)

    all_elim_ord = list(nx.algorithms.dag.all_topological_sorts(bn.dag))
    print("Total elimination oderings: " + str(len(all_elim_ord)))
    for elim_ord in all_elim_ord:
        elim_ord = list(reversed(elim_ord))
        rec_bn = reconstruct_bn(bn, elim_ord, plot=True)
