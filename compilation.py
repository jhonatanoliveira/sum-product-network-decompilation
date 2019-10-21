import itertools
import networkx as nx
from functools import reduce
from models import ArithmeticCircuit, Factor


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

    return ArithmeticCircuit(ac_dag, bn.var_cardinalities)


class EliminationOrdering:

    @staticmethod
    def moral_graph(dag):
        """
        Source: https://networkx.github.io/documentation/stable/
        _modules/networkx/algorithms/moral.html#moral_graph
        """
        H = dag.to_undirected()
        for preds in dag.pred.values():
            predecessors_combinations = itertools.combinations(preds, r=2)
            H.add_edges_from(predecessors_combinations)
        return H

    @staticmethod
    def get_elimination_ordering(bn, ord_type):
        dag = bn.dag
        var_cardinalities = bn.var_cardinalities
        moral_graph = EliminationOrdering.moral_graph(dag)
        ordering = None
        if ord_type == "rev":
            ordering = list(reversed(list(nx.topological_sort(dag))))
        elif ord_type == "top":
            ordering = list(nx.topological_sort(dag))
        elif ord_type == "mn":
            ordering = EliminationOrdering.greedy_ordering(
                moral_graph, EliminationOrdering.heuristic_min_neighbour,
                var_cardinalities)
        elif ord_type == "mw":
            ordering = EliminationOrdering.greedy_ordering(
                moral_graph, EliminationOrdering.heuristic_min_weight,
                var_cardinalities)
        elif ord_type == "all-rev":
            all_ordering = list(
                nx.algorithms.dag.all_topological_sorts(bn.dag))
            ordering = [list(reversed(elim_ord)) for elim_ord in all_ordering]
        elif "," in ord_type:
            ordering = [v.strip() for v in ord_type.split(",")]
        else:
            raise NotImplementedError

        return ordering

    @staticmethod
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

    @staticmethod
    def heuristic_min_neighbour(node, moral_graph, cardinalities):
        return len(list(moral_graph.neighbors(node)))

    @staticmethod
    def heuristic_min_weight(node, moral_graph, cardinalities):
        neighbours = list(moral_graph.neighbors(node))
        return reduce(
                (lambda x, y: x * y),
                [cardinalities[n] for n in neighbours])\
            if len(neighbours) > 0 else 1
