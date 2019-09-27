import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from functools import reduce
import argparse
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
                factor.variables[0] + str(assignment[0]) + ")"
            param_var = "P(" +\
                factor.variables[0] + str(assignment[0]) + "|" +\
                ",".join([factor.variables[k] + str(assignment[k])
                          for k in range(1, factor.total_variables)]) + ")"\
                if factor.total_variables > 1\
                else "P(" + factor.variables[0] + str(assignment[0]) + ")"
            ac_dag.add_node(indc_var)
            ac_dag.add_node(param_var)
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


def draw_graph(graph):
    pos = graphviz_layout(graph, prog='dot')
    nx.draw(graph, pos, with_labels=True, arrows=True)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile a BN into an AC")
    parser.add_argument("--bn", type=str, help="Path to '.bn' file.")

    bn = get_bn_from_file("hmm.bn")

    compile_variable_elimination(bn, ["X3", "H3", "X2", "H2", "X1", "H1"])
    # draw_bn(bn.dag)
