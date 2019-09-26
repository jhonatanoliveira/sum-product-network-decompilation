import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from functools import reduce
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
        total_vals = reduce(
            (lambda x, y: x * y), [cardinalities[k] for k in cardinalities])
        values = [0 for _ in range(0, total_vals)]
        factor = Factor(variables, cardinalities, values)
        factors[node] = factor
    return BayesianNetwork(bn_dag, var_cardinalities, factors)


def ac_leaves(ac_graph, bn, node_counters):
    for variable in bn.factors:
        factor = bn.factors[variable]
        # j = 0
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
            ac_graph.add_node(indc_var)
            ac_graph.add_node(param_var)
            n_prod = "*" + str(node_counters[1])
            node_counters[1] += 1  # increment product counter
            factor.values[i] = n_prod
            ac_graph.add_node(n_prod)
            ac_graph.add_edge(n_prod, indc_var)
            ac_graph.add_edge(n_prod, param_var)
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
    draw_graph(ac_dag)
    # for variable in elim_ord:
    #     var_factors = list(
    #         filter(lambda factor: variable in factor, curr_bn_factors))
    #     prod_factor = None
    #     for factor in var_factors:
    #         if prod_factor is None:
    #             prod_factor = factor
    #         else:
    #             pass
    #             prod_factor = ac_factor_product(
    #                 prod_factor, factor,
    #                 bn.var_cardinalities)


def draw_graph(graph):
    pos = graphviz_layout(graph, prog='dot')
    nx.draw(graph, pos, with_labels=True, arrows=True)
    plt.show()


if __name__ == "__main__":
    bn = get_bn_from_file("toy.bn")
    # draw_bn(bn.dag)
    compile_variable_elimination(bn, ["A"])
