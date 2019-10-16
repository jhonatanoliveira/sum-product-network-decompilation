import math
from functools import reduce
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt  # noqa: E402


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

    @staticmethod
    def construct_factors(bn_dag, var_cardinalities):
        factors = {}
        for node in bn_dag.nodes():
            variables = [node] + list(bn_dag.predecessors(node))
            cardinalities = {var: var_cardinalities[var] for var in variables}
            total_values = reduce(
                (lambda x, y: x * y),
                [cardinalities[k] for k in cardinalities])
            values = [0 for _ in range(0, total_values)]
            factor = Factor(variables, cardinalities, values)
            factors[node] = factor
        return factors

    def stride(self, variable):
        return self._strides[variable] if variable in self._strides else 0


class ProbabilisticGraphicalModel:

    def draw(self, show=True):
        raise NotImplementedError

    def draw_graph(self, graph, show, color_map=None, labels=None):
        pos = graphviz_layout(graph, prog='dot')
        nx.draw(graph, pos, with_labels=True, arrows=True,
                node_color=color_map, node_size=400,
                font_size=8, labels=labels)
        if show:
            plt.show()

    def copy(self):
        raise NotImplementedError


class BayesianNetwork(ProbabilisticGraphicalModel):

    def __init__(self, dag, var_cardinalities, factors):
        # nx Graph
        self.dag = dag
        # Dict {'A': 2, 'B': 3}
        self.variables = list(var_cardinalities.keys())
        self.var_cardinalities = var_cardinalities
        self.factors = factors

    @staticmethod
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
        factors = Factor.construct_factors(bn_dag, var_cardinalities)
        return BayesianNetwork(bn_dag, var_cardinalities, factors)

    def draw(self, show=True):
        NODE_COLOR = "#DA667B"
        color_map = []
        labels = {}
        for node in self.dag.nodes():
            color_map.append(NODE_COLOR)
            labels[node] = node
        self.draw_graph(self.dag, show, color_map, labels)


class ArithmeticCircuit(ProbabilisticGraphicalModel):

    def __init__(self, dag, var_cardinalities):
        self.dag = dag
        self.variables = list(var_cardinalities.keys())
        self.var_cardinalities = var_cardinalities

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
        self.draw_graph(self.dag, show, color_map, labels)


class SumProductNetwork(ProbabilisticGraphicalModel):

    def __init__(self, dag, var_cardinalities):
        self.dag = dag
        self.variables = list(var_cardinalities.keys())
        self.var_cardinalities = var_cardinalities

    def scopes(self):
        # compute node scopes
        scopes = {}
        for node in reversed(list(nx.topological_sort(self.dag))):
            if "+" in node or "*" in node:
                scope = set()
                for child in self.dag.successors(node):
                    scope = scope.union(scopes[child])
                scopes[node] = scope
            elif "T" in node:
                tmp = node.replace("T(", "").replace(")", "")
                scope = set([tmp[:tmp.index("-")]])
                scopes[node] = scope
            else:
                raise NotImplementedError("Node scope not defined")
        return scopes

    def draw(self, show=True):
        SUM_COLOR, PROD_COLOR, IND_COLOR, PROB_COLOR, TERM_COLOR =\
            "#DEE5E5", "#009DDC", "#F26430", "#6761A8", "#009B72"
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
        self.draw_graph(self.dag, show, color_map, labels)

    def copy(self):
        return SumProductNetwork(self.dag.copy(), self.var_cardinalities)


class SubplotDrawer:

    def __init__(self, main_title=None):
        self.graphs = []
        self.subtitles = []
        self.main_title = main_title if main_title else "Graphs"

    def add(self, graph, subtitle=None):
        self.graphs.append(graph)
        _sub = subtitle if subtitle else "Graph " + str(len(self.graphs))
        self.subtitles.append(_sub)

    def plot(self):
        total_graphs = len(self.graphs)
        subplot_amt = math.ceil(math.sqrt(total_graphs))
        subplot_rows = (subplot_amt - 1)\
            if (subplot_amt * subplot_amt - total_graphs) >= subplot_amt\
            else subplot_amt
        plt.figure()
        for i, graph in enumerate(self.graphs):
            plt.subplot(subplot_rows, subplot_amt, i + 1)
            graph.draw(show=False)
            if self.subtitles:
                plt.title(self.subtitles[i])
            plt.grid(True)
        if self.main_title:
            plt.suptitle(self.main_title)
        plt.show()
