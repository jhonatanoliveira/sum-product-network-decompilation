from region_graph import RegionGraph
from main import ArithmeticCircuit, get_bn, spn_assign_lv_sum_depth,\
    draw_subplot_graphs
import itertools
import argparse
import networkx as nx


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Compile P&D into an AC and decompile to BN")
    parser.add_argument("--size", type=int, help="Consider image size x size")
    args = parser.parse_args()

    size = args.size
    height, width, delta = size, size, 1
    input_leaves_amt = 2
    sum_nodes_amt = 2

    pd = RegionGraph(list(range(height * width)))
    pd.make_poon_structure(width, height, 1)
    layers = pd.make_layers()

    dag = nx.DiGraph()

    input_layer = layers[0]
    distributions = {}
    leaf_counter = 0
    for input_region in input_layer:
        distributions[input_region] = []
        for idx in range(input_leaves_amt):
            leaf_name = "T(" + str(leaf_counter) + "-" + str(idx) + ")"
            distributions[input_region].append(leaf_name)
            dag.add_node(leaf_name)
        leaf_counter += 1

    is_partition = True
    prod_counter = 0
    sum_counter = 0
    for layer_idx, layer in enumerate(layers[1:]):
        layer_idx += 1
        if is_partition:
            is_partition = False
            for partition in layer:
                distributions[partition] = []
                sum_nodes = [distributions[region] for region in partition]
                for prod_children in itertools.product(*sum_nodes):
                    prod_children = list(prod_children)
                    prod_name = "*" + str(prod_counter)
                    prod_counter += 1
                    distributions[partition].append(prod_name)
                    dag.add_node(prod_name)
                    dag.add_edges_from(
                        list(itertools.product([prod_name], prod_children)))
        else:
            is_partition = True
            if len(layer) > 1:
                for region in layer:
                    distributions[region] = []
                    for idx_amt in range(sum_nodes_amt):
                        sum_name = "+" + str(sum_counter)
                        sum_counter += 1
                        distributions[region].append(sum_name)
                        dag.add_node(sum_name)
                        for partition in layers[layer_idx - 1]:
                            rec_region = set().union(*partition)
                            if rec_region == set(region):
                                sum_children = distributions[partition]
                                dag.add_edges_from(
                                    list(itertools.product(
                                        [sum_name], sum_children)))
            else:
                region = layer[0]
                sum_name = "+" + str(sum_counter)
                sum_counter += 1
                distributions[region] = [sum_name]
                dag.add_node(sum_name)
                for partition in layers[layer_idx - 1]:
                    sum_children = distributions[partition]
                    dag.add_edges_from(list(itertools.product(
                        [sum_name], sum_children)))

    ac = ArithmeticCircuit(dag)
    bn = get_bn(ac, spn_assign_lv_sum_depth)

    graphs = [ac, bn]
    graphs_subtitles = ["P&D SPN", "Decompiled BN"]
    main_title = "Decompiling P&D SPN"
    draw_subplot_graphs(graphs, graphs_subtitles, main_title)
