import networkx as nx
import argparse
from main import get_bn_from_file, reconstruct_bn

parser = argparse.ArgumentParser(description="Compile a BN into an AC")
parser.add_argument("--bn", type=str, help="Path to '.bn' file.")
args = parser.parse_args()

if __name__ == "__main__":

    bn = get_bn_from_file(args.bn)
    all_elim_ord = list(nx.algorithms.dag.all_topological_sorts(bn.dag))
    print("Total elimination oderings: " + str(len(all_elim_ord)))
    for elim_ord in all_elim_ord:
        elim_ord = list(reversed(elim_ord))
        rec_bn = reconstruct_bn(bn, elim_ord, plot=True)
