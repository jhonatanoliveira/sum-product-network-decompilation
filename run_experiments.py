from models import BayesianNetwork
from compilation import compile_variable_elimination, EliminationOrdering
from convert_ac2spn import convert_ac2spn


if __name__ == "__main__":

    bn = BayesianNetwork.get_bn_from_file("bns/hmm.bn")
    elim_ord = EliminationOrdering.get_elimination_ordering(bn, "rev")
    ac = compile_variable_elimination(bn, elim_ord)
    spn = convert_ac2spn(ac)
    spn.draw()
