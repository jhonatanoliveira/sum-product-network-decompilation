from models import BayesianNetwork, SubplotDrawer
from compilation import compile_variable_elimination, EliminationOrdering
from convert_ac2spn import convert_ac2spn


if __name__ == "__main__":

    plotter = SubplotDrawer("Compiling and Decompiling BNs")

    bn = BayesianNetwork.get_bn_from_file("bns/hmm.bn")
    plotter.add(bn, "BN")

    elim_ord = EliminationOrdering.get_elimination_ordering(bn, "rev")
    print("Elimination ordering: " + str(elim_ord))

    ac = compile_variable_elimination(bn, elim_ord)
    plotter.add(ac, "Compiled AC")

    spn = convert_ac2spn(ac, plotter)

    plotter.plot()
