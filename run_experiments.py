from models import BayesianNetwork, SubplotDrawer
from compilation import compile_variable_elimination, EliminationOrdering
from convert_ac2spn import convert_ac2spn
from comp_marg_spn import compile_marginalized_spn


if __name__ == "__main__":

    plotter = SubplotDrawer("Compiling and Decompiling BNs")

    bn = BayesianNetwork.get_bn_from_file("bns/diamond2.bn")
    plotter.add(bn, "BN")

    elim_ord = EliminationOrdering.get_elimination_ordering(bn, "rev")
    print("Elimination ordering: " + str(elim_ord))

    # Compilation
    ac = compile_variable_elimination(bn, elim_ord)
    plotter.add(ac, "Compiled AC")

    # Convert to SPN
    spn = convert_ac2spn(ac, plotter)

    # Marginalized SPN
    internal_nodes = [var for var in bn.dag.nodes()
                      if len(list(bn.dag.successors(var))) > 0]
    spn = compile_marginalized_spn(spn, internal_nodes,
                                   collapse_sums=False, subplot_drawer=None)
    plotter.add(spn, "Compiled Marginalized SPN")

    # Show subplots
    plotter.plot()
