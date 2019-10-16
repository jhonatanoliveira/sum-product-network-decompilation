from models import BayesianNetwork, SubplotDrawer
from compilation import compile_variable_elimination, EliminationOrdering
from convert_ac2spn import convert_ac2spn
from comp_marg_spn import compile_marginalized_spn
from decompilation import decompile


if __name__ == "__main__":

    bn_file_name = "bns/hmm.bn"
    is_collapse_sums = False
    elim_ord_type = "rev"
    comp_marg_spn = "internal"

    plotter = SubplotDrawer("Compiling and Decompiling BNs")

    bn = BayesianNetwork.get_bn_from_file(bn_file_name)
    plotter.add(bn, "BN")

    elim_ord = EliminationOrdering.get_elimination_ordering(bn, elim_ord_type)
    print("Elimination ordering: " + str(elim_ord))

    # Compilation
    ac = compile_variable_elimination(bn, elim_ord)
    plotter.add(ac, "Compiled AC")

    # Convert to SPN
    spn = convert_ac2spn(ac, plotter)

    # Marginalized SPN
    marg_nodes = []
    if comp_marg_spn == "internal":
        marg_nodes = [var for var in bn.dag.nodes()
                      if len(list(bn.dag.successors(var))) > 0]
    spn = compile_marginalized_spn(spn, marg_nodes,
                                   collapse_sums=is_collapse_sums,
                                   subplot_drawer=None)
    comp_marg_subplot_title = "Compiled Marginalized SPN - Collpase Sum"\
        if is_collapse_sums else "Compiled Marginalized SPN"

    plotter.add(spn, comp_marg_subplot_title)

    # Decompilation
    bn = decompile(spn, comp_assumption="ve")
    plotter.add(bn, "Decompiled BN")

    # Show subplots
    plotter.plot()
