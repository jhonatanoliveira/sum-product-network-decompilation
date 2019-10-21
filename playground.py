from models import BayesianNetwork, SubplotDrawer
from compilation import compile_variable_elimination, EliminationOrdering
from convert_ac2spn import convert_ac2spn
from comp_marg_spn import compile_marginalized_spn
from decompilation import decompile


if __name__ == "__main__":

    bn_file_name = "bns/complex_diamond.bn"
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
    comp_marg_subplot_title = "Compiled Marginalized SPN"
    if is_collapse_sums:
        comp_marg_subplot_title += " - Collpase Sum"
    if comp_marg_spn == "none":
        comp_marg_subplot_title += " - None Marg"
    elif comp_marg_spn == "internal":
        comp_marg_subplot_title += " - Marg Internal"

    plotter.add(spn, comp_marg_subplot_title)

    # Decompilation
    decomp_bn = decompile(spn, comp_assumption="ve")
    plotter.add(decomp_bn, "Decompiled BN")

    print("------ Report ------")
    print("-> Ori BN #nodes: " + str(len(bn.dag.nodes())))
    print("-> Decomp BN #nodes: " + str(len(decomp_bn.dag.nodes())))
    print("-> Ori BN #edges: " + str(len(bn.dag.edges())))
    print("-> Decomp BN #edges: " + str(len(decomp_bn.dag.edges())))
    # mor_dag = EliminationOrdering.moral_graph(bn.dag)
    # mor_bn = BayesianNetwork(mor_dag, bn.var_cardinalities, bn.factors)
    # plotter.add(mor_bn, "Moral BN")
    # print("----> Moral amt edges: " + str(len(mor_dag.edges())))

    # Show subplots
    plotter.plot()
