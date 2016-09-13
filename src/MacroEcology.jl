module MacroEcology

import DataFrames
import Bio.Phylo
import NamedArrays

include("DataTypes.jl")
include("Constructors.jl")
include("PlotRecipes.jl")

export Assmbl, Assemblage, PhyloAssemblage, OccMatrix, AbundanceMatrix, PAMatrix #types and their constructors

end # module
