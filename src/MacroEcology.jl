module MacroEcology

import DataFrames
import Bio.Phylo
import NamedArrays
import Shapefile
import Base: getindex, setindex!, size

include("DataTypes.jl")
include("Constructor_helperfunctions.jl")
include("Constructors.jl")
include("Commatrix_functions.jl")
include("PlotRecipes.jl")

export SiteData, ComMatrix, Assemblage #types and their constructors
export nspecies, nsites, occupancy, richness

end # module
