__precompile__()
module MacroEcology

import DataFrames
import Bio.Phylo
import NamedArrays
import Shapefile
import Base: getindex, setindex!, size, show, summary

include("DataTypes.jl")
include("Constructor_helperfunctions.jl")
include("Constructors.jl")
include("Commatrix_functions.jl")
include("Gettersandsetters.jl")
include("Subsetting.jl")
include("PlotRecipes.jl")

export SiteData, ComMatrix, Assemblage, coordtype #types and their constructors
export nspecies, nsites, occupancy, richness, records, sitenames, specnames
export setindex!, getindex, size, show, summary
export coords, subset!, subset, addshape!, deleteshape!
export DataFrames, NamedArrays

end # module
