__precompile__()
module SpatialEcology

using Reexport
@reexport using DataFrames
@reexport using NamedArrays

#import Bio.Phylo
import RecipesBase
import Shapefile
import RCall: @R_str, rcopy
import Base: getindex, setindex!, size, show, summary

include("DataTypes.jl")
include("Constructor_helperfunctions.jl")
include("Constructors.jl")
include("Commatrix_functions.jl")
include("GetandSetdata.jl")
include("Gridfunctions.jl")
include("Subsetting.jl")
include("RObjects.jl")
#include("PlotRecipes.jl")

export SiteData, ComMatrix, Assemblage, coordtype #types and their constructors
export nspecies, nsites, occupancy, richness, records, sitenames, specnames
export setindex!, getindex, size, show, summary
export coords, subset!, subset, addshape!, deleteshape!
export xcells, ycells, cells, xmin, xmax, ymin, ymax, xrange, yrange, xcellsize, ycellsize, cellsize, boundingbox #it is possible that I will export none of these

end # module
