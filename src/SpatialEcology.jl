__precompile__()
module SpatialEcology

using Reexport
@reexport using DataFrames
@reexport using NamedArrays

#import Bio.Phylo
import RecipesBase
import RCall: @R_str, rcopy
import PlotUtils: @colorant_str, register_gradient_colors, sample_evenly
import Base: copy, getindex, setindex!, size, show, summary, view


export SiteData, ComMatrix, Assemblage, coordtype, DispersionField #types and their constructors
export nspecies, nsites, occupancy, richness, records, sitenames, specnames, coordinates
export copy, setindex!, getindex, size, show, summary, view
export coords, subset
export xcells, ycells, cells, xmin, xmax, ymin, ymax, xrange, yrange, xcellsize, ycellsize, cellsize, boundingbox #it is possible that I will export none of these
export getRobject
export registercolors

include("DataTypes.jl")
include("Constructor_helperfunctions.jl")
include("Constructors.jl")
include("Commatrix_functions.jl")
include("GetandSetdata.jl")
include("Gridfunctions.jl")
include("Subsetting.jl")
include("RObjects.jl")
include("PlotRecipes.jl")
include("Colorgradients.jl")
include("DispersionFields.jl")

end # module
