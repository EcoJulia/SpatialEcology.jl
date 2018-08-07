__precompile__(true)
module SpatialEcology

using DataFrames
using SparseArrays
using Statistics

import EcoBase
import RecipesBase
import Base: copy, getindex, setindex!, size, show, summary, view, Meta.isexpr
import StatsBase: nquantile

export SiteData, ComMatrix, Assemblage, coordtype, DispersionField #types and their constructors
export AbstractComMatrix
export nspecies, nsites, occupancy, richness, nrecords, sitenames, specnames, coordinates
export occurring, noccurring, occupied, noccupied, occurrences, cooccurring
export traits, sitestats, sitestatnames, traitnames, commatrix
export addtraits!, addsitestats!
export asquantiles, asquantiles!
export coordstype, subset
export xcells, ycells, cells, xmin, xmax, ymin, ymax, xrange, yrange, xcellsize, ycellsize, cellsize, boundingbox #it is possible that I will export none of these
export sitetotals, speciestotals, getspecies, getsite
export groupspecies, groupsites

include("DataTypes.jl")
include("Constructor_helperfunctions.jl")
include("Constructors.jl")
include("Sparse_matrixfunctions.jl")
include("Conveniencefunctions.jl")
include("ComMatrix.jl")
include("Subsetting.jl")
include("GetandSetdata.jl")
include("Gridfunctions.jl")
include("Grouping.jl")
include("DispersionFields.jl")

end # module
