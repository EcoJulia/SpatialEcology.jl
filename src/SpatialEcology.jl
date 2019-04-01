__precompile__(true)
module SpatialEcology

using DataFrames
import DataFrames: aggregate
import DataFramesMeta: @with
using SparseArrays
using Statistics
using RandomBooleanMatrices
import RandomNumbers.Xorshifts: Xoroshiro128Plus
import RandomBooleanMatrices: matrixrandomizer, MatrixGenerator

import EcoBase
import EcoBase: nthings, nplaces, occupancy, richness, nrecords, placenames, thingnames,
        occurring, noccurring, occupied, noccupied, occurrences,
        placeoccurrences, thingoccurrences, cooccurring, places, things,
        indices, coordinates, xcells, ycells, cells, xmin, xmax, ymin, ymax,
        xrange, yrange, xcellsize, ycellsize, cellsize, getcoords

import RecipesBase
import Base: copy, getindex, setindex!, size, show, summary, view, Meta.isexpr
import StatsBase: nquantile

export SiteData, ComMatrix, Assemblage, coordtype #types and their constructors
export AbstractComMatrix
export nspecies, nsites, occupancy, richness, nrecords, sitenames, speciesnames, coordinates
export occurring, noccurring, occupied, noccupied, occurrences, cooccurring
export traits, sitestats, sitestatnames, traitnames, commatrix
export addtraits!, addsitestats!
export asquantiles, asquantiles!
export coordstype, subset
export xcells, ycells, cells, xmin, xmax, ymin, ymax, xrange, yrange, xcellsize, ycellsize, cellsize, boundingbox #it is possible that I will export none of these
export sitetotals, speciestotals, getspecies, getsite
export groupspecies, groupsites
export aggregate
export @with, @traits, @sitestats
export matrixrandomizer, matrixrandomizations
export dispersionfield

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
include("Operations.jl")
include("Randomizations.jl")
include("PlotRecipes.jl")

end # module
