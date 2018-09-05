
@enum coordstype auto griddata pointdata

const OccTypes = Union{Bool, Int, Float64}

abstract type OccData end
abstract type SpatialData end
abstract type Assmbl <: SpatialData  end #Not sure about this structure - so far no type inherits from occdata. Perhaps SimpleTraits.jl is/has a solution
# this is here because we also need phylogeny assemblages
abstract type AbstractAssemblage <: Assmbl end
abstract type AbstractOccFields{T<:OccTypes} end
abstract type AbstractComMatrix{T<:OccTypes} end
abstract type SiteFields end

# I could implement sitestats as a Dict with several DataFrames to make space for big data sets, but I prefer to not do this now. Example below.

# I could do a lot more with immutable types if I had a clearer view/copy implementation
mutable struct GridTopology
    xmin::Number
    xcellsize::Number
    xcells::Int
    ymin::Number
    ycellsize::Number
    ycells::Int
end

mutable struct Bbox
    xmin::Number
    xmax::Number
    ymin::Number
    ymax::Number
end

abstract type AbstractPointData <: EcoBase.AbstractLocations end

# Do I need sitenames here? I think so, they should match those in sitestats, and be separate
mutable struct PointData <: AbstractPointData
    coords::Matrix{Float64}
    sitestats::DataFrames.DataFrame
    # inner constructor
    function PointData(coords, sitestats = DataFrames.DataFrame(id = 1:size(coords,1)))

        DataFrames.nrow(sitestats) == size(coords, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(coords, sitestats)
    end
end

mutable struct GridData <: EcoBase.AbstractGrid
    indices::Matrix{Int}
    grid::GridTopology
    sitestats::DataFrames.DataFrame

    function GridData(indices, grid, sitestats = DataFrames.DataFrame(id = 1:size(coords,1)))

        DataFrames.nrow(sitestats) == size(indices, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(indices, grid, sitestats)
    end
end

mutable struct ComMatrix{T} <: AbstractComMatrix{T}
    occurrences::SparseMatrixCSC{T}
    specnames::Vector{String}
    sitenames::Vector{String}
    ComMatrix{T}(occ::SparseMatrixCSC{T}, spn::Vector{String}, sin::Vector{String}) where {T} = new(dropzeros!(occ), spn, sin)
end

# likewise, do I need a specnames here? Should traits have a :series field (like now) or all matching be done on the specnames?
mutable struct OccFields{T <: OccTypes} <: AbstractOccFields{T}
    commatrix::ComMatrix{T}
    traits::DataFrames.DataFrame
    function OccFields{T}(commatrix::ComMatrix{T}, traits::DataFrames.DataFrame) where T <: OccTypes
        DataFrames.nrow(traits) ==  nspecies(commatrix) || throw(DimensionMismatch("Wrong number of species in traits"))
        new(commatrix, traits)
    end
end


abstract type AbstractSiteData <: SpatialData end

# Not really sure what this type is for
mutable struct SiteData{S} <: AbstractSiteData where S <: SiteFields
    site::S
end

abstract type SEAssemblage{D<:Real, P<:SiteFields} <: EcoBase.AbstractAssemblage{D, OccFields, P} end
    site::S
mutable struct Assemblage{D<:Real, P<:SiteFields} <: SEAssemblage{D, P} # A type to keep subtypes together, ensuring that they are all aligned at all times

    # inner constructor
    function Assemblage{S, T}(site::S, occ::OccFields{T}) where {S <: SiteFields, T <: OccTypes}
        size(occ.commatrix.occurrences, 2) == size(coordinates(site), 1) || error("Length mismatch between occurrence matrix and coordinates")
        #TODO activate this # sitenames(occ) == sitenames(site) || error("sitenames do not match") #I need a constructor that matches them up actively
        new(site, occ)
    end
end
