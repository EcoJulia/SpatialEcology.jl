
@enum coordstype auto griddata pointdata

# I could implement sitestats as a Dict with several DataFrames to make space for big data sets, but I prefer to not do this now. Example below.

abstract type AbstractComMatrix{D<:Real} end

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

mutable struct ComMatrix{D} <: AbstractComMatrix{D}
    occurrences::SparseMatrixCSC{D}
    specnames::Vector{String}
    sitenames::Vector{String}
    ComMatrix{T}(occ::SparseMatrixCSC{T}, spn::Vector{String}, sin::Vector{String}) where {T} = new(dropzeros!(occ), spn, sin)
end

# likewise, do I need a specnames here? Should traits have a :series field (like now) or all matching be done on the specnames?
mutable struct Species <: AbstractThings
    speciesnames::Vector{<:Union{String, Symbol}}
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

mutable struct Assemblage{S, T} <: AbstractAssemblage where {S <: SiteFields, T <: OccTypes} # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::S
    occ::OccFields{T}

    # inner constructor
    function Assemblage{D,  P}(commat::ComMatrix{D}, species::Species, site::P) where {D <: Real,  P <: AbstractPlaces}
        size(occurrences(commatrix), 2) == size(coordinates(site), 1) || error("Length mismatch between occurrence matrix and coordinates")
        size(occurrences(commatrix), 1) == nspecies(Species) || error("Length mismatch between occurrence matrix and coordinates")
        new(commat, species, site)
    end
end
