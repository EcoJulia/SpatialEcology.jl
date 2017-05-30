
@enum coordstype auto griddata pointdata
#@enum inputdatatype auto phylocom worldmapfile benholtmatrix

abstract type OccData end
abstract type SpatialData end
abstract type Assmbl <: SpatialData  end #Not sure about this structure - so far no type inherits from occdata. Perhaps SimpleTraits.jl is/has a solution
# this is here because we also need phylogeny assemblages
abstract type AbstractAssemblage <: Assmbl end
abstract type AbstractOccFields{T} end
abstract type AbstractComMatrix{T} end
abstract type SiteFields end

# I could implement sitestats as a Dict with several DataFrames to make space for big data sets, but I prefer to not do this now. Example below.

# I could do a lot more with immutable types if I had a clearer view/copy implementation
type GridTopology
    xmin::Number
    xcellsize::Number
    xcells::Int
    ymin::Number
    ycellsize::Number
    ycells::Int
end

type Bbox
    xmin::Number
    xmax::Number
    ymin::Number
    ymax::Number
end

abstract type AbstractPointData <: SiteFields end

# Do I need sitenames here? I think so, they should match those in sitestats, and be separate
type PointData <: AbstractPointData
    coords::Matrix{Float64}
    sitestats::DataFrames.DataFrame
    # inner constructor
    function PointData(coords, sitestats = DataFrames.DataFrame(id = 1:size(coords,1)))

        DataFrames.nrow(sitestats) == size(coords, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(coords, sitestats)
    end
end

abstract type AbstractGridData <: SiteFields end

type GridData <: AbstractGridData
    indices::Matrix{Int}
    grid::GridTopology
    sitestats::DataFrames.DataFrame

    function GridData(indices, grid, sitestats = DataFrames.DataFrame(id = 1:size(coords,1)))

        DataFrames.nrow(sitestats) == size(indices, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(indices, grid, sitestats)
    end
end

type ComMatrix{T <: Union{Bool, Int}} <: AbstractComMatrix{T}
    occurrences::SparseMatrixCSC{T}
    specnames::Vector{String}
    sitenames::Vector{String}
end

# likewise, do I need a specnames here? Should traits have a :series field (like now) or all matching be done on the specnames?
type OccFields{T <: Union{Bool, Int}} <: AbstractOccFields{T}
    commatrix::ComMatrix{T}
    traits::DataFrames.DataFrame

    function OccFields(commatrix, traits)
        DataFrames.nrow(traits) ==  nspecies(commatrix) || throw(DimensionMismatch("Wrong number of species in traits"))
        new(commatrix, traits)
    end
end

abstract type AbstractSiteData <: SpatialData end

# Not really sure what this type is for
type SiteData{S <: SiteFields} <: AbstractSiteData
    site::S
end

type Assemblage{S <: SiteFields, T <: Union{Bool, Int}} <: AbstractAssemblage # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::S
    occ::OccFields{T}

    # inner constructor
    function Assemblage(site, occ)
        size(occ.commatrix.occurrences, 1) == size(coordinates(site), 1) || error("Length mismatch between occurrence matrix and coordinates")
        #TODO activate this # sitenames(occ) == sitenames(site) || error("sitenames do not match") #I need a constructor that matches them up actively
        new(site, occ)
    end
end
