

@enum coordstype auto griddata pointdata
#@enum inputdatatype auto phylocom worldmapfile benholtmatrix

abstract OccData
abstract SpatialData
abstract Assmbl <: SpatialData  #Not sure about this structure - so far no type inherits from occdata. Perhaps SimpleTraits.jl is/has a solution
# this is here because we also need phylogeny assemblages
abstract AbstractAssemblage <: Assmbl
abstract AbstractOccFields{T}
abstract AbstractComMatrix{T}
abstract SiteFields

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


abstract AbstractPointData <: SiteFields

type PointData <: AbstractPointData
    coords::NamedArrays.NamedMatrix{Float64}
    sitestats::DataFrames.DataFrame
    # inner constructor
    function PointData(coords, sitestats = DataFrames.DataFrame(id = 1:size(coords,1)))

        DataFrames.nrow(sitestats) == size(coords, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(coords, sitestats)
    end
end

abstract AbstractGridData <: SiteFields

type GridData <: AbstractGridData
    indices::NamedArrays.NamedMatrix{Int}
    grid::GridTopology
    sitestats::DataFrames.DataFrame

    function GridData(indices, grid, sitestats = DataFrames.DataFrame(id = 1:size(coords,1)))

        DataFrames.nrow(sitestats) == size(indices, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(indices, grid, sitestats)
    end
end



type ComMatrix{T <: Union{Bool, Int}} <: AbstractComMatrix{T}
    occurrences::NamedArrays.NamedArray{T, 2} #this is sparse
end

type OccFields{T <: Union{Bool, Int}} <: AbstractOccFields{T}
    commatrix::ComMatrix{T}
    traits::DataFrames.DataFrame

    function OccFields(commatrix, traits)
        DataFrames.nrow(traits) ==  nspecies(commatrix) || throw(DimensionMismatch("Wrong number of species in traits"))
        new(commatrix, traits)
    end
end

abstract AbstractSiteData <: SpatialData
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
        new(site, occ)
    end
end
