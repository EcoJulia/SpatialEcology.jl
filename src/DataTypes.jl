

@enum coordstype auto grid points
#@enum inputdatatype auto phylocom worldmapfile benholtmatrix

abstract OccData
abstract SpatialData
abstract Assmbl <: SpatialData  #Not sure about this structure - so far no type inherits from occdata. Perhaps SimpleTraits.jl is/has a solution

# I could implement sitestats as a Dict with several DataFrames to make space for big data sets, but I prefer to not do this now. Example below.



type GridTopology
    xmin::Number
    xcellsize::Number
    xcells::Int
    ymin::Number
    ycellsize::Number
    ycells::Int
end


abstract SiteFields

type PointData <: SiteFields
    coords::NamedArrays.NamedMatrix{Float64}
    sitestats::DataFrames.DataFrame
    shape::Nullable{Shapefile.Handle}
    # inner constructor
    function PointData(coords, sitestats = DataFrames.DataFrame(id = 1:size(coords,1)), shape = Nullable{Shapefile.Handle}())

        DataFrames.nrow(sitestats) == size(coords, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(coords, sitestats, shape)
    end
end

type GridData <: SiteFields
    indices::NamedArrays.NamedMatrix{Int}
    grid::GridTopology
    sitestats::DataFrames.DataFrame
    shape::Nullable{Shapefile.Handle}

    function GridData(indices, grid, sitestats = DataFrames.DataFrame(id = 1:size(coords,1)), shape = Nullable{Shapefile.Handle}())

        DataFrames.nrow(sitestats) == size(indices, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(indices, grid, sitestats, shape)
    end
end




type ComMatrix{T <: Union{Bool, Int}}
    occurrences::NamedArrays.NamedArray{T, 2}
end

type OccFields{T <: Union{Bool, Int}}
    commatrix::ComMatrix{T}
    traits::DataFrames.DataFrame

    function OccFields(commatrix, traits)
        DataFrames.nrow(traits) ==  nspecies(commatrix) || throw(DimensionMismatch("Wrong number of species in traits"))
        new(commatrix, traits)
    end
end



type SiteData{S <: SiteFields} <: SpatialData
    site::S
end


type Assemblage{S <: SiteFields, T <: Union{Bool, Int}} <: Assmbl # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::S
    occ::OccFields{T}

    # inner constructor
    function Assemblage(site, occ)
        size(occ.commatrix.occurrences, 1) == size(site.coords, 1) || error("Length mismatch between occurrence matrix and coordinates")
        new(site, occ)
    end
end

type Bbox
    xmin::Number
    xmax::Number
    ymin::Number
    ymax::Number
end
