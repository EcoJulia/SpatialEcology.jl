

@enum coordstype auto grid points
#@enum inputdatatype auto phylocom worldmapfile benholtmatrix

abstract OccData
abstract SpatialData
abstract Assmbl <: SpatialData

# I could implement sitestats as a Dict with several DataFrames to make space for big data sets, but I prefer to not do this now. Example below.


type SiteFields
    coords::NamedArrays.NamedMatrix{Float64}  #This should be spatialpoints - not yet implemented?
    cdtype::coordstype
    sitestats::DataFrames.DataFrame
    shape::Nullable{Shapefile.Handle}

    # inner constructor
    function SiteFields(coords, cdtype = auto,
            sitestats = DataFrames.DataFrame(id = 1:size(coords,1)),
            shape = Nullable{Shapefile.Handle}())

        nrow(sitestats) == size(coords, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        if cdtype == auto
            cdtype = isgrid(coords) ? grid : points
        end
        new(coords, cdtype, sitestats, shape)
    end
end

type ComMatrix{T}
    occurrences::NamedArrays.NamedArray{T, 2}
end

type OccFields{T}
    commatrix::ComMatrix{T}
    traits::DataFrames.DataFrame

    function OccFields(commatrix, traits = DataFrames.DataFrame(id = 1:Nspecies(commatrix)))
        nrow(traits) ==  nspecies(commatrix) || throw(DimensionMismatch("Wrong number of species in traits"))
        new(occurrences, traits)
    end
end


type SiteData <: SpatialData
    site::SiteFields
end

type Assemblage{T} <: Assmbl # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::SiteFields
    occ::OccFields{T}

    # inner constructor
    function Assemblage(site::SiteFields, occ::OccFields)
        size(occ.commatrix.occurrences, 1) == size(site.coords, 1) || error("Length mismatch between occurrence matrix and coordinates")
        new(site, occ)
    end
end





#############################################################


# type SiteData <: SpatialData
#     coords::Matrix{Float64}  #This should be spatialpoints - not yet implemented?
#     cdtype::coordstype
#     sitestats::Dict{Symbol, DataFrames.DataFrame}
#     shape::Nullable{Shapefile.Handle}
#
#     # inner constructor
#     function SiteData(coords, cdtype,
#             sitestats = Dict(:sites => DataFrames.DataFrame(site = 1:size(coords,1))),
#             shape = Nullable{Shapefile.Handle}())
#
#         [nrow(v) == size(coords, 1) || error("Wrong number of rows in $k") for (k,v) in sitestats] # a little check for the right number
#         new(coords, cdtype, sitestats, shape)
# end
