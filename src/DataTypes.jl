

@enum coordstype auto grid points
@enum inputdatatype auto phylocom worldmapfile benholtmatrix

abstract SpatialData
abstract Assmbl <: SpatialData

# I could implement sitestats as a Dict with several DataFrames to make space for big data sets, but I prefer to not do this now. Example below.

type SiteFields
    coords::NamedArrays.NamedMatrix{Float64}  #This should be spatialpoints - not yet implemented?
    cdtype::coordstype
    sitestats::DataFrames.DataFrame
    shape::Nullable{ShapeFiles.ShapeFile}

    # inner constructor
    function SiteFields(coords, cdtype,
            sitestats = DataFrames.DataFrame(site = 1:size(coords,1)),
            shape = Nullable{ShapeFiles.ShapeFile}())

        nrow(sitestats) == size(coords, 1) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(coords, cdtype, sitestats, shape)
    end
end

type OccFields{T}
    occurrences::NamedArrays.NamedArray{T, 2}
    traits::DataFrames.DataFrame

    function OccFields(occurrences, traits)
        nrow(traits) == size(occurrences, 2) || throw(DimensionMismatch("Wrong number of species in traits"))
        new(occurrences, traits)
    end
end

type PhyloFields
    phylo::Phylo.Phylogeny #This makes it a special type, because so many functions are only defined when there is a phylogeny
    nodespecies::Matrix{Bool}

    function PhyloFields(phylo, nodespecies)
        (Ntip(phylo) == size(nodespecies, 1) && Nnode(phylo) == size(nodespecies, 2)) || throw(DimensionMismatch("Dimension mismatch between nodespecies matrix and phylogeny"))

type SiteData <: SpatialData
    site::SiteFields
end

type Assemblage{T} <: Assmbl # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::SiteFields
    occ::OccFields{T}

    # inner constructor
    function Assemblage(site, occ)
        size(occ.occurrences, 1) == size(site.coords, 1) || error("Length mismatch between occurrence matrix and coordinates")
        new(site, occ)
    end
end

type PhyloAssemblage{T} <: Assmbl # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::SiteFields
    occ::OccFields{T}
    phy::PhyloFields

    function PhyloAssemblage(site, occ, phy)
        size(occ.occurrences, 1) == size(site.coords, 1) || throw(DimensionMismatch("Length mismatch between occurrence matrix and coordinates"))
        Ntip(phy.Phylo) == size(occ.occurrences, 2) || throw(DimensionMismatch("Occurrence matrix and phylogeny do not match in species numbers"))
        new(site, occ, phy)
    end
end





#############################################################


# type SiteData <: SpatialData
#     coords::Matrix{Float64}  #This should be spatialpoints - not yet implemented?
#     cdtype::coordstype
#     sitestats::Dict{Symbol, DataFrames.DataFrame}
#     shape::Nullable{ShapeFiles.ShapeFile}
#
#     # inner constructor
#     function SiteData(coords, cdtype,
#             sitestats = Dict(:sites => DataFrames.DataFrame(site = 1:size(coords,1))),
#             shape = Nullable{ShapeFiles.ShapeFile}())
#
#         [nrow(v) == size(coords, 1) || error("Wrong number of rows in $k") for (k,v) in sitestats] # a little check for the right number
#         new(coords, cdtype, sitestats, shape)
# end
