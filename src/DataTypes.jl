
@enum coordstype auto griddata pointdata

# Insertion point for defining my own functions. They are mainly abstract
# over the full type and it's view, so another possibility would be to
# implement with unions
abstract type SELocations{T<:EcoBase.AbstractLocationData} <: EcoBase.AbstractPlaces{T} end
abstract type SEThings{D <: Real} <: EcoBase.AbstractThings end
abstract type SEGrid <: EcoBase.AbstractRegularGrid end
abstract type SEPoints <: EcoBase.AbstractPoints end


# I could implement sitestats as a Dict with several DataFrames to make space for big data sets, but I prefer to not do this now. Example below.

# I could do a lot more with immutable types if I had a clearer view/copy implementation
mutable struct GridTopology <: EcoBase.AbstractRegularGrid
    xs::StepRangeLen{Float64}
    ys::StepRangeLen{Float64}
end

mutable struct Bbox
    xmin::Number
    xmax::Number
    ymin::Number
    ymax::Number
end

# Do I need sitenames here? I think so, they should match those in sitestats, and be separate
mutable struct PointData <: SEPoints
    coords::Matrix{Float64}
end

mutable struct GridData <: SEGrid
    indices::Matrix{Int}
    grid::GridTopology
end

# RasterData: a SEGrid backed by a Rasters.jl Bool mask raster.
#
# `mask` is the rectangular domain (dims, CRS, extent) — the raster analogue of
# GridData's `grid::GridTopology`. `cellinds` is the per-site index vector — one
# CartesianIndex per site, the analogue of GridData's `indices` matrix. For a
# full RasterData, `cellinds == findall(mask)`; for a site-subset view
# (`SubRasterData`) it is a SubArray over a selection of those cells. This is
# what lets views select an arbitrary, possibly reordered subset of sites — a
# Bool mask alone cannot represent that.
#
# INVARIANT: `mask` is always stored in canonical `(X, Y)` dimension order, so a
# cell's `CartesianIndex` reads as `(x_index, y_index)`. The constructors in the
# Rasters extension enforce this by `permutedims`-ing inputs before building a
# RasterData; downstream code (`coordinates`, `xcells`, …) relies on it.
#
# The `mask` field is left untyped-over-`M` rather than hardcoding `Raster{Bool}`
# so the core package defines this type without Rasters loaded; the bridging
# constructors and methods live in `ext/RastersExt.jl` (loaded with Rasters).
mutable struct RasterData{M} <: SEGrid
    mask::M
    cellinds::Vector{CartesianIndex{2}}
end
RasterData(mask) = RasterData(mask, findall(mask))

mutable struct Locations{T<:Union{GridData, PointData, RasterData}} <: SELocations{T}
    coords::T
    sitestats::DataFrames.DataFrame
    function Locations{T}(coords, sitestats = DataFrames.DataFrame(sites = string.(1:nsites(coords)))) where T
        DataFrames.nrow(sitestats) == nsites(coords) || throw(DimensionMismatch("Wrong number of rows in sitestat")) # a little check for the right number
        new(coords, sitestats)
    end
end

Locations(coords, sitestats) = Locations{typeof(coords)}(coords, sitestats)

#it's a question how many of these structs need to be mutable, as opposed to when you want to allocate a new object
"""
    Species, Sites

The two dimensions of a `ComMatrix`. Their lookups hold the species and site names,
so `com[Species = At("sp_a")]`, `com[At("sp_a"), At("site_1")]` and
`dims(com, Sites)` work as for any DimensionalData array.
"""
DD.@dim Species "Species"
DD.@dim Sites "Sites"

# Names are stored as String whatever AbstractString they arrive as (CSV.jl reads
# DataStrings), as EcoBase's thingnames/placenames promise. The lookup is always
# `Unordered`, even when the names happen to be sorted, so that selecting by name
# never relies on an order the names do not promise.
_namelookup(names) = DD.Lookups.Categorical(String.(names); order = DD.Lookups.Unordered())

"""
    ComMatrix

A species × sites matrix of presences or abundances, with the names of the species
and sites. It is a DimensionalData `AbstractDimArray` with dimensions `Species` and
`Sites`, so it indexes like a matrix: `com[i, j]` is a value, `com[1:3, :]` a
`ComMatrix`, `com[1:3, 2]` a `DimVector` (the dropped site is kept in its `refdims`),
and names are selected with `At`, e.g. `com[At("sp_a"), :]`. A `view` of a
`ComMatrix` is a `SubComMatrix`, which shares its data.

The occurrences are stored in an `OccMatrix`, which keeps their transpose for fast
per-species queries and stays consistent under writes (`com[i, j] = x`, `.=`, and
writes through views).
"""
struct ComMatrix{D<:Real, Dm<:Tuple, R<:Tuple, A<:AbstractMatrix{D}, Na, Me} <: DD.AbstractDimArray{D, 2, Dm, A}
    data::A
    dims::Dm
    refdims::R
    name::Na
    metadata::Me
    function ComMatrix(data::AbstractMatrix, dims::Tuple, refdims::Tuple, name, metadata)
        s = _comstorage(data)
        new{eltype(s), typeof(dims), typeof(refdims), typeof(s), typeof(name), typeof(metadata)}(s, dims, refdims, name, metadata)
    end
end

# A fresh sparse matrix (from slicing, broadcasting, copying, or a constructor) is
# taken over as OccMatrix storage, without a copy; its stored zeros are dropped in
# place. An OccMatrix is shared, as DimensionalData's `rebuild` shares data; a view
# of one is a SubComMatrix; anything else (e.g. dense results of `sum(com, dims = 2)`)
# is kept as it is.
_comstorage(d::SparseMatrixCSC{<:Any, Int}) = OccMatrix(d)
_comstorage(d::SparseMatrixCSC{T}) where T = OccMatrix(SparseMatrixCSC{T, Int}(d))
_comstorage(d::AbstractMatrix) = d

"""
    ComMatrix{D}(occurrences::AbstractSparseMatrixCSC{D}, speciesnames, sitenames)

Build a `ComMatrix` from a species × sites sparse matrix and the names. A
`SparseMatrixCSC{D, Int}` is taken over without a copy (its stored zeros are dropped
in place); any other sparse matrix, including an `OccMatrix`, is copied.
"""
function ComMatrix{D}(occ::SparseArrays.AbstractSparseMatrixCSC{D}, spn::AbstractVector{<:AbstractString}, sin::AbstractVector{<:AbstractString}) where {D}
    length(spn) == size(occ, 1) || throw(DimensionMismatch("$(length(spn)) species names for $(size(occ, 1)) species"))
    length(sin) == size(occ, 2) || throw(DimensionMismatch("$(length(sin)) site names for $(size(occ, 2)) sites"))
    data = occ isa SparseMatrixCSC ? occ : SparseMatrixCSC{D, Int}(occ)
    com = ComMatrix(data, (Species(_namelookup(spn)), Sites(_namelookup(sin))), (), DD.NoName(), DD.NoMetadata())
    transposed(parent(com))   # build the transpose up front, so the first species query is fast
    com
end

# DimensionalData rebuilds arrays after slicing, broadcasting, `copy`, `similar` etc.,
# positionally through `rebuild` or by keyword through ConstructionBase. Either way a
# species × sites result is a ComMatrix again, and anything else (a single species or
# site, a transpose) a plain DimArray.
_iscomdims(dims) = length(dims) == 2 && dims[1] isa Species && dims[2] isa Sites
_comorarray(data, dims, refdims, name, metadata) =
    _iscomdims(dims) && data isa AbstractMatrix ?
        ComMatrix(data, dims, refdims, name, metadata) :
        DD.DimArray(data, dims; refdims, name, metadata)
DD.rebuild(A::ComMatrix, data, dims::Tuple, refdims, name, metadata) =
    _comorarray(data, dims, refdims, name, metadata)
ConstructionBase.constructorof(::Type{<:ComMatrix}) = _comorarray

# DimensionalData turns a Bool mask into a `Base.LogicalIndex` before indexing the
# parent array, but SparseArrays' getindex indexes into its index vectors by
# position, which a LogicalIndex does not support. Hand the storage positions.
_positions(I::Base.LogicalIndex) = collect(I)
_positions(I) = I
Base.@propagate_inbounds function DD.rebuildsliced(f::Function, A::ComMatrix, I::Tuple, name = DD.name(A))
    I1 = map(_positions, to_indices(A, I))
    DD.rebuildsliced(f, A, _sliceparent(f, parent(A), I1), I1, name)
end
# Indexing a view of the storage would go element by element; index the storage
# itself at the view's positions instead, which is a sparse getindex.
_sliceparent(::typeof(getindex), v::SubArray{<:Any, 2, <:OccMatrix}, I) =
    parent(v)[Base.reindex(parentindices(v), I)...]
_sliceparent(f, data, I) = f(data, I...)

"""
    SubComMatrix

A `view` of a `ComMatrix`: a `ComMatrix` whose data is a `SubArray` of its parent's
storage. Writes through it reach the parent.
"""
const SubComMatrix{D<:Real} = ComMatrix{D, <:Tuple, <:Tuple, <:SubArray}

"""
    AbstractComMatrix

Any `ComMatrix`, including views (`SubComMatrix`).
"""
const AbstractComMatrix{D<:Real} = ComMatrix{D}

# A ComMatrix holding its own OccMatrix storage (not a view), which keeps the
# transpose used by the fast per-species methods.
const StoredComMatrix{D<:Real} = ComMatrix{D, <:Tuple, <:Tuple, <:OccMatrix}

# likewise, do I need a speciesnames here? Should traits have a :series field (like now) or all matching be done on the speciesnames?
mutable struct SpeciesData{D} <: SEThings{D}
    commatrix::ComMatrix{D}
    traits::DataFrames.DataFrame
    function SpeciesData{D}(commatrix::ComMatrix{D}, traits::DataFrames.DataFrame) where D <: Real
        DataFrames.nrow(traits) ==  nspecies(commatrix) || throw(DimensionMismatch("Wrong number of species in traits"))
        new(commatrix, traits)
    end
end
SpeciesData(commatrix::ComMatrix{D}, traits) where D<:Real = SpeciesData{D}(commatrix, traits)

abstract type SEAssemblage{D<:Real, T<:SEThings, P<:SELocations} <: EcoBase.AbstractAssemblage{D, T, P} end

"""
    Assemblage{D<:Real, P<:Locations}
"""
mutable struct Assemblage{D<:Real, P<:Locations} <: SEAssemblage{D, SpeciesData{D}, P} # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::P
    occ::SpeciesData{D}

    # inner constructor
    function Assemblage{D, P}(site::P, occ::SpeciesData{D}) where {P <: SELocations, D <: Real}
        size(occurrences(occ), 2) == size(coordinates(site), 1) || error("Length mismatch between occurrence matrix and coordinates")
        #TODO activate this # sitenames(occ) == sitenames(site) || error("sitenames do not match") #I need a constructor that matches them up actively
        new(site, occ)
    end
end
