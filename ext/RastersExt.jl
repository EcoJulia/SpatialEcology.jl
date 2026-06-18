# RastersExt.jl — the Rasters.jl integration for SpatialEcology.
#
# Loaded automatically when a session has both SpatialEcology and Rasters
# imported. Everything here needs Rasters symbols (Raster, RasterSeries,
# RasterStack, lookup, dims, aggregate, …); the RasterData / SubRasterData
# *types* live in core (parameterized over the mask type) so the rest of the
# package builds without Rasters.
module RastersExt

using SpatialEcology
using Rasters
using DataFrames
using SparseArrays
import EcoBase

using SpatialEcology: RasterData, SubRasterData, AnyRasterData, ComMatrix, SpeciesData,
    Locations, SubLocations, Assemblage, SEAssemblage, Bbox,
    getcoords, places, occurrences, speciesnames, nspecies, nplaces, traits

const SE = SpatialEcology

# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------

# Reorder a raster to canonical (X, Y) dimension order — the RasterData
# invariant (see DataTypes.jl). Errors if the raster is not 2-D over X and Y.
_canon(r::AbstractRaster) = permutedims(r, (X, Y))

# A Bool view of the mask in canonical order, for use as a boolean index.
_boolmask(mask::AbstractRaster) = Bool.(_canon(mask))

# Default domain when the caller gives no mask: the cells that carry valid
# (non-missing) data in *every* species layer. For the common case of plain
# Bool rasters with no missing values this is the whole grid.
_default_mask(layers) = reduce((a, b) -> a .& b, (boolmask(l) for l in layers))

# ----------------------------------------------------------------------------
# grid interface (the methods removed from Gridfunctions.jl)
# ----------------------------------------------------------------------------

SE.xrange(rd::AnyRasterData) = parent(lookup(rd.mask, X()))
SE.yrange(rd::AnyRasterData) = parent(lookup(rd.mask, Y()))
SE.xmin(rd::AnyRasterData)      = minimum(SE.xrange(rd))
SE.xmax(rd::AnyRasterData)      = maximum(SE.xrange(rd))
SE.ymin(rd::AnyRasterData)      = minimum(SE.yrange(rd))
SE.ymax(rd::AnyRasterData)      = maximum(SE.yrange(rd))
SE.xcellsize(rd::AnyRasterData) = abs(step(lookup(rd.mask, X())))
SE.ycellsize(rd::AnyRasterData) = abs(step(lookup(rd.mask, Y())))
SE.xcells(rd::AnyRasterData)    = size(rd.mask, 1)   # dim 1 is X by the invariant
SE.ycells(rd::AnyRasterData)    = size(rd.mask, 2)
SE.boundingbox(rd::AnyRasterData) = Bbox(SE.xmin(rd), SE.xmax(rd), SE.ymin(rd), SE.ymax(rd))

# Location-level forwards (these replace the `@forward_func` lines for the
# raster types that previously sat in Gridfunctions.jl).
for f in (:xmin, :xmax, :ymin, :ymax, :xcellsize, :ycellsize, :xcells, :ycells, :boundingbox)
    @eval SE.$f(l::Locations{<:RasterData})       = SE.$f(l.coords)
    @eval SE.$f(l::SubLocations{<:SubRasterData}) = SE.$f(l.coords)
end

Base.show(io::IO, rd::RasterData) =
    println(io, "RasterData $(SE.xcells(rd))×$(SE.ycells(rd)) grid, $(nplaces(rd)) sites")
Base.show(io::IO, rd::SubRasterData) =
    println(io, "SubRasterData view: $(nplaces(rd)) of $(count(rd.mask)) sites on a $(SE.xcells(rd))×$(SE.ycells(rd)) grid")

# ----------------------------------------------------------------------------
# constructors:  rasters  ->  Assemblage
# ----------------------------------------------------------------------------

# Shared builder. `layers` is a vector of per-species Bool rasters, `sp_names`
# their names. Both the mask and every layer are normalized to canonical (X, Y)
# order, and each layer is validated to share the mask's dims/extent.
function _assemblage(layers, sp_names, mask;
                     sitestats::DataFrame, traits::DataFrame)
    maskc  = _boolmask(mask)
    n_site = count(maskc)

    cols = map(layers) do r
        rc = _canon(r)
        dims(rc) == dims(maskc) ||
            throw(DimensionMismatch(
                "a species raster's dims/extent do not match the mask raster"))
        SparseArrays.sparse(vec(Bool.(rc[maskc])))
    end
    occ = SparseArrays.sparse(reduce(hcat, cols)')   # n_species × n_sites

    site_ids = string.(1:n_site)
    commat   = ComMatrix{Bool}(occ, string.(sp_names), site_ids)
    spdata   = SpeciesData(commat,
                   isempty(traits) ? DataFrame(name = string.(sp_names)) : traits)

    ss = isempty(sitestats) ? DataFrame() : copy(sitestats)
    "sites" in names(ss) || (ss[!, :sites] = site_ids)

    Assemblage(Locations(RasterData(maskc), ss), spdata)
end

"""
    Assemblage(series::RasterSeries, [mask]; sitestats = DataFrame(), traits = DataFrame())
    Assemblage(stack::RasterStack,  [mask]; sitestats = DataFrame(), traits = DataFrame())

Build a presence-absence `Assemblage` from Rasters.jl data, one Bool raster per
species (a `RasterSeries`, or the layers of a `RasterStack`). Species names are
taken from the series index / stack layer names.

`mask` is a Bool raster marking which cells are sites. If omitted, it defaults to
the cells that hold valid (non-missing) data in every species layer. Inputs in
either `(X, Y)` or `(Y, X)` dimension order are accepted and normalized.

Per-site data (e.g. environmental covariates) can be passed as `sitestats`, and
species traits as `traits`.
"""
function SE.Assemblage(series::RasterSeries, mask = nothing;
                       sitestats::DataFrame = DataFrame(), traits::DataFrame = DataFrame())
    layers   = collect(series)
    sp_names = string.(parent(lookup(series, 1)))
    _assemblage(layers, sp_names, isnothing(mask) ? _default_mask(layers) : mask;
                sitestats = sitestats, traits = traits)
end

function SE.Assemblage(stack::RasterStack, mask = nothing;
                       sitestats::DataFrame = DataFrame(), traits::DataFrame = DataFrame())
    ks       = keys(stack)
    layers   = [stack[k] for k in ks]
    sp_names = string.(collect(ks))
    _assemblage(layers, sp_names, isnothing(mask) ? _default_mask(layers) : mask;
                sitestats = sitestats, traits = traits)
end

# ----------------------------------------------------------------------------
# back-conversion:  Assemblage  ->  rasters
# ----------------------------------------------------------------------------

# A per-site vector scattered back onto the (full) domain. Works for full
# assemblages and site-subset views alike — non-selected cells get missingval.
SE.to_raster(values::AbstractVector, asm::SEAssemblage) = SE.to_raster(values, getcoords(places(asm)))

function SE.to_raster(values::AbstractVector, rd::AnyRasterData)
    length(values) == length(rd.cellinds) ||
        throw(DimensionMismatch("got $(length(values)) values for $(length(rd.cellinds)) sites"))
    T   = eltype(values)
    mv  = T <: AbstractFloat ? T(NaN) : zero(T)
    out = Raster(fill(mv, dims(rd.mask)); missingval = mv)
    out[rd.cellinds] .= values
    out
end

function SE.to_rasterseries(asm::SEAssemblage)
    rd = getcoords(places(asm))
    rd isa AnyRasterData || error("to_rasterseries requires a raster-backed assemblage")
    occ = occurrences(asm)
    sp  = speciesnames(asm)
    rasters = map(1:nspecies(asm)) do i
        r = Raster(falses(size(rd.mask)), dims(rd.mask); missingval = false, name = Symbol(sp[i]))
        r[rd.cellinds] .= Vector(occ[i, :])
        r
    end
    RasterSeries(rasters, (; name = collect(sp)))
end

SE.richness_raster(asm::SEAssemblage) = SE.to_raster(vec(sum(occurrences(asm), dims = 1)), asm)

# ----------------------------------------------------------------------------
# coarsen (the raster-backed grid) — uses Rasters.aggregate under the hood
# ----------------------------------------------------------------------------

"""
    coarsen(asm::Assemblage{Bool, …RasterData…}, factor; fun = maximum)

Coarsen a raster-backed assemblage by `factor` (an `Integer`, or a
`(x, y)` tuple). Each species is occupied in a coarse cell when `fun` over the
covered fine cells is true — the default `maximum` means *any* fine cell, the
presence-absence analogue of the grid `coarsen`. The coarse domain keeps any
fine cell that was in the original domain. Species traits are carried over.
"""
function SE.coarsen(asm::Assemblage{Bool, <:Locations{<:RasterData}},
                    factor::Union{Integer, Tuple{Integer, Integer}}; fun = maximum)
    rd        = getcoords(places(asm))
    series    = SE.to_rasterseries(asm)
    agglayers = [Rasters.aggregate(fun, r, factor) for r in series]
    aggmask   = Rasters.aggregate(maximum, rd.mask, factor)
    aggseries = RasterSeries(agglayers, (; name = collect(speciesnames(asm))))
    SE.Assemblage(aggseries, aggmask; traits = traits(asm))
end

# Make Rasters' own `aggregate` verb work on SpatialEcology assemblages, routing
# to `coarsen`. Defining methods on Rasters.aggregate for *our* types is not
# piracy. The upshot: with Rasters loaded, `aggregate(asm, 2)` works unqualified
# (there is no name clash, since SpatialEcology no longer exports `aggregate`),
# in addition to the always-available `coarsen(asm, 2)`.
Rasters.aggregate(fun, asm::SEAssemblage, factor) = SE.coarsen(asm, factor; fun = fun)
Rasters.aggregate(asm::SEAssemblage, factor; fun = maximum) = SE.coarsen(asm, factor; fun = fun)

end # module
