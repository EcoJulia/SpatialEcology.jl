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
# Lookup traits (sampling/order/intervalbounds) are not exported by Rasters, and
# `order` would clash with DataFrames, so reach them through the module.
const LU = Rasters.DimensionalData.Lookups

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

# Cell boundaries read off the raster's own lookup: n + 1 of them, in the
# lookup's own direction, so a ReverseOrdered Y - how a north-up raster is
# normally stored - yields descending edges and everything EcoBase derives
# from them descends with it.
#
# Without these, EcoBase falls back to `xedges(::AbstractGrid)`, which rebuilds
# the edges as `xmin:xcellsize:...`. That is always ascending and needs a
# constant step, so it disagreed with this grid's own `xrange`/`yrange` on both
# counts. Pairing that ascending range with row indices read from a descending
# lookup is what drew north-up rasters upside down.
function _lookupedges(l)
    v = parent(l)
    n = length(v)
    if LU.sampling(l) isa LU.Intervals
        # Real cell bounds. Each pair is (low, high) in value order whatever
        # the lookup's direction, so take whichever end the lookup meets first.
        b = LU.intervalbounds(l)
        lead = LU.order(l) isa LU.ReverseOrdered ? last : first
        trail = LU.order(l) isa LU.ReverseOrdered ? first : last
        return vcat([lead(x) for x in b], trail(b[end]))
    end
    # Points sampling, which is what a plain `X(1.0:4.0)` gives and what
    # `xrange`/`yrange` above report: each value sits inside its cell, so cut
    # halfway to each neighbour and extend the two outer half-cells by the
    # spacing beside them.
    n == 1 && return [v[1], v[1]]
    mids = [(v[i] + v[i + 1]) / 2 for i in 1:(n - 1)]
    return vcat(v[1] - (v[2] - v[1]) / 2, mids,
                v[end] + (v[end] - v[end - 1]) / 2)
end

EcoBase.xedges(rd::AnyRasterData) = _lookupedges(lookup(rd.mask, X()))
EcoBase.yedges(rd::AnyRasterData) = _lookupedges(lookup(rd.mask, Y()))

# Where in its cell a reported coordinate sits, read off the lookup's sampling
# rather than assumed. Rasters read off disk are commonly Intervals(Start()),
# whose coordinates are the cell's lower corner; EcoBase defaults to
# CellCentre(), so leaving this undeclared put every such coordinate half a
# cell out of place.
#
# DimensionalData's three loci are exact and direction-independent: the index
# value sits at the minimum, the midpoint and the maximum of its own interval
# for Start, Center and End respectively, whichever way the lookup runs.
# Points sampling has no locus, and _lookupedges above puts each value at the
# centre of the cell it builds, so it is CellCentre().
function _lookupanchor(l)
    LU.sampling(l) isa LU.Intervals || return EcoBase.CellCentre()
    loc = LU.locus(l)
    loc isa LU.Start && return EcoBase.CellCorner()
    loc isa LU.Center && return EcoBase.CellCentre()
    # An End locus puts the coordinate at the cell's upper edge. EcoBase has
    # CellCentre and CellCorner, the latter being the *lower* corner, so there
    # is nothing honest to return; say so rather than be half a cell wrong.
    return error("this raster's $(nameof(typeof(loc))) locus puts each " *
                 "coordinate at its cell's upper edge, which EcoBase's " *
                 "AbstractCellAnchor cannot express — shift the raster to " *
                 "Start or Center sampling first")
end

function EcoBase.cellanchor(rd::AnyRasterData)
    ax = _lookupanchor(lookup(rd.mask, X()))
    ay = _lookupanchor(lookup(rd.mask, Y()))
    ax === ay || error("this raster's X and Y lookups disagree about where " *
                       "in the cell their coordinates sit ($ax vs $ay), and " *
                       "EcoBase declares one anchor for the whole grid")
    return ax
end

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
SE.to_raster(values::AbstractVector, asm::SEAssemblage; kwargs...) = SE.to_raster(values, getcoords(places(asm)); kwargs...)

# `missingval` is the value off-domain (non-site) cells take in the output
# raster. It defaults to `NaN` for floats and `zero` otherwise (the historical
# behaviour); pass e.g. `missingval = NaN` to scatter an integer count vector
# onto a float raster whose off-domain cells render as missing/transparent.
function SE.to_raster(values::AbstractVector, rd::AnyRasterData;
                      missingval = eltype(values) <: AbstractFloat ? eltype(values)(NaN) : zero(eltype(values)))
    length(values) == length(rd.cellinds) ||
        throw(DimensionMismatch("got $(length(values)) values for $(length(rd.cellinds)) sites"))
    T   = promote_type(eltype(values), typeof(missingval))
    mv  = convert(T, missingval)
    out = Raster(fill(mv, dims(rd.mask)); missingval = mv)
    out[rd.cellinds] .= values
    out
end

function SE.to_rasterseries(asm::SEAssemblage)
    rd = getcoords(places(asm))
    rd isa AnyRasterData || error("to_rasterseries requires a raster-backed assemblage")
    sp = speciesnames(asm)
    rasters = map(1:nspecies(asm)) do i
        r = Raster(falses(size(rd.mask)), dims(rd.mask); missingval = false, name = Symbol(sp[i]))
        r[rd.cellinds[occupied(asm, i)]] .= true
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
