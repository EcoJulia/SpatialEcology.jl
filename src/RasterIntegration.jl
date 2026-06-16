# RasterIntegration.jl
# Constructors and converters that bridge Rasters.jl and SpatialEcology.
# Requires Rasters to be loaded by the calling package/script.

# -------------------------------------------------------------------------
# Assemblage constructor from a RasterSeries{Bool} + a Bool mask raster.
#
# The RasterSeries provides species × sites occurrence data; the mask
# defines which cells are valid sites. Per-site environmental data (e.g.
# PCA axes) can be passed in `sitestats`; species traits in `traits`.
# -------------------------------------------------------------------------

function Assemblage(series::Rasters.RasterSeries, mask;
                    sitestats::DataFrames.DataFrame = DataFrames.DataFrame(),
                    traits::DataFrames.DataFrame   = DataFrames.DataFrame())

    sp_names = string.(parent(Rasters.lookup(series, 1)))
    n_sp     = length(series)
    n_site   = count(mask)

    # Build sparse Bool ComMatrix: one column per occupied site, one row per species.
    occ = SparseArrays.sparse(
        reduce(hcat, [SparseArrays.sparse(vec(Bool.(r[mask]))) for r in series])')
    # occ is n_species × n_sites

    site_ids = string.(1:n_site)
    commat   = ComMatrix{Bool}(occ, sp_names, site_ids)
    spdata   = SpeciesData(commat,
                   isempty(traits) ? DataFrames.DataFrame(name = sp_names) : traits)

    # sitestats must carry a :sites column (used by sitenames/placenames). Copy
    # the user's frame so we don't mutate it, and add :sites if absent.
    ss = isempty(sitestats) ? DataFrames.DataFrame() : copy(sitestats)
    "sites" in names(ss) || (ss[!, :sites] = site_ids)

    locdata = Locations{RasterData}(RasterData(mask), ss)
    Assemblage(locdata, spdata)
end

# -------------------------------------------------------------------------
# Back-conversion: per-site vector → Raster (for plotting / spatial ops)
# -------------------------------------------------------------------------
# Works for full assemblages and site-subset views alike: values are scattered
# at `cellinds` onto the (full) domain mask, so a subset plots on the same map
# with non-selected cells set to missingval.

to_raster(values::AbstractVector, asm::SEAssemblage) = to_raster(values, getcoords(places(asm)))

function to_raster(values::AbstractVector, rd::AnyRasterData)
    length(values) == length(rd.cellinds) ||
        throw(DimensionMismatch("got $(length(values)) values for $(length(rd.cellinds)) sites"))
    T   = eltype(values)
    mv  = T <: AbstractFloat ? T(NaN) : zero(T)
    out = Rasters.Raster(fill(mv, Rasters.dims(rd.mask)); missingval = mv)
    out[rd.cellinds] .= values
    out
end

# -------------------------------------------------------------------------
# Back-conversion: Assemblage → RasterSeries (full round-trip)
# -------------------------------------------------------------------------

function to_rasterseries(asm::SEAssemblage)
    rd = getcoords(places(asm))
    rd isa AnyRasterData || error("to_rasterseries requires a raster-backed assemblage")
    occ = occurrences(asm)
    sp  = speciesnames(asm)
    rasters = map(1:nspecies(asm)) do i
        r = Rasters.Raster(falses(size(rd.mask)), Rasters.dims(rd.mask);
                           missingval = false, name = Symbol(sp[i]))
        r[rd.cellinds] .= Vector(occ[i, :])
        r
    end
    Rasters.RasterSeries(rasters, (; name = collect(sp)))
end

# -------------------------------------------------------------------------
# Convenience: richness as a Raster
# -------------------------------------------------------------------------

richness_raster(asm::SEAssemblage) = to_raster(vec(sum(occurrences(asm), dims = 1)), asm)
