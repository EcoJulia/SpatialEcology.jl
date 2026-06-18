# Raster data (Rasters.jl integration)

SpatialEcology can build an `Assemblage` directly from
[Rasters.jl](https://rafaqz.github.io/Rasters.jl/stable/) data and convert
results back to rasters for mapping. This integration ships as a *package
extension*: it activates automatically once Rasters is loaded alongside
SpatialEcology, and adds no load-time cost for users who don't need it.

```@example raster
using SpatialEcology, Rasters
```

## Building an assemblage from rasters

A presence–absence dataset is one `Bool` raster per species over a shared grid,
together with a `Bool` *mask* marking which cells are sites (the study domain).
The species layers can be a `RasterStack` (one layer per species) or a
`RasterSeries`; layer names become species names.

```@example raster
x, y = X(1.0:6.0), Y(1.0:5.0)          # X is longitude, Y is latitude

species = RasterStack((
    sparrow = Raster([i <= 3        for i in 1:6, j in 1:5], (x, y)),
    finch   = Raster([iseven(i + j) for i in 1:6, j in 1:5], (x, y)),
    owl     = Raster([j >= 3        for i in 1:6, j in 1:5], (x, y)),
))

asm = Assemblage(species)
```

With no `mask` given, the domain defaults to the cells that hold valid
(non-missing) data in every layer — here, the whole grid. Pass an explicit mask
to restrict the assemblage to a study region:

```@example raster
mask = Raster([i + j <= 9 for i in 1:6, j in 1:5], (x, y))
asm = Assemblage(species, mask)
```

Inputs in either `(X, Y)` or `(Y, X)` dimension order are accepted — rasters
read from disk (e.g. GeoTIFFs) are often `(Y, X)` — and are normalized
internally, so coordinates always come out correctly. Each species layer is
checked to share the mask's dimensions and extent.

Per-site covariates and species traits go in `sitestats` and `traits`:

```@example raster
using DataFrames
asm = Assemblage(species, mask;
                 traits = DataFrame(name = ["sparrow", "finch", "owl"], mass = [25, 18, 900]))
```

## Working with the assemblage

The result is an ordinary SpatialEcology `Assemblage` — every operation
(richness, occupancy, views, grouping, randomization, …) works unchanged:

```@example raster
richness(asm)
```

```@example raster
v = view(asm, species = ["sparrow", "owl"])
nspecies(v), nsites(v)
```

## Converting back to rasters

`to_raster` scatters a per-site vector back onto the full grid (cells outside the
domain become `missing`); `richness_raster` is a shortcut for site richness, and
`to_rasterseries` round-trips the occurrences to one raster per species:

```@example raster
richness_raster(asm)
```

```@example raster
to_raster(Float64.(richness(asm)), asm)
```

These work on views too — a site subset is drawn on the same map, with
unselected cells left missing.

## Coarsening the grid

`coarsen` lumps a raster-backed assemblage to a coarser grain by an integer
factor (or an `(x, y)` tuple), using `Rasters.aggregate` underneath. The default
`fun = maximum` gives presence-absence "*any* fine cell occupied" semantics.

```@example raster
coarse = coarsen(asm, 2)
nsites(asm), nsites(coarse)
```

`coarsen` is the same verb used for ordinary gridded assemblages, so it always
works whether or not Rasters is loaded. As a convenience, with Rasters loaded
the extension also teaches `Rasters.aggregate` to accept assemblages, so the
familiar Rasters verb works too — `aggregate(asm, 2)` is equivalent to the call
above. (SpatialEcology does not export `aggregate` itself, so there is no name
clash with Rasters.)
