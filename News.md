# version 0.10.0
- add a [Rasters.jl](https://github.com/rafaqz/Rasters.jl) integration: build an `Assemblage` from a `RasterStack` or `RasterSeries` of per-species `Bool` rasters (plus an optional domain mask), convert results back with `to_raster`, `to_rasterseries` and `richness_raster`, and `aggregate` raster-backed assemblages
- Rasters is a weak dependency loaded through a package extension, so it is only pulled in when you `using Rasters` yourself — it adds no cost to non-raster workflows
- minimum Julia version is now 1.9 (package extensions)

# version 0.4.1
- removed `full` again
- greatly expanded test suite

# version 0.4.0
- julia-0.7 compatibility

# version 0.3.1
- several smaller patch commits

# version 0.3.0
- depend on DataFrames 0.11 and Missings

# version 0.2.0
- change from species-as-columns to sites-as-columns, which transposes all functions on ComMatrix
- add sitecolumns keyword argument to transpose input matrices not following this
- Remove RCall dependency and functionality
- various improvements to speed
- add functions
- - occurrences
- - cooccurences
- - asquantiles
- - full
- - sitetotals
- - speciestotals
