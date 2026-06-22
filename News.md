# version 0.10.0
- **breaking:** the grid-coarsening function `aggregate` is renamed to `coarsen`, and `aggregate` is no longer exported. This avoids a name clash with `Rasters.aggregate` (which SpatialEcology does not own). Replace `aggregate(assemblage, …)` with `coarsen(assemblage, …)`. With Rasters loaded, `aggregate(assemblage, …)` keeps working, because the extension teaches `Rasters.aggregate` to accept assemblages.
- add a [Rasters.jl](https://github.com/rafaqz/Rasters.jl) integration: build an `Assemblage` from a `RasterStack` or `RasterSeries` of per-species `Bool` rasters (plus an optional domain mask), convert results back with `to_raster`, `to_rasterseries` and `richness_raster`, and `coarsen` raster-backed assemblages
- Rasters is a weak dependency loaded through a package extension, so it is only pulled in when you `using Rasters` yourself — it adds no cost to non-raster workflows
- minimum Julia version is now 1.9 (package extensions)
- **breaking:** removed the vestigial, exported `SiteData` type (and the internal `SubSiteData`/`SESpatialData`). Its `show` and accessors were broken and nothing in the normal workflow produced one.
- fixed `getindex(::Locations, inds)`, which previously errored; it now returns a correct site-subset of the locations.
- fixed dropping of empty sites: `Assemblage(...; dropemptysites = true)` now works for both grid and point data, and the `dropemptyspecies`/`dropemptysites` keywords are now honoured when constructing from a matrix or `ComMatrix` + coordinates.

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
