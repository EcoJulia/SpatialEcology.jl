# version 0.11.0
- **breaking:** requires [EcoBase 0.2](https://github.com/EcoJulia/EcoBase.jl/pull/25). EcoBase 0.2 inserts `AbstractGridded` and `AbstractAreas` between `AbstractGrid` and `AbstractLocationData`, adds cell-anchor and coordinate-order arguments to `indices` and `coordinates`, and answers gridded questions for the places and assemblages that hold a grid — so this package no longer forwards `xmin`, `indices`, `xrange` and friends itself. `SEGrid` still reaches `AbstractGrid`, so no SpatialEcology type changes. Note that other EcoBase consumers may still be pinned to EcoBase 0.1.
- **raster-backed assemblages can now be plotted.** `indices` reached for a field a `RasterData` does not have (it stores `cellinds`), so plotting one threw a `FieldError`; it now has its own method, built from `cellinds` and so in the same site order as `coordinates`.
- **fixed:** north-up rasters were drawn upside down. `xedges`/`yedges` for a `RasterData` are now read off the raster's own lookup, in the lookup's own direction, rather than being rebuilt from a start and a step. The rebuilt edges always ascended, so pairing them with row indices from a `ReverseOrdered` Y lookup — how anything read off disk with Rasters normally comes — flipped the map vertically. As a side effect the grid interface now also works over an `Irregular` lookup, where the old constant-step reconstruction threw.
- **fixed:** a raster whose lookup is `Intervals(Start())` — how rasters read off disk commonly come — reports coordinates at each cell's lower corner, but SpatialEcology left EcoBase to assume they were cell centres, putting every coordinate half a cell out. `cellanchor` is now read off the lookup's locus, so `coordinates(rd, order, CellCentre())` shifts them and `CellCorner()` does not. A `Center()` locus and plain (`Points`) sampling are `CellCentre()` as before, so nothing changes for the common case; an `End()` locus, which EcoBase's `AbstractCellAnchor` cannot express, now errors rather than being silently half a cell wrong.
- `indices(g::SEGrid, idx)` now types its selector as `idx::Integer`. Left untyped it was mutually ambiguous with EcoBase's `indices(::AbstractGridded, ::AbstractCoordinateOrder)`, so asking a grid for a coordinate order threw instead of answering. `indices(grid, 1)` is unaffected.

# version 0.10.1
- add public accessors `domain` and `cellindices` for raster-backed assemblages (and their `Locations` / `RasterData`): `domain` returns the `Bool` domain raster, `cellindices` the per-site `CartesianIndex`es into it. These replace reaching into `getcoords(places(asm)).mask` / `.cellinds`.
- `to_raster` gains a `missingval` keyword controlling the value off-domain cells take in the output raster (default unchanged: `NaN` for floats, `zero` otherwise). Pass `missingval = NaN` to scatter an integer vector onto a float raster whose off-domain cells render as missing.
- **performance**: `ComMatrix` now stores `occurrences_t` — the CSC transpose of the occurrence matrix (`n_sites × n_species`). Because the main matrix is `n_species × n_sites` CSC, per-species queries previously required an O(nnz_total) row scan; they are now O(nnz_per_species) column reads. Affected operations (all on the full `ComMatrix`; `SubComMatrix` views keep the prior path):
  - `occupied(com, species)` — the main beneficiary; called per-species in range-shape, centroid, hull-area, and ellipse-fitting loops
  - `getspecies(com, species)` — returns a column view of `occurrences_t` (lazy, O(1) pointer)
  - `noccupied(com, species)` — now O(1) via `length(nzrange(...))`
  - `occupancy(com)` and `speciestotals(com)` for `Bool` assemblages — O(n_species) colptr diff instead of O(nnz_total) rowval scan

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
