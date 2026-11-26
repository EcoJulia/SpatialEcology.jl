# version 0.12.0
- **fixed:** indexing a `ComMatrix` ([#22](https://github.com/EcoJulia/SpatialEcology.jl/issues/22)). `com[i, j]` returns the value; any other selection returns a `ComMatrix` that keeps the names of its species and sites (a single species or site stays a dimension, so `com[1:3, 1]` is 3×1). Positions, ranges, `:`, `end`, Bool masks and species/site names all work as indices, and an unknown name is an error. `@view com[i, j]` gives a `SubComMatrix`. Before, single values and vectors threw, and other selections renamed everything `species1…`/`site1…`.
- **fixed:** `Assemblage(...; dropemptyspecies = true)` renamed the species to `species1…`, so they no longer matched the traits. `dropemptysites = true` subset the occurrences but not the community matrix's site names, so name lookups such as `occupied(cm, "sp")` answered for the wrong sites; it also threw for raster-backed sites.
- **fixed:** the transposed occurrence matrix kept since 0.10.1 for fast species queries went stale after `setindex!`, a write through a view, or `dropemptysites`, so `occupied`, `noccupied`, `occupancy` and `getspecies` answered from the old data. It now lives with the occurrences in a new storage type, `OccMatrix`, which drops it on every write. Writing a zero also left a stored zero that `occupied`, `noccupied`, `occupancy` and `nrecords` counted as a presence; writes now remove it.
- **breaking:** `occurrences(com)` returns that `OccMatrix` (an `AbstractSparseMatrixCSC`) instead of a `SparseMatrixCSC`, and `getspecies(com, sp)` returns a `SparseVector` copy instead of a view: the view was into the stored transpose, so it went stale after writes, and writing through it changed the transpose but not the occurrences.
- `coarsen` accepts a non-integer cellsize for point data ([#21](https://github.com/EcoJulia/SpatialEcology.jl/issues/21)), e.g. `coarsen(points, 0.5)` or `coarsen(points, (0.5, 0.25))`. A non-integer factor for gridded data now throws an `ArgumentError` rather than a `MethodError`.
- **fixed:** gridding point data no longer errors when the largest x or y coordinate falls on a cell edge. In longitude/latitude data a point at latitude 90 is placed in the cell below the pole.
- **breaking:** for gridded data, `coarsen(x, n)` now merges blocks of `n × n` cells, as documented, also for a `GridTopology`. Since 2019 `n` had been used as the new cellsize in coordinate units, which is only the same thing for a cellsize of 1; other grids got the wrong grain or an error.
- **breaking:** coarsened cells sit at the centre of the cells they merge. For grids with a cellsize of 1 and cell centres at whole numbers, this moves the coordinates by half a cell.
- **breaking:** coarsening a raster-backed assemblage now marks a species as present in a coarse cell if it is present in *any* of its cells, as documented. Before, a coarse cell was a site only if all its cells were in the domain, and held a species only if the species filled it, whatever `fun` was. `minimum` gives that answer for species.
- **fixed:** coarsening a raster-backed assemblage dropped the incomplete block at the far edge when the size was not a multiple of the factor.
- **fixed:** coarsening gridded data onto a `GridTopology` with a different origin could put cells in the wrong place. A new grid (or `xmin`/`ymin`) that leaves out sites now throws an `ArgumentError`.
- **fixed:** `coarsen` works on views of sites or species, `coarsen(places(x), factor)` works for point data and tuple factors, and `Rasters.aggregate` works again for gridded and point assemblages (it threw a `MethodError`). A coarsened assemblage no longer shares its traits table with the original.
- `addtraits!` and `addsitestats!` with a `DataFrame` are safe to re-run: columns the assemblage already has are replaced, with an `@info` naming them (`makeunique = true` keeps both), and a failing join leaves the assemblage unchanged.

# version 0.11.3
- **fixed:** copying a site view of a gridded `Assemblage`, or dropping its empty sites, threw an `InexactError` for cell sizes that are not exact in binary, such as 1/12 degree.
- **fixed:** gridding coordinates used half the x cellsize as the y offset, giving wrong y indices when the two cellsizes differ.
- requires EcoBase 0.2.2 and RandomBooleanMatrices 0.2.1.

# version 0.11.2
- supports RandomBooleanMatrices 0.2. `matrixrandomizer` on a Boolean `ComMatrix` or `Assemblage` now wraps RandomBooleanMatrices' own generator, and forwards `method` and `trades`.
- **fixed:** random draws ignored the `rng` given to `matrixrandomizer`, and `rand!` left the transposed occurrence matrix stale.

# version 0.11.1
- **fixed:** species and site names are now always stored as `String`, as EcoBase's `thingnames`/`placenames` promise. Names read with CSV.jl 1.x are `DataStrings.DataString`s, which the constructors' `string.()` passed through unchanged, so name selectors such as `view(asm, species = ["sp_a"])` threw a `MethodError` on CSV-read assemblages.
- requires EcoBase 0.2.1, whose `asindices` matches names across string types, so selectors that are not `String`s (a `SubString`, or a CSV column itself) find `String` names too.

# version 0.11.0
- **breaking:** requires [EcoBase 0.2](https://github.com/EcoJulia/EcoBase.jl/pull/25). EcoBase 0.2 inserts `AbstractGridded` and `AbstractAreas` between `AbstractGrid` and `AbstractLocationData`, adds cell-anchor and coordinate-order arguments to `indices` and `coordinates`, and answers gridded questions for the places and assemblages that hold a grid — so this package no longer forwards `xmin`, `indices`, `xrange` and friends itself. `SEGrid` and `GridTopology` now subtype `EcoBase.AbstractRegularGrid`, which is what EcoBase 0.2 renames `AbstractGrid` to — the old name is kept as a deprecated binding for the same type, so this is a spelling change and no SpatialEcology type moved. Note that other EcoBase consumers may still be pinned to EcoBase 0.1.
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
