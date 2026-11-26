#------------------
# coarsen

"""
    coarsen(object, grid [, fun]; xmin, ymin)

Coarsen `object` (either an `Assemblage` or `Locations` type) to `grid`. Point
data are gridded, and gridded data are lumped into larger cells.

`grid` can be a `GridTopology`, a number, or an `(x, y)` tuple of numbers:

- for point data a number is the cellsize in coordinate units, and may be any
  positive number. Cell edges fall on multiples of the cellsize, and a point on
  an edge belongs to the cell above or to the right. The exception is the North
  Pole: if the coordinates can be longitudes and latitudes (x within -180 to 360
  and y within -90 to 90), a point at latitude 90 belongs to the cell below it.
- for gridded data a number is how many cells to merge along each axis, and must
  be a positive `Integer`. The blocks are anchored to the coordinate origin, not
  to the corner of the grid (if the old cell edges fall on multiples of the old
  cellsize, the new ones fall on multiples of the new cellsize), so the outer
  cells may be only partly covered by the old grid.

With a `GridTopology`, each point, or the centre of each cell, goes to the cell
of the new grid that holds it, and it is an error if some site falls outside.

`fun` is an optional function specifying how to lump occurrences. If not
specified the default is `any` for Boolean Assemblages and `sum` for Integer
ones. With a numeric `grid`, `xmin` and `ymin` optionally set the centre of the
first cell.

Raster-backed assemblages have their own method for an integer factor, which
counts the blocks from the first cell of the raster.
"""
coarsen(lo::SELocations, factor; kw...) = Locations{GridData}(coarsen(lo.coords, factor; kw...))
coarsen(gr::GridTopology, factor) = _grid_from_factor(gr, factor)
function coarsen(coords::Union{SEGrid, SEPoints}, factor::Union{Real, Tuple{Real, Real}}; xmin = nothing, ymin = nothing)
    coarsen(coords, _grid_from_factor(coords, factor; xmin = xmin, ymin = ymin))
end
function coarsen(coords::Union{SEGrid, SEPoints}, newgrid::GridTopology)
    ind = _apply_grid(coords, newgrid)
    GridData(sortslices(unique(ind, dims = 1), dims = 1), newgrid)
end

function coarsen(asm::SEAssemblage{D, T, P}, factor::Union{Real, Tuple{Real, Real}}, fun = _default_fun(D); xmin = nothing, ymin = nothing) where D where T where P <: SELocations
    gt = _grid_from_factor(asm.site.coords, factor; xmin = xmin, ymin = ymin)
    coarsen(asm, gt, fun)
end

function coarsen(asm::SEAssemblage{D}, gt::GridTopology, fun = _default_fun(D)) where D
    # create the new grid
    tmpgrid = _apply_grid(asm.site.coords, gt)
    tmpsites = sortslices(unique(tmpgrid, dims = 1), dims = 1)

    # fill a new occurrence matrix
    retmat = spzeros(D, nspecies(asm), size(tmpsites, 1))

    # get a dictionary of the indices of each new group
    groupinds = Dict{typeof(first(eachrow(tmpgrid))), Vector{Int}}()
    for (i, row) in enumerate(eachrow(tmpgrid))
        if haskey(groupinds, row)
            push!(groupinds[row], i)
        else
            groupinds[row] = [i]
        end
    end

    for newcell in axes(tmpsites, 1)
        inds = groupinds[tmpsites[newcell, :]]
        if fun === Base.any
            retmat[nzrows(view(occurrences(asm), :, inds)), newcell] .= true
        elseif fun === Base.sum
            retmat[:, newcell] .= vec(rowsum(view(occurrences(asm), :, inds)))
        else
            retmat[:, newcell] .= vec(mapslices(fun, occurrences(asm)[:, inds], dims = 2))
        end
    end

    #build the basic objects
    newsitenames = string.(["agg_"], 1:size(retmat, 2))
    sit = Locations{GridData}(GridData(tmpsites, gt), DataFrame(sites = newsitenames))
    retcommat = ComMatrix{D}(retmat, speciesnames(asm), newsitenames)
    sp = SpeciesData(retcommat, copy(traits(asm)))
    Assemblage(sit,sp)
end

#-----------
# `coarsen` helper functions

const _TOL = 2 * sqrt(eps()) # in units of cells, the tolerance `getindices` uses

# the new grid must have a cell for every site
function _check_covered(ind, grid)
    all(i -> 1 <= i <= xcells(grid), view(ind, :, 1)) && all(i -> 1 <= i <= ycells(grid), view(ind, :, 2)) ||
        throw(ArgumentError("the new grid does not cover all sites"))
    ind
end

# each cell goes to the new cell that holds its centre
function _apply_grid(gr::SEGrid, newgrid::GridTopology)
    centres = coordinates(gr, EcoBase.XThenY(), EcoBase.CellCentre())
    _check_covered(getindices(Matrix{Float64}(centres), newgrid), newgrid)
end
function _apply_grid(pt::SEPoints, newgrid::GridTopology)
    coords = coordinates(pt)
    ind = getindices(coords, newgrid)
    if _islonlat(coords)
        # latitude 90 is on the top edge of the world, so it belongs to the cell below
        ind[:, 2] .= min.(ind[:, 2], _cells_below(90, ymin(newgrid), ycellsize(newgrid)))
    end
    _check_covered(ind, newgrid)
end

# whether the coordinates can be longitudes and latitudes
_islonlat(coords) = all(x -> -180 <= x <= 360, view(coords, :, 1)) && all(y -> -90 <= y <= 90, view(coords, :, 2))

# the number of cells, counted from the one centred at `start`, that begin below `top`
_cells_below(top, start, cellsize) = ceil(Int, (top - start) / cellsize + 0.5 - _TOL)

# The cell centres covering the coordinates between `mi` and `ma`. Cell edges fall
# on multiples of the cellsize (unless `inmin` places the first centre), cells are
# counted with the same rule as `getindices` so that every point has a cell, and
# no cell begins at or above `top`.
function _range_for_points(mi, ma, cellsize; inmin = nothing, top = Inf)
    cellsize = Float64(cellsize)
    newmin = isnothing(inmin) ?
        (min(floor(mi / cellsize + _TOL), ceil(top / cellsize - _TOL) - 1) + 0.5) * cellsize : inmin
    ncells = 1 + floor(Int, (ma - newmin) / cellsize + 0.5 + _TOL)
    isfinite(top) && (ncells = min(ncells, _cells_below(top, newmin, cellsize)))
    # never an empty range: a first centre beyond all the data is caught by `_check_covered`
    range(newmin, step = cellsize, length = max(ncells, 1))
end

# The cell centres of the grid merging blocks of `n` cells of `r`. The blocks are
# anchored to the coordinate origin rather than to the first cell, so the outer
# ones may be only partly covered, and grids sharing their cells coarsen onto the
# same coarse grid. `inmin` places the first centre instead.
function _merged_range(r, n; inmin = nothing)
    h = step(r)
    isnothing(inmin) || return _range_for_points(first(r), last(r), n * h; inmin = inmin)
    q = floor(Int, first(r) / h + _TOL)     # number of the first cell, counted from the origin
    offset = first(r) - q * h               # where the centres sit within their cells
    b1, b2 = fld(q, n), fld(q + length(r) - 1, n)
    range(offset + (b1 * n + (n - 1) / 2) * h, step = n * h, length = b2 - b1 + 1)
end

function _grid_from_factor(pt::SEPoints, cellsize; xmin = nothing, ymin = nothing)
    all(>(0), cellsize) || throw(ArgumentError("the cellsize must be positive, got $cellsize"))
    coords = coordinates(pt)
    (xmi, xma), (ymi, yma) = mapslices(extrema, coords, dims = 1)
    top = _islonlat(coords) ? 90 : Inf
    GridTopology(_range_for_points(xmi, xma, first(cellsize); inmin = xmin),
        _range_for_points(ymi, yma, last(cellsize); inmin = ymin, top = top))
end

_check_factor(factor) = all(f -> f isa Integer && f > 0, factor) || throw(ArgumentError(
    "coarsening gridded data needs a positive integer factor (or a tuple of two), got $factor; a non-integer cellsize is only supported for point data"))

function _grid_from_factor(gd::SEGrid, factor; xmin = nothing, ymin = nothing)
    _check_factor(factor) # before reaching for `gd.grid`, which raster data do not have
    _grid_from_factor(gd.grid, factor; xmin = xmin, ymin = ymin)
end
function _grid_from_factor(gd::GridTopology, factor; xmin = nothing, ymin = nothing)
    _check_factor(factor)
    GridTopology(_merged_range(gd.xs, first(factor); inmin = xmin), _merged_range(gd.ys, last(factor); inmin = ymin))
end

_default_fun(::Type{Bool}) = Base.any
_default_fun(::Type{<:Integer}) = Base.sum
_default_fun(::Any) = error("Default aggregation functions are only defined for Assemblage{Bool} (presence-absence) and Assemblage{Int}")
