#------------------
# coarsen

"""
    coarsen(object, grid [, fun]; xmin, ymin)

Coarsen `object` (either an `Assemblage` or `Locations` type) to `grid`.
If `object` is an assemblage of point data this grids all points and returns a
gridded assemblage.

`grid` can be a `GridTopology`, a number, or an `(x, y)` tuple of numbers:

- for point data it is the cellsize in coordinate units, and may be any positive
  number. Cell edges fall on multiples of the cellsize, and a point on an edge
  belongs to the cell above or to the right. The exception is the North Pole: if
  the coordinates can be longitudes and latitudes (x within -180 to 360 and y
  within -90 to 90), a point at latitude 90 belongs to the cell below it.
- for gridded data it is the coarsening factor, and must be an `Integer`.

`fun` is an optional function specifying how to lump occurrences. If not
specified the default is `any` for Boolean Assemblages and `sum` for Integer
ones. With a numeric `grid`, `xmin` and `ymin` optionally set the centre of the
first cell.
"""
coarsen(lo::SELocations, factor) = Locations{GridData}(coarsen(lo.coords, factor))
coarsen(gr::GridTopology, factor) = _grid_from_factor(gr, factor)
coarsen(gr::SEGrid, factor::Integer) = (g = _grid_from_factor(gr.grid, factor); coarsen(gr, g))
coarsen(gr::SEGrid, newgrid::GridTopology) = (ind = _apply_grid(gr, newgrid); GridData(sortslices(unique(ind, dims = 1), dims = 1), newgrid))

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
    sp = SpeciesData(retcommat, traits(asm))
    Assemblage(sit,sp)
end

#-----------
# `coarsen` helper functions

_lowerleft(gr) = (xmin(gr), ymin(gr)) .- (0.5 .* cellsize(gr))

function _apply_grid(gr::SEGrid, newgrid::GridTopology)
    factor = [cellsize(newgrid)...] ./ cellsize(gr)
    shift = (_lowerleft(gr) .- _lowerleft(newgrid)) .% factor ./cellsize(gr)
    newind = ceil.(Int, (shift' .+ gr.indices) ./ factor')
    newind
end
function _apply_grid(pt::SEPoints, newgrid::GridTopology)
    coords = coordinates(pt)
    ind = getindices(coords, newgrid)
    if _islonlat(coords)
        # latitude 90 is on the top edge of the world, so it belongs to the cell below
        ind[:, 2] .= min.(ind[:, 2], _cells_below(90, ymin(newgrid), ycellsize(newgrid)))
    end
    ind
end

const _TOL = 2 * sqrt(eps()) # in units of cells, the tolerance `getindices` uses

# whether the coordinates can be longitudes and latitudes
_islonlat(coords) = all(x -> -180 <= x <= 360, view(coords, :, 1)) && all(y -> -90 <= y <= 90, view(coords, :, 2))

# the number of cells, counted from the one centred at `start`, that begin below `top`
_cells_below(top, start, cellsize) = ceil(Int, (top - start) / cellsize + 0.5 - _TOL)

# The cell centres covering the points between `mi` and `ma`. Cell edges fall on
# multiples of the cellsize (unless `inmin` places the first centre), cells are
# counted with the same rule as `getindices` so that every point has a cell, and
# no cell begins at or above `top`.
function _range_for_points(mi, ma, cellsize; inmin = nothing, top = Inf)
    cellsize = Float64(cellsize)
    newmin = isnothing(inmin) ?
        (min(floor(mi / cellsize + _TOL), ceil(top / cellsize - _TOL) - 1) + 0.5) * cellsize : inmin
    ncells = 1 + floor(Int, (ma - newmin) / cellsize + 0.5 + _TOL)
    isfinite(top) && (ncells = min(ncells, _cells_below(top, newmin, cellsize)))
    range(newmin, step = cellsize, length = ncells)
end

function _range_from_factor(mi, ma, factor; inmin = nothing)
    newmin = isnothing(inmin) ? (floor(mi / factor) + 0.5) * factor : inmin
    newmax = isnothing(inmin) ? (ceil(ma/factor) - 0.5) * factor : floor((ma - inmin) / factor) * factor + inmin
    range(newmin, stop = newmax, step = factor)
end

function _grid_from_factor(pt::SEPoints, cellsize; xmin = nothing, ymin = nothing)
    all(>(0), cellsize) || throw(ArgumentError("the cellsize must be positive, got $cellsize"))
    coords = coordinates(pt)
    (xmi, xma), (ymi, yma) = mapslices(extrema, coords, dims = 1)
    top = _islonlat(coords) ? 90 : Inf
    GridTopology(_range_for_points(xmi, xma, first(cellsize); inmin = xmin),
        _range_for_points(ymi, yma, last(cellsize); inmin = ymin, top = top))
end
function _grid_from_factor(gd::SEGrid, factor; xmin = nothing, ymin = nothing)
    all(f -> f isa Integer, factor) || throw(ArgumentError(
        "coarsening gridded data needs an integer factor (or a tuple of two), got $factor; a non-integer cellsize is only supported for point data"))
    _grid_from_factor(gd.grid, factor; xmin = xmin, ymin = ymin)
end
_grid_from_factor(gd::GridTopology, factor; xmin = nothing, ymin = nothing) = GridTopology(_range_from_factor(extrema(gd.xs)..., first(factor); inmin = xmin), _range_from_factor(extrema(gd.ys)..., last(factor); inmin = ymin))

_default_fun(::Type{Bool}) = Base.any
_default_fun(::Type{<:Integer}) = Base.sum
_default_fun(::Any) = error("Default aggregation functions are only defined for Assemblage{Bool} (presence-absence) and Assemblage{Int}")
