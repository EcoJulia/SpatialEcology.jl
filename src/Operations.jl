#------------------
# aggregate

"""
    aggregate(object, grid [, fun])

Aggregate `object` (either an `Assemblage` or `Locations` type) to `grid`.
If `object` is an `Assemblage{PointData}` this will grid all points and return an
`Assemblage{GridData}`. `grid` can be a `GridTopology` or a single `Integer`
signifying the aggregation factor for already gridded data, the cellsize for
point data. `fun` is an optional function specifying how to lump occurrences. If
not specified the default function is `any` for Boolean Assemblages and `sum` for
Integer ones.
"""
aggregate(lo::SELocations, factor) = Locations{GridData}(aggregate(lo.coords, factor))
aggregate(gr::GridTopology, factor) = _grid_from_factor(gr, factor)
aggregate(gr::SEGrid, factor::Integer) = (g = _grid_from_factor(gr.grid, factor); aggregate(gr, g))
aggregate(gr::SEGrid, newgrid::GridTopology) = (ind = _apply_grid(gr, newgrid); GridData(sortslices(unique(ind, dims = 1), dims = 1), newgrid))

function aggregate(asm::SEAssemblage{D, T, P}, factor::Union{Integer, Tuple{Integer, Integer}}, fun = _default_fun(D)) where D where T where P <: SELocations

    gt = _grid_from_factor(asm.site.coords, factor)
    aggregate(asm, gt, fun)
end

function aggregate(asm::SEAssemblage{D}, gt::GridTopology, fun = _default_fun(D)) where D where T where P <: SELocations

    # create the new grid
    tmpgrid = _apply_grid(asm.site.coords, gt)
    tmpsites = sortslices(unique(tmpgrid, dims = 1), dims = 1)

    # fill a new occurrence matrix
    retmat = spzeros(D, nspecies(asm), size(tmpsites, 1))
    for newcell in axes(tmpsites, 1)
        inds = findall(x->tmpgrid[x, :] == tmpsites[newcell, :], 1:size(tmpgrid, 1))
        retmat[:, newcell] = mapslices(fun, view(occurrences(asm), :, inds), dims = 2)
    end

    #build the basic objects
    newsitenames = string.(["agg_"], 1:size(retmat, 2))
    sit = Locations{GridData}(GridData(tmpsites, gt), DataFrame(sites = newsitenames))
    retcommat = ComMatrix{D}(retmat, speciesnames(asm), newsitenames)
    sp = SpeciesData(retcommat, traits(asm))
    Assemblage(sit,sp)
end

#-----------
# `aggregate` helper functions

_lowerleft(gr) = (xmin(gr), ymin(gr)) .- (0.5 .* cellsize(gr))

function _apply_grid(gr::SEGrid, newgrid::GridTopology)
    factor = [cellsize(newgrid)...] ./ cellsize(gr)
    shift = (_lowerleft(gr) .- _lowerleft(newgrid)) .% factor ./cellsize(gr)
    newind = ceil.(Int, (shift' .+ gr.indices) ./ factor')
    newind
end
_apply_grid(pt::SEPoints, newgrid::GridTopology) = getindices(coordinates(pt), newgrid)

_range_from_factor(mi, ma, factor) =
    range((floor(mi / factor) + 0.5) * factor, stop = (ceil(ma/factor) - 0.5) * factor, step = factor)

function _grid_from_factor(pt::SEPoints, factor)
    (xmi, xma), (ymi, yma) = mapslices(extrema, coordinates(pt), dims = 1)
    GridTopology(_range_from_factor(xmi, xma, first(factor)), _range_from_factor(ymi, yma, last(factor)))
end
_grid_from_factor(gd::SEGrid, factor) = _grid_from_factor(gd.grid, factor)
_grid_from_factor(gd::GridTopology, factor) = GridTopology(_range_from_factor(extrema(gd.xs)..., first(factor)), _range_from_factor(extrema(gd.ys)..., last(factor)))

_default_fun(::Type{Bool}) = Base.any
_default_fun(::Type{Integer}) = Base.sum
_default_fun(::Any) = error("Default aggregation functions are only defined for Assemblage{Bool} (presence-absence) and Assemblage{Int}")
