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
aggregate(gr::GridTopology, factor) = GridTopology(gr.xs[1]:step(gr.xs)*factor:gr.xs[end], gr.ys[1]:step(gr.ys)*factor:gr.ys[end])
aggregate(gr::SEGrid, factor::Integer) = (g = aggregate(gr.grid, factor); aggregate(gr, g))
aggregate(gr::SEGrid, newgrid::GridTopology) = (ind = apply_grid(gr, newgrid); GridData(sortslices(unique(ind, dims = 1), dims = 1), newgrid))

function aggregate(asm::SEAssemblage{D, T, P}, factor::Integer, fun = _default_fun(D)) where D where T where P <: SELocations{<: SEGrid}

    gt = aggregate(asm.site.coords.grid, factor)
    aggregate(asm, gt, fun)
end

function aggregate(asm::SEAssemblage{D}, gt::GridTopology, fun = _default_fun(D)) where D

    # create the new grid
    tmpgrid = apply_grid(asm.site.coords, gt)
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

function apply_grid(gr::SEGrid, newgrid::GridTopology)
    factor = [cellsize(newgrid)...] ./ cellsize(gr)
    shift = [xmin(newgrid) - xmin(gr), ymin(newgrid) - ymin(gr)] .% factor ./cellsize(gr)
    newind = ceil.(Int, (shift' .+ gr.indices) ./ factor')
    newind
end
apply_grid(pt::SEPoints, newgrid::GridTopology) = getindices(coordinates(pt), newgrid)

_default_fun(::Type{Bool}) = Base.any
_default_fun(::Type{Integer}) = Base.sum
_default_fun(::Any) = error("Default aggregation functions are only defined for Assemblage{Bool} (presence-absence) and Assemblage{Int}")
