
#aggregate(sd::SpeciesData{D}) where D =
#function aggregate(com::ComMatrix{T}, factor, fun = :auto) where {T}

aggregate(lo::SELocations, factor) = Locations{GridData}(aggregate(lo.coords, factor))
aggregate(gr::GridTopology, factor) = GridTopology(gr.xs[1]:step(gr.xs)*factor:gr.xs[end], gr.ys[1]:step(gr.ys)*factor:gr.ys[end])


#aggregate(asm::SEAssemblage{D, T, P}, factor::Integer, fun = :auto; data_function = :auto) where D where T where P <: SEGrid =
#    Assemblage{D, GridData}(aggregate(asm.site, factor, fun; data_function), aggregate_grid(asm.occ, factor, fun; data_function))

aggregate(gr::SEGrid, factor::Integer) = (g = aggregate(gr.grid, factor); aggregate(gr, g))
aggregate(gr::SEGrid, newgrid::GridTopology) = (ind = apply_grid(gr, newgrid); GridData(sortslices(unique(ind, dims = 1), dims = 1), newgrid))

function apply_grid(gr::SEGrid, newgrid::GridTopology)
    factor = [cellsize(newgrid)...] ./ cellsize(gr)
    shift = [xmin(newgrid) - xmin(gr), ymin(newgrid) - ymin(gr)] .% factor ./cellsize(gr)
    newind = ceil.(Int, (shift' .+ gr.indices) ./ factor')
    newind
end

replace_auto(x, ::Type{Bool}) = x == :auto ? Base.any : x
replace_auto(x, ::Type{Integer}) = x == :auto ? Base.sum : x

function aggregate(asm::SEAssemblage{D, T, P}, factor::Integer, fun = :auto) where D where T where P <: SELocations{<: SEGrid}
    _fun = replace_auto(fun, D)
    _fun == :auto && error("Default aggregation functions are only defined for Bool (presence-absence) and Int (abundance) commatrices")

    # create the new grid
    gt = aggregate(asm.site.coords.grid, factor)
    tmpgrid = apply_grid(asm.site.coords, gt)
    tmpsites = sortslices(unique(tmpgrid, dims = 1), dims = 1)

    # fill a new occurrence matrix
    retmat = spzeros(D, nspecies(asm), size(tmpsites, 1))
    for newcell in axes(tmpsites, 1)
        inds = findall(x->tmpgrid[x, :] == tmpsites[newcell, :], 1:size(tmpgrid, 1))
        retmat[:, newcell] = mapslices(_fun, view(occurrences(asm), :, inds), dims = 2)
    end

    #build the basic objects
    newsitenames = string.(["agg_"], 1:size(retmat, 2))
    sit = Locations{GridData}(GridData(tmpsites, gt), DataFrame(sites = newsitenames))
    retcommat = ComMatrix{D}(retmat, speciesnames(asm), newsitenames)
    sp = SpeciesData(retcommat, traits(asm))
    Assemblage(sit,sp)
end
