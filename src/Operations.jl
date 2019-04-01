
#aggregate(sd::SpeciesData{D}) where D =
#function aggregate(com::ComMatrix{T}, factor, fun = :auto) where {T}

aggregate(lo::SELocations, factor) = Locations{GridData}(aggregate(lo.coords, factor))
aggregate(gr::GridTopology, factor) = GridTopology(gr.xs[1]:step(gr.xs)*factor:gr.xs[end], gr.ys[1]:step(gr.ys)*factor:gr.ys[end])


#aggregate(asm::SEAssemblage{D, T, P}, factor::Integer, fun = :auto; data_function = :auto) where D where T where P <: SEGrid =
#    Assemblage{D, GridData}(aggregate(asm.site, factor, fun; data_function), aggregate_grid(asm.occ, factor, fun; data_function))

aggregate(gr::SEGrid, factor::Integer) = GridData(ceil.(Int, gr.indices ./ factor), aggregate(gr.grid, factor))
function aggregate(gr::SEGrid, newgrid::GridTopology)
    factor = [cellsize(newgrid)...] ./ cellsize(gr)
    shift = [xmin(newgrid) - xmin(gr), ymin(newgrid) - ymin(gr)] .% factor ./cellsize(gr)
    newind = ceil.(Int, (shift' .+ gr.indices) ./ factor')
    newind, newgrid
end


replace_auto(x, ::Type{Bool}) = x == :auto ? Base.any : x
replace_auto(x, ::Type{Integer}) = x == :auto ? Base.sum : x

function aggregate(asm::SEAssemblage{D, T, P}, factor::Integer, fun = :auto) where D where T where P <: SELocations{<: SEGrid}
    _fun = replace_auto(fun, D)
    _fun == :auto && error("Default aggregation functions are only defined for Bool (presence-absence) and Int (abundance) commatrices")
    newx, newy = ceil.(Int,cells(asm)./factor)
    newsites = newx*newy
    retmat = spzeros(D, nspecies(asm), newsites)
    #_fun(x) = _fun(filter(y -> isfinite(y) && !ismissing(y), x)) not sure I have missing values

    for j in 1:newx-1
       for i in 1:newy-1
          newi, newj = (i-1)*factor .+ (1:factor), (j-1)*factor .+ (1:factor)
          window = findall(x->((ni,nj) = indices(asm, x); ni ∈ newi && nj ∈ newj), 1:nsites(asm))
          length(window) > 0 &&
            (retmat[:, (j-1)*newy + i] = mapslices(_fun, view(occurrences(asm), :, window), dims = 2))
       end
       newi, newj = (newx-1)*factor .+ (1:xcells(asm)), (j-1)*factor .+ (1:factor)
       window = findall(x->((ni,nj) = indices(asm, x); ni ∈ newi && nj ∈ newj), 1:nsites(asm))
       length(window) > 0 &&
         (retmat[:, j*newy] = mapslices(_fun, view(occurrences(asm), :, window), dims = 2))
    end

    for i in 1:newy-1
      newi, newj = (newx-1)*factor + 1:xcells(asm), (newy-1)*factor .+ (1:ycells(asm))
      window = findall(x->((ni,nj) = indices(asm, x); ni ∈ newi && nj ∈ newj), 1:nsites(asm))
      length(window) > 0 &&
        (retmat[:, (newx-1)*newy + i] = mapslices(_fun, view(occurrences(asm), :, window), dims = 2))
    end

    newsitenames = string.(["agg_"], 1:newsites)
    @show size(retmat)
    retcommat = ComMatrix{D}(retmat, speciesnames(asm), newsitenames)
    sp = SpeciesData(retcommat, traits(asm))
    sit = aggregate(asm.site, factor)
    Assemblage(sit,sp)
end
