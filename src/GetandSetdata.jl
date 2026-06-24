"""
    coordinates(asm)
"""
coordinates(asm::SEAssemblage) = coordinates(asm.site)
coordinates(pd::SEPoints) = pd.coords
coordinates(l::SELocations) = coordinates(getcoords(l))
function coordinates(gd::SEGrid)
    index1 = xrange(gd.grid)[gd.indices[:,1]]
    index2 = yrange(gd.grid)[gd.indices[:,2]]
    hcat(index1, index2)
end

function coordinates(rd::AnyRasterData)
    # Returns n_sites × 2 [x y], one row per site in `cellinds` (in site order).
    # xrange/yrange dispatch defined in Gridfunctions.jl for the raster types.
    xs = EcoBase.xrange(rd)
    ys = EcoBase.yrange(rd)
    hcat([xs[ci[1]] for ci in rd.cellinds], [ys[ci[2]] for ci in rd.cellinds])
end

getcoords(l::SELocations) = l.coords

"""
    domain(x)

The `Bool` domain raster that a raster-backed `Assemblage` (or its `Locations` /
`RasterData`) is defined on — the raster analogue of a grid's `GridTopology`.
Only defined for raster-backed objects. For a site-subset view the full parent
domain is returned (a subset is still drawn on the same map).
"""
domain(asm::SEAssemblage) = domain(getcoords(places(asm)))
domain(l::SELocations) = domain(getcoords(l))
domain(rd::AnyRasterData) = rd.mask

"""
    cellindices(x)

The vector of `CartesianIndex`es into the [`domain`](@ref) raster, one per site
in site order, for a raster-backed `Assemblage` (or its `Locations` /
`RasterData`). This is what lets a view select an arbitrary, possibly reordered
subset of sites. Only defined for raster-backed objects.
"""
cellindices(asm::SEAssemblage) = cellindices(getcoords(places(asm)))
cellindices(l::SELocations) = cellindices(getcoords(l))
cellindices(rd::AnyRasterData) = rd.cellinds

"""
    traits(x)

Return the dataframe of species traits defined in x (which is usually
an `Assemblage` object)
"""
traits(occ::SEThings) = occ.traits
traits(asm::SEAssemblage) = traits(asm.occ)

places(asm::SEAssemblage) = asm.site
things(asm::SEAssemblage) = asm.occ

"""
    sitestats(asm)
"""
sitestats(asm::SEAssemblage) = asm.site.sitestats

"""
    traitnames(asm)
"""
traitnames(asm::SEAssemblage) = names(traits(asm))

"""
    sitestatnames(asm)
"""
sitestatnames(asm::SEAssemblage) = names(sitestats(asm))

"""
    commatrix(asm)
"""
commatrix(asm::SEAssemblage) = commatrix(asm.occ)
commatrix(occ::SEThings) = occ.commatrix

"""
    occurrences(asm)
"""
occurrences(asm::Union{SEAssemblage, SEThings}) = occurrences(commatrix(asm))
occurrences(cm::AbstractComMatrix) = cm.occurrences

"""
    addtraits!(asm::Assemblage, newtraits::DataFrames.DataFrame, species::Symbol; validate = true, tolerance = 0.5, makeunique = false)
"""
function addtraits!(asm::Assemblage, newtraits::DataFrames.DataFrame, species::Symbol; validate = true, tolerance = 0.5, makeunique = false)
    if validate
        dif, left, right = length(intersect(newtraits[!,species], speciesnames(asm))), nspecies(asm), size(newtraits,1)
        max(dif/left, dif/right) == 0 && error("No match between species names, aborting join")
        println("$dif matching species names,\n",
                "\t$(round(100*dif/left, sigdigits = 3))% of $left species in the Assemblage\n",
                "\t$(round(100*dif/right, sigdigits = 3))% of $right species in the new traits data\n")
        max(dif/left, dif/right) < tolerance && error("Aborting join, as fit was smaller than the tolerance of $tolerance . To perform the join decrease the tolerance value")
    end

    jointraits = rename(newtraits, species => :name)
    asm.occ.traits.__run_number = 1:nrow(asm.occ.traits)
    asm.occ.traits = leftjoin(asm.occ.traits, jointraits, on = :name, makeunique = makeunique)
    sort!(asm.occ.traits, :__run_number)
    select!(asm.occ.traits, Not(:__run_number))
    nothing
end

function addtraits!(asm::Assemblage, newtraits::AbstractVector, name::Union{String, Symbol})
    length(newtraits) == nspecies(asm) || error("Cannot add a vector of length $(length(newtraits)) to an Assemblage with $(nspecies(asm)) species")
    asm.occ.traits[!,Symbol(name)] = newtraits
end

"""
    addsitestats!(asm::Assemblage, newsites::DataFrames.DataFrame, sites::Symbol; validate = true, tolerance = 0.5, makeunique = false)
"""
function addsitestats!(asm::Assemblage, newsites::DataFrames.DataFrame, sites::Symbol; validate = true, tolerance = 0.5, makeunique = false)
    if validate
        dif, left, right = length(intersect(newsites[!,sites], sitenames(asm))) , nsites(asm), size(newsites,1)
        max(dif/left, dif/right) == 0 && error("No match between site names, aborting join")
        println("$dif matching site names,\n",
                "\t$(round(100*dif/left, sigdigits = 3))% of $left sites in the Assemblage\n",
                "\t$(round(100*dif/right, sigdigits = 3))% of $right sites in the new sitestats data\n")
        max(dif/left, dif/right) < tolerance && error("Aborting join, as fit was smaller than the tolerance of $tolerance . To perform the join decrease the tolerance value")
    end

    joinsites = rename(newsites, sites => :sites)
    asm.site.sitestats.__run_number = 1:nrow(asm.site.sitestats)
    asm.site.sitestats = leftjoin(asm.site.sitestats, joinsites, on = :sites, makeunique = makeunique)
    sort!(asm.site.sitestats, :__run_number)
    select!(asm.site.sitestats, Not(:__run_number))
    nothing
end

function addsitestats!(asm::Assemblage, newsites::AbstractVector, name::Union{String, Symbol})
    length(newsites) == nsites(asm) || error("Cannot add a vector of length $(length(newsites)) to an Assemblage with $(nsites(asm)) sites")
    asm.site.sitestats[!,Symbol(name)] = newsites
end
#
# function assemblagejoin!(df1::AbstractDataFrame, df2::AbstractDataFrame, on_left::Symbol, on_right::Symbol)
#     right = DataFrame(Dict(names(df2)[i] => missings(eltypest(df2)[i], nrow(df1)) for i in 1:ncol(df2)))
#     @inbounds for (i, j) in enumerate(indexin(df2[on_right], df1[on_left]))
#         if ! (j == 0)
#             right[j,:] = df2[i,:]
#         end
#     end
#     right = DataFrames.without(right, on_right)
#     DataFrames.hcat!(df1, right)
# end

"""
    @traits(x, expr)
"""
macro traits(x, expr)
    esc(:(@with(traits($x), $expr)))
end

"""
    @sitestats(x, expr)
"""
macro sitestats(x, expr)
    esc(:(@with(sitestats($x), $expr)))
end

"""
    dispersionfield(asm, site)
"""
dispersionfield(asm::EcoBase.AbstractAssemblage, site) =
    occurrences(asm)' * placeoccurrences(asm, site)
