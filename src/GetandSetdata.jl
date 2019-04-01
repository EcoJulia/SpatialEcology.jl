coordinates(asm::SEAssemblage) = coordinates(asm.site)
coordinates(sd::SiteData) = coordinates(sd.site)
coordinates(pd::SEPoints) = pd.coords
coordinates(l::SELocations) = coordinates(getcoords(l))
function coordinates(gd::SEGrid)
    index1 = xrange(gd.grid)[gd.indices[:,1]]
    index2 = yrange(gd.grid)[gd.indices[:,2]]
    hcat(index1, index2)
end

getcoords(l::SELocations) = l.coords

traits(occ::SEThings) = occ.traits
traits(asm::SEAssemblage) = traits(asm.occ)
places(asm::SEAssemblage) = asm.site
things(asm::SEAssemblage) = asm.occ

sitestats(asm::SEAssemblage) = asm.site.sitestats

traitnames(asm::SEAssemblage) = names(traits(asm))
sitestatnames(asm::SEAssemblage) = names(sitestats(asm))

commatrix(asm::SEAssemblage) = commatrix(asm.occ)
commatrix(occ::SEThings) = occ.commatrix

occurrences(asm::Union{SEAssemblage, SEThings}) = occurrences(commatrix(asm))
occurrences(cm::AbstractComMatrix) = cm.occurrences

function addtraits!(asm::Assemblage, newtraits::DataFrames.DataFrame, species::Symbol; validate = true, tolerance = 0.5, makeunique = false)
    if validate
        dif, left, right = length(intersect(newtraits[species], speciesnames(asm))), nspecies(asm), size(newtraits,1)
        max(dif/left, dif/right) == 0 && error("No match between species names, aborting join")
        println("$dif matching species names,\n",
                "\t$(signif(100*dif/left,3))% of $left species in the Assemblage\n",
                "\t$(signif(100*dif/right,3))% of $right species in the new traits data\n")
        max(dif/left, dif/right) < tolerance && error("Aborting join, as fit was smaller than the tolerance of $tolerance . To perform the join decrease the tolerance value")
    end

    nm = names(newtraits)
    rename!(newtraits, species => :name)
    asm.occ.traits = join(asm.occ.traits, newtraits, kind = :left, on = :name, makeunique = makeunique)
    names!(newtraits, nm)
    #assemblagejoin!(asm.occ.traits, newtraits, :name, species)
    nothing
end

function addtraits!(asm::Assemblage, newtraits::AbstractVector, name::Union{String, Symbol})
    length(newtraits) == nspecies(asm) || error("Cannot add a vector of length $(length(newtraits)) to an Assemblage with $(nspecies(asm)) species")
    asm.occ.traits[Symbol(name)] = newtraits
end

function addsitestats!(asm::Assemblage, newsites::DataFrames.DataFrame, sites::Symbol; validate = true, tolerance = 0.5, makeunique = false)
    if validate
        dif, left, right = length(intersect(newsites[sites], sitenames(asm))) , nsites(asm), size(newsites,1)
        max(dif/left, dif/right) == 0 && error("No match between site names, aborting join")
        println("$dif matching site names,\n",
                "\t$(round(100*dif/left, sigdigits = 3))% of $left sites in the Assemblage\n",
                "\t$(round(100*dif/right, sigdigits = 3))% of $right sites in the new sitestats data\n")
        max(dif/left, dif/right) < tolerance && error("Aborting join, as fit was smaller than the tolerance of $tolerance . To perform the join decrease the tolerance value")
    end

    #assemblagejoin!(asm.site.sitestats, newsites, :sites, sites) #TODO this should instead be on the sitenames of the objects and adjusted below
    nm = names(newsites)
    rename!(newsites, sites => :sites)
    asm.site.sitestats = join(asm.site.sitestats, newsites, kind = :left, on = :sites, makeunique = makeunique)
    names!(newsites, nm)
    nothing
end

function addsitestats!(asm::Assemblage, newsites::AbstractVector, name::Union{String, Symbol})
    length(newsites) == nsites(asm) || error("Cannot add a vector of length $(length(newsites)) to an Assemblage with $(nsites(asm)) sites")
    asm.site.sitestats[Symbol(name)] = newsites
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

macro traits(x, expr)
    :(@with(traits($x), $expr))
end

macro sitestats(x, expr)
    :(@with(sitestats($x), $expr))
end

dispersionfield(asm::EcoBase.AbstractAssemblage, site) =
    occurrences(asm)' * placeoccurrences(asm, site)
