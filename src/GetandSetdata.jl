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

The species × sites occurrence matrix of an `Assemblage` or `ComMatrix`, as the
sparse `OccMatrix` the `ComMatrix` stores (a view of it for `SubComMatrix`), not a
copy. Writing into it with indexing or broadcasting keeps the `ComMatrix`
consistent; see `SpatialEcology.OccMatrix`.
"""
occurrences(asm::Union{SEAssemblage, SEThings}) = occurrences(commatrix(asm))
occurrences(cm::AbstractComMatrix) = cm.occurrences

#--------------------------------------------------------------------------
# Renaming species and sites (EcoBase's thingnames!/placenames!)
#
# An Assemblage holds each set of names twice: species in the ComMatrix and in the
# `name` column of the traits, sites in the ComMatrix and in the `sites` column of
# the sitestats. Renaming changes both, and gives each its own vector.

const speciesnames! = thingnames!
const sitenames! = placenames!

# The new names as a fresh Vector{String}, after checking that there is one per
# species (or site), that none is missing and that no two are alike. Every
# rename checks here before writing anything, so a failed one changes nothing.
function _newnames(newnames::AbstractVector, n::Integer, what::AbstractString)
    length(newnames) == n ||
        throw(ArgumentError("expected $n $what names, got $(length(newnames))"))
    any(ismissing, newnames) && throw(ArgumentError("$what names cannot be missing"))
    names = String[string(x) for x in newnames]
    seen, dups = Set{String}(), String[]
    for nm in names
        nm in seen ? push!(dups, nm) : push!(seen, nm)
    end
    isempty(dups) ||
        throw(ArgumentError("$what names must be unique, but $(join(repr.(unique(dups)), ", ")) appear more than once"))
    names
end

# The ComMatrix's vectors are overwritten in place, not replaced, because a
# SubComMatrix holds views of them: views taken earlier then show the new names.
function thingnames!(com::ComMatrix, newnames::AbstractVector)
    com.speciesnames .= _newnames(newnames, nspecies(com), "species")
    com
end

function placenames!(com::ComMatrix, newnames::AbstractVector)
    com.sitenames .= _newnames(newnames, nsites(com), "site")
    com
end

# The DataFrame columns are replaced instead. A SubDataFrame reads its parent's
# columns when asked, so views still see the new names, and a `sites` column of
# another type (an Int one, say) becomes a String one.
function thingnames!(occ::SpeciesData, newnames::AbstractVector)
    names = _newnames(newnames, nspecies(occ), "species")
    thingnames!(occ.commatrix, names)      # copied into the ComMatrix's own vector
    occ.traits[!, :name] = names
    occ
end

placenames!(occ::SpeciesData, newnames::AbstractVector) = (placenames!(occ.commatrix, newnames); occ)

function placenames!(lo::Locations, newnames::AbstractVector)
    lo.sitestats[!, :sites] = _newnames(newnames, nsites(lo), "site")
    lo
end

"""
    speciesnames!(asm, newnames)
    speciesnames!(asm, old => new, ...)
    speciesnames!(asm, dict)
    speciesnames!(f, asm)

Rename the species of an `Assemblage` (or a `ComMatrix`) in place, and return it.

`newnames` gives a new name for every species, in order. `old => new` pairs, or a
`Dict` of them, rename only the species they name; all are looked up before any is
renamed, so `"a" => "b", "b" => "a"` swaps two names. A function `f` is applied to
every name. Names are stored as `String`s, must be unique and cannot be `missing`;
if they are not, nothing is renamed.

The names in the community matrix and the `name` column of the traits change
together, and views of `asm` taken earlier show the new names. A view cannot be
renamed itself: rename its parent.

```julia
speciesnames!(asm, "sp1" => "Abies alba")
speciesnames!(n -> replace(n, "_" => " "), asm)
```
"""
thingnames!(asm::Assemblage, newnames::AbstractVector) = (thingnames!(asm.occ, newnames); asm)

"""
    sitenames!(asm, newnames)
    sitenames!(asm, old => new, ...)
    sitenames!(asm, dict)
    sitenames!(f, asm)

Rename the sites of an `Assemblage` (or a `ComMatrix`) in place, and return it.
The forms and rules are those of [`speciesnames!`](@ref). The names in the
community matrix and the `sites` column of the sitestats change together.

```julia
sitenames!(s -> "plot_" * s, asm)
```
"""
function placenames!(asm::Assemblage, newnames::AbstractVector)
    _newnames(newnames, nsites(asm), "site")       # check before changing either copy
    placenames!(asm.occ, newnames)
    placenames!(asm.site, newnames)
    asm
end

const _HasThingNames = Union{ComMatrix, SpeciesData, Assemblage}
const _HasPlaceNames = Union{ComMatrix, SpeciesData, Locations, Assemblage}

thingnames!(f::Function, x::_HasThingNames) = thingnames!(x, map(f, _currentnames(thingnames(x))))
placenames!(f::Function, x::_HasPlaceNames) = placenames!(x, map(f, _currentnames(placenames(x))))

thingnames!(x::_HasThingNames, p::Pair, ps::Pair...) = thingnames!(x, _renamed(thingnames(x), (p, ps...), "species"))
placenames!(x::_HasPlaceNames, p::Pair, ps::Pair...) = placenames!(x, _renamed(placenames(x), (p, ps...), "site"))
thingnames!(x::_HasThingNames, d::AbstractDict) = thingnames!(x, _renamed(thingnames(x), d, "species"))
placenames!(x::_HasPlaceNames, d::AbstractDict) = placenames!(x, _renamed(placenames(x), d, "site"))

# A user-supplied `sites` column need not hold Strings.
_currentnames(names) = String[string(n) for n in names]

# The names after renaming each `old => new` in `ps`. Every old name is looked up
# among the current names before any is changed, so a swap works.
function _renamed(current, ps, what::AbstractString)
    current = _currentnames(current)
    names = collect(Any, current)      # `_newnames` checks and converts the new ones
    for (old, new) in ps
        old isa Union{AbstractString, Symbol} ||
            throw(ArgumentError("rename $what by name, as old => new; got $(repr(old)) => $(repr(new))"))
        i = findfirst(==(String(old)), current)
        i === nothing && throw(ArgumentError("no $what named $(repr(String(old)))"))
        names[i] = new
    end
    names
end

# Views refuse, so a rename cannot reach the parent's names behind its back.
const _ThingNameViews = Union{SubComMatrix, SubSpeciesData, SubAssemblage}
const _PlaceNameViews = Union{SubComMatrix, SubSpeciesData, SubLocations, SubAssemblage}

_noviewrename(what) = throw(ArgumentError("cannot rename the $what of a view; rename them in its parent instead (pass old => new pairs to rename only some)"))
thingnames!(::_ThingNameViews, args...) = _noviewrename("species")
thingnames!(::Function, ::_ThingNameViews) = _noviewrename("species")
placenames!(::_PlaceNameViews, args...) = _noviewrename("sites")
placenames!(::Function, ::_PlaceNameViews) = _noviewrename("sites")

# Left-joins `new` onto `old` by `key`, keeping the row order of `old`. Columns of `old`
# that `new` also has are replaced by the new ones, unless `makeunique` asks for suffixed
# copies instead. Returns a new DataFrame and leaves `old` untouched, so a join that
# throws changes nothing.
function replacingjoin(old::DataFrames.DataFrame, new::DataFrames.DataFrame, key::Symbol, what::String; makeunique = false)
    replaced = makeunique ? String[] : setdiff(intersect(names(old), names(new)), [string(key)])
    isempty(replaced) || @info "Replacing existing $what column$(length(replaced) == 1 ? "" : "s") $(join(replaced, ", ")) (pass `makeunique = true` to keep both)"
    left = select(old, Not(replaced), copycols = false)
    left.__run_number = 1:nrow(left)
    joined = leftjoin(left, new, on = key, makeunique = makeunique)
    sort!(joined, :__run_number)
    select!(joined, Not(:__run_number))
end

"""
    addtraits!(asm::Assemblage, newtraits::DataFrames.DataFrame, species::Symbol; validate = true, tolerance = 0.5, makeunique = false)

Left-join `newtraits` onto the traits of `asm`, matching its `species` column to the species names.
Trait columns that `asm` already has are replaced by the new ones, so the call can be repeated;
pass `makeunique = true` to keep both, with the new columns suffixed.
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
    asm.occ.traits = replacingjoin(asm.occ.traits, jointraits, :name, "trait", makeunique = makeunique)
    nothing
end

function addtraits!(asm::Assemblage, newtraits::AbstractVector, name::Union{String, Symbol})
    length(newtraits) == nspecies(asm) || error("Cannot add a vector of length $(length(newtraits)) to an Assemblage with $(nspecies(asm)) species")
    asm.occ.traits[!,Symbol(name)] = newtraits
end

"""
    addsitestats!(asm::Assemblage, newsites::DataFrames.DataFrame, sites::Symbol; validate = true, tolerance = 0.5, makeunique = false)

Left-join `newsites` onto the site statistics of `asm`, matching its `sites` column to the site names.
Columns that `asm` already has are replaced by the new ones, so the call can be repeated;
pass `makeunique = true` to keep both, with the new columns suffixed.
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
    asm.site.sitestats = replacingjoin(asm.site.sitestats, joinsites, :sites, "site statistic", makeunique = makeunique)
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
