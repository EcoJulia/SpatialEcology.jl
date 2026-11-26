
macro forward_func(ex, fs)
    T, field = ex.args[1], ex.args[2].value

    T = esc(T)
    fs = Meta.isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
    :($([:(@inline $f(x::$T, args...) = $f(x.$field, args...))
        for f in fs]...);
    nothing)
end

@forward_func SEAssemblage.occ nthings, nplaces, occupancy, richness, nrecords, thingnames
@forward_func SEThings.commatrix nthings, nplaces, thingnames, placenames, occupancy, richness, nrecords, occurring, occupied
@forward_func SEAssemblage.site placenames


#--------------------------------------------------------------------------
# Basic summary functions

"""
    occupancy(com::AbstractComMatrix)
"""
occupancy(com::AbstractComMatrix) = occupancy(occurrences(com))

"""
    richness(com::AbstractComMatrix)
"""
richness(com::AbstractComMatrix) = richness(occurrences(com))

"""
    occurring(com::AbstractComMatrix)
"""
occurring(com::AbstractComMatrix) = occurring(occurrences(com))
occurring(com::AbstractComMatrix, idx) = occurring(occurrences(com), asindices(idx, sitenames(com)))

"""
    occupied(com::AbstractComMatrix)
"""
occupied(com::AbstractComMatrix) = occupied(occurrences(com))
occupied(com::AbstractComMatrix, idx) = occupied(occurrences(com), asindices(idx, speciesnames(com)))

# Fast per-species site lookup via the transpose kept by the OccMatrix storage.
#
# The occurrence matrix is n_species × n_sites CSC — column access (site queries)
# is O(nnz_per_site), but row access (species queries) requires scanning all
# columns: O(nnz_total). The transpose (n_sites × n_species CSC) flips this:
# per-species site lookup is now a CSC column read, O(nnz_per_species).
#
# Only `ComMatrix` keeps a transpose; `SubComMatrix` keeps the slow fallback.
function occupied(com::ComMatrix, idx::Integer)
    t = transposed(com.occurrences)
    rowvals(t)[nzrange(t, idx)]
end
occupied(com::ComMatrix, idx::AbstractVector{<:Integer}) =
    nzrows(view(transposed(com.occurrences), :, idx))
occupied(com::ComMatrix, idx) = occupied(com, asindices(idx, speciesnames(com)))

const nspecies = nthings

"""
    nthings(com::AbstractComMatrix)
"""
nthings(com::AbstractComMatrix) = size(com.occurrences, 1)

const nsites = nplaces

"""
    nsites(com)
"""
nplaces(com::AbstractComMatrix) = size(com.occurrences, 2)
nplaces(sd::SELocations) = DataFrames.nrow(sd.sitestats)
nplaces(gr::GridData) = size(gr.indices, 1)
nplaces(pd::PointData) = size(pd.coords, 1)
nplaces(rd::RasterData) = length(rd.cellinds)  # SubRasterData handled in Gridfunctions.jl

"""
    nrecords(com)
"""
nrecords(com::AbstractComMatrix) = _nnz(occurrences(com))

const getspecies = thingoccurrences

"""
    getspecies(com)
"""
thingoccurrences(com::AbstractComMatrix, idx) = thingoccurrences(occurrences(com), asindices(idx, thingnames(com)))

# Fast single-species occurrences: a column of the transpose is the species' row of
# the occurrences, read in O(nnz_per_species) via the CSC column pointer rather than
# O(nnz_total) for a CSC row scan. It is returned as a copy (a SparseVector): a view
# into the transpose would go stale after the next write, and writing through it
# would change the transpose but not the occurrences.
thingoccurrences(com::ComMatrix, idx::Integer) = transposed(com.occurrences)[:, idx]
thingoccurrences(com::ComMatrix, idx) =
    thingoccurrences(com, asindices(idx, thingnames(com)))

const getsite = placeoccurrences

"""
    getsite(com)
"""
placeoccurrences(com::AbstractComMatrix, idx) = placeoccurrences(occurrences(com), asindices(idx, placenames(com)))

const speciesnames = thingnames

"""
    speciesnames(com)
"""
thingnames(com::AbstractComMatrix) = com.speciesnames

const sitenames = placenames

"""
    sitenames(com)
"""
placenames(com::AbstractComMatrix) = com.sitenames
placenames(sd::SELocations) = sd.sitestats[:,:sites]

"""
    sitetotals(com)
"""
sitetotals(com::AbstractComMatrix) = vec(colsum(com.occurrences))

"""
    speciestotals(com)
"""
speciestotals(com::AbstractComMatrix) = vec(rowsum(com.occurrences))

# Fast range-size count without allocating the sites vector: O(1) per species.
noccupied(com::ComMatrix, idx::Integer) = length(nzrange(transposed(com.occurrences), idx))
noccupied(com::ComMatrix, idx) = noccupied(com, asindices(idx, speciesnames(com)))

# For Bool assemblages, per-species counts are the colptr differences of the
# transpose (O(n_species) read), not a rowval scan of the occurrences (O(nnz_total)).
# This counts stored entries, which is right because OccMatrix never stores zeros.
occupancy(com::ComMatrix{Bool}) =
    Vector{Int}(diff(SparseArrays.getcolptr(transposed(com.occurrences))))
speciestotals(com::ComMatrix{Bool}) = occupancy(com)

size(com::AbstractComMatrix) = size(occurrences(com))
size(com::AbstractComMatrix, dims...) = size(occurrences(com), dims...)


#----------------------------------------------------------------------------------------------
# Mutating and concatenating




#-------------------------------------------------------------------------------------------------
# show methods

summary(com::AbstractComMatrix) = "$(nsites(com))x$(nspecies(com)) $(typeof(com))"

function createsummaryline(vec::AbstractVector{T}) where T<:AbstractString
    linefunc(vec) = mapreduce(x->x*", ", *, vec[1:(end-1)])*vec[end]
    length(vec) == 1 && return vec[1]
    length(vec) < 6 && return linefunc(vec)
    linefunc(vec[1:3])*"..."*linefunc(vec[(end-1):end])
end

function show(io::IO, com::ComMatrix)
    sp = createsummaryline(speciesnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "ComMatrix with $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, com::Assemblage)
    sp = createsummaryline(speciesnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "Assemblage with $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

#TODO also create render functions for Juno

#-------------------------------------------------------------------------------------------------
# Indexing

# Turn a species or site selector into positions. Positions, ranges and `:` pass
# through, a Bool mask selects where it is true, and names (Strings or Symbols,
# single or in a vector) are looked up in `names`. Unlike `asindices`, an unknown
# name is an error rather than silently skipped.
_toindices(::Colon, names) = Base.OneTo(length(names))
_toindices(i::Integer, names) = i
_toindices(I::AbstractVector{<:Integer}, names) = I
function _toindices(I::AbstractVector{Bool}, names)
    length(I) == length(names) ||
        throw(BoundsError(names, I))
    findall(I)
end
_toindices(I::Union{AbstractString, Symbol}, names) = (_checknames((I,), names); asindices(I, names))
_toindices(I::AbstractVector{<:Union{Missing, AbstractString, Symbol}}, names) = (_checknames(I, names); asindices(I, names))

# A `missing` is not a name either.
function _checknames(I, names)
    known = Set(String.(names))
    unknown = [n === missing ? "missing" : repr(String(n)) for n in I if n === missing || !(String(n) in known)]
    isempty(unknown) ||
        throw(ArgumentError("no species or site named $(join(unknown, ", "))"))
    nothing
end

# A single position keeps its dimension: com[1:3, 1] is a 3×1 ComMatrix.
_keepdim(i::Integer) = i:i
_keepdim(I) = I

"""
    com[species, sites]

Index a `ComMatrix` (or `SubComMatrix`) like a matrix with species as rows and
sites as columns. Each index can be a position, a range or vector of positions,
`:`, a Bool mask, or species/site names (Strings or Symbols), alone or in a vector.

A single species in a single site returns the value. Any other selection returns a
new `ComMatrix` with the selected species and sites and their names, so
`com[1:3, 1]` is a 3×1 `ComMatrix`; for a plain vector use `getspecies`/`getsite`.
Use `view` (or `@view com[species, sites]`) to get a `SubComMatrix` instead of a copy.
An unknown name throws an `ArgumentError`.
"""
function getindex(com::AbstractComMatrix, I, J)
    i = _toindices(I, speciesnames(com))
    j = _toindices(J, sitenames(com))
    _getindex(com, i, j)
end
_getindex(com::AbstractComMatrix, i::Integer, j::Integer) = occurrences(com)[i, j]
function _getindex(com::AbstractComMatrix, i, j)
    ii, jj = _keepdim(i), _keepdim(j)
    occ = _slice(occurrences(com), ii, jj)
    ComMatrix{eltype(occ)}(occ, speciesnames(com)[ii], sitenames(com)[jj])
end

# Slice the CSC matrix itself, so the result is a fresh SparseMatrixCSC built by
# SparseArrays' sparse getindex; for a view, index its parent at the view's positions.
_slice(W::OccMatrix, ii, jj) = W.csc[ii, jj]
function _slice(v::SubArray{<:Any, 2, <:OccMatrix}, ii, jj)
    pi = parentindices(v)
    parent(v).csc[pi[1][ii], pi[2][jj]]
end

"""
    com[species, sites] = x

Write into a `ComMatrix`, or through a `SubComMatrix` into its parent. Indices are as
for `getindex`. Writing a zero removes the entry, so it is not counted as a presence.
"""
function setindex!(com::AbstractComMatrix, X, I, J)
    setindex!(occurrences(com), X, _toindices(I, speciesnames(com)), _toindices(J, sitenames(com)))
    com
end

Base.axes(com::AbstractComMatrix) = axes(occurrences(com))
Base.axes(com::AbstractComMatrix, d) = axes(occurrences(com), d)
Base.firstindex(com::AbstractComMatrix, d) = 1
Base.lastindex(com::AbstractComMatrix, d) = size(com, d)

getindex(site::SELocations, inds) = copy(view(site, inds))

function getindex(com::SEAssemblage, ind::Symbol)
    if ind in propertynames(com.site.sitestats)
        return com.site.sitestats[:,ind]
    elseif ind in propertynames(com.occ.traits)
        return com.occ.traits[:,ind]
    else
        error("No such name in traits or sitestats")
    end
end

getindex(com::SEAssemblage, ::Colon, ind::Symbol) = getindex(com, ind)

function getindex(com::SEAssemblage, ::typeof(!), ind::Symbol)
    if ind in propertynames(com.site.sitestats)
        return com.site.sitestats[!,ind]
    elseif ind in propertynames(com.occ.traits)
        return com.occ.traits[!,ind]
    else
        error("No such name in traits or sitestats")
    end
end



# Extend Distances interface for pairwise distances

pairwise(m::PreMetric, com::AbstractComMatrix; dims=2) = pairwise(m, collect(occurrences(com)), dims=dims)
