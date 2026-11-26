
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
# (a ComMatrix is an AbstractMatrix, for which EcoBase defines vector selectors too)
occurring(com::AbstractComMatrix, idx::AbstractVector) = occurring(occurrences(com), asindices(idx, sitenames(com)))

"""
    occupied(com::AbstractComMatrix)
"""
occupied(com::AbstractComMatrix) = occupied(occurrences(com))
occupied(com::AbstractComMatrix, idx) = occupied(occurrences(com), asindices(idx, speciesnames(com)))
occupied(com::AbstractComMatrix, idx::AbstractVector) = occupied(occurrences(com), asindices(idx, speciesnames(com)))

# Fast per-species site lookup via the transpose kept by the OccMatrix storage.
#
# The occurrence matrix is n_species × n_sites CSC — column access (site queries)
# is O(nnz_per_site), but row access (species queries) requires scanning all
# columns: O(nnz_total). The transpose (n_sites × n_species CSC) flips this:
# per-species site lookup is now a CSC column read, O(nnz_per_species).
#
# Only a ComMatrix with its own storage keeps a transpose; a `SubComMatrix` (view)
# keeps the slow fallback.
function occupied(com::StoredComMatrix, idx::Integer)
    t = transposed(parent(com))
    rowvals(t)[nzrange(t, idx)]
end
occupied(com::StoredComMatrix, idx::AbstractVector{<:Integer}) =
    nzrows(view(transposed(parent(com)), :, idx))
occupied(com::StoredComMatrix, idx) = occupied(com, asindices(idx, speciesnames(com)))
occupied(com::StoredComMatrix, idx::AbstractVector) = occupied(com, asindices(idx, speciesnames(com)))

const nspecies = nthings

"""
    nthings(com::AbstractComMatrix)
"""
nthings(com::AbstractComMatrix) = size(parent(com), 1)

const nsites = nplaces

"""
    nsites(com)
"""
nplaces(com::AbstractComMatrix) = size(parent(com), 2)
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
thingoccurrences(com::StoredComMatrix, idx::Integer) = transposed(parent(com))[:, idx]
thingoccurrences(com::StoredComMatrix, idx) =
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
thingnames(com::AbstractComMatrix) = parent(DD.lookup(com, Species))

const sitenames = placenames

"""
    sitenames(com)
"""
placenames(com::AbstractComMatrix) = parent(DD.lookup(com, Sites))
placenames(sd::SELocations) = sd.sitestats[:,:sites]

"""
    sitetotals(com)
"""
sitetotals(com::AbstractComMatrix) = vec(colsum(parent(com)))

"""
    speciestotals(com)
"""
speciestotals(com::AbstractComMatrix) = vec(rowsum(parent(com)))

# Fast range-size count without allocating the sites vector: O(1) per species.
noccupied(com::StoredComMatrix, idx::Integer) = length(nzrange(transposed(parent(com)), idx))
noccupied(com::StoredComMatrix, idx) = noccupied(com, asindices(idx, speciesnames(com)))

# For Bool assemblages, per-species counts are the colptr differences of the
# transpose (O(n_species) read), not a rowval scan of the occurrences (O(nnz_total)).
# This counts stored entries, which is right because OccMatrix never stores zeros.
occupancy(com::StoredComMatrix{Bool}) =
    Vector{Int}(diff(SparseArrays.getcolptr(transposed(parent(com)))))
speciestotals(com::StoredComMatrix{Bool}) = occupancy(com)


#----------------------------------------------------------------------------------------------
# Mutating and concatenating




#-------------------------------------------------------------------------------------------------
# show methods

summary(com::AbstractComMatrix) =
    "$(nspecies(com))×$(nsites(com)) $(com isa SubComMatrix ? "SubComMatrix" : "ComMatrix"){$(eltype(com))}"

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
#
# Indexing, views and broadcasting come from DimensionalData (see `ComMatrix`): a
# single value is returned as is, other selections are ComMatrices (or DimVectors
# when a dimension is dropped), and names are selected with `At`. Writes go to the
# OccMatrix storage, which keeps itself consistent.

# The keyword `view(com; species, sites)` keeps both dimensions: a single species or
# site gives a one-row or one-column SubComMatrix.
_keepdim(i::Integer) = i:i
_keepdim(I) = I

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
