
macro forward_func(ex, fs)
    T, field = ex.args[1], ex.args[2].value

    T = esc(T)
    fs = Meta.isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
    :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
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
nplaces(sd::SiteData) = size(coordinates(sd.site), 1)
nplaces(sd::SELocations) = DataFrames.ncol(sd.sitestats)
nplaces(gr::GridData) = size(gr.indices, 1)
nplaces(pd::PointData) = size(pd.coords, 1)

"""
    nrecords(com)
"""
nrecords(com::AbstractComMatrix) = _nnz(occurrences(com))

const getspecies = thingoccurrences

"""
    getspecies(com)
"""
thingoccurrences(com::AbstractComMatrix, idx) = thingoccurrences(occurrences(com), asindices(idx, thingnames(com)))

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
placenames(sd::SiteData) = sitenames(sd.site)
placenames(sd::SELocations) = sd.sitestats[:,:sites]

"""
    sitetotals(com)
"""
sitetotals(com::AbstractComMatrix) = vec(colsum(com.occurrences))

"""
    speciestotals(com)
"""
speciestotals(com::AbstractComMatrix) = vec(rowsum(com.occurrences))

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

function show(io::IO, sd::SiteData)
    println(io, "Spatial data set with $(nsites(sd)) sites\n\nSite names:\n$(createsummaryline(sitenames(sd)))")
end

#TODO also create render functions for Juno

getindex(com::AbstractComMatrix, inds...) = ComMatrix(getindex(com.occurrences, inds...))

setindex!(com::AbstractComMatrix, X, inds...) = setindex!(com.occurrences, X, inds...)

function getindex(site::S, inds) where S<:SELocations
  S(coordinates(site)[inds,:], site.sitestats[inds,:])
end

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
