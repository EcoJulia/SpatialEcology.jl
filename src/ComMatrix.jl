
macro forward_func(ex, fs)
    T, field = ex.args[1], ex.args[2].value

    T = esc(T)
    fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
    :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
        for f in fs]...);
    nothing)
end

@forward_func Assmbl.occ nspecies, nsites, occupancy, richness, nrecords, occurring, occupied, specnames
@forward_func AbstractOccFields.commatrix nspecies, nsites, specnames, sitenames, occupancy, richness, nrecords, occurring, occupied
@forward_func Assmbl.site sitenames


#--------------------------------------------------------------------------
# Basic summary functions

occurring(com::AbstractComMatrix) = nzrows(com.occurrences)
occupied(com::AbstractComMatrix) = nzcols(com.occurrences)
occupied(com::AbstractComMatrix, idx) = findall(!iszero, com.occurrences[idx,:])
occurring(com::AbstractComMatrix, idx) = findall(!iszero, com.occurrences[:,idx])

noccurring(x) = length(occurring(x))
noccupied(x) = length(occupied(x))
noccurring(x, idx) = length(occurring(x, idx))
noccupied(x, idx) = length(occupied(x, idx))

nspecies(com::AbstractComMatrix) = size(com.occurrences, 1)
nsites(com::AbstractComMatrix) = size(com.occurrences, 2)

nsites(sd::SpatialData) = size(coordinates(sd.site), 1)
nsites(sd::SiteFields) = DataFrames.ncol(sd.sitestats)

getspecies(com::AbstractComMatrix{T}, idx) where T = view(com.occurrences, idx, :)
getsite(com::AbstractComMatrix{T}, idx) where T = view(com.occurrences, :, idx)

specnames(com::AbstractComMatrix) = com.specnames

sitenames(com::AbstractComMatrix) = com.sitenames
sitenames(sd::SpatialData) = sitenames(sd.site)
sitenames(sd::SiteFields) = collect(sd.sitestats[:sites])

sitetotals(com::AbstractComMatrix) = vec(colsum(com.occurrences))
speciestotals(com::AbstractComMatrix) = vec(rowsum(com.occurrences))

richness(com::AbstractComMatrix{T}) where T<:Bool = vec(colsum(com.occurrences))
richness(com::AbstractComMatrix{T}) where T<:Real = vec(mapslices(nnz, com.occurrences, dims = 1))

occupancy(com::AbstractComMatrix{T}) where T<:Bool = vec(rowsum(com.occurrences))
occupancy(com::AbstractComMatrix{T}) where T<:Real = vec(mapslices(nnz, com.occurrences, dims = 2))

nrecords(com::AbstractComMatrix) = _nnz(occurrences(com))

size(com::AbstractComMatrix) = size(occurrences(com))
size(com::AbstractComMatrix, dims...) = size(occurrences(com), dims...)

"""
    cooccurring(com, inds...)

Ret
"""
cooccurring(com::AbstractComMatrix, inds...) = cooccurring(com, [inds...])
function cooccurring(com::AbstractComMatrix, inds::AbstractVector)
    sub = view(com, species = inds)
    richness(sub) .== nspecies(sub)
end


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
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "ComMatrix with $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, com::Assemblage)
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "Assemblage with $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, sd::SiteData)
    println(io, "Spatial data set with $(nsites(sd)) sites\n\nSite names:\n$(createsummarylines(sitenames(sd)))")
end

#TODO also create render functions for Juno

getindex(com::AbstractComMatrix, inds...) = ComMatrix(getindex(com.occurrences, inds...))

setindex!(com::AbstractComMatrix, X, inds...) = setindex!(com.occurrences, X, inds...)

function getindex(site::S, inds) where S<:SiteFields
  S(coordinates(site)[inds,:], site.sitestats[inds,:])
end

function getindex(com::Assmbl, ind::Symbol)
    if ind in names(com.site.sitestats)
        return com.site.sitestats[ind]
    elseif ind in names(com.occ.traits)
        return com.occ.traits[ind]
    else
        error("No such name in traits or sitestats")
    end
end
