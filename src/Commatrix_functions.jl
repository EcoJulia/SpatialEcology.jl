# the forward macro was copied in from Lazy.jl at the suggestion of @MikeInnes
macro forward(ex, fs)
    T, field = ex.args[1], ex.args[2].args[1]
    T = esc(T)
    fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
    :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
        for f in fs]...);
    nothing)
end

@my_forward Assmbl.occ nspecies, nsites, occupancy, richness, records, occurring, occupied, specnames, sitenames
@my_forward AbstractOccFields.commatrix nspecies, nsites, specnames, sitenames, occupancy, richness, records, occurring, occupied
@my_forward Assmbl.site sitenames

occurring(com::AbstractComMatrix) = nzcols(com.occurrences)
occupied(com::AbstractComMatrix) = nzrows(com.occurrences)
occurring(com::AbstractComMatrix{T}, idx) where T<:Bool = find(com.occurrences[idx,:])
occupied(com::AbstractComMatrix{T}, idx) where T<:Bool = find(com.occurrences[:,idx])
occurring(com::AbstractComMatrix{T}, idx) where T<:Int = find(com.occurrences[idx,:] .> 0)
occupied(com::AbstractComMatrix{T}, idx) where T<:Int = find(com.occurrences[:,idx] .> 0)

noccurring(x) = length(occurring(x))
noccupied(x) = length(occupied(x))
noccurring(x, idx) = length(occurring(x, idx))
noccupied(x, idx) = length(occupied(x, idx))

nspecies(com::AbstractComMatrix) = size(com.occurrences, 2)
nsites(com::AbstractComMatrix) = size(com.occurrences, 1)

nsites(sd::SpatialData) = size(coordinates(sd.site), 1)
nsites(sd::SiteFields) = DataFrames.nrow(sd.sitestats)

specnames(com::AbstractComMatrix) = com.specnames

sitenames(com::AbstractComMatrix) = com.sitenames
sitenames(sd::SpatialData) = sitenames(sd.site)
sitenames(sd::SiteFields) = collect(sd.sitestats[:sites])

occupancy(com::AbstractComMatrix{T}) where T<:Bool = vec(colsum(com.occurrences))
occupancy(com::AbstractComMatrix{T}) where T<:Int = vec(mapslices(x->sum(i > 0 for i in x), com.occurrences, 1))

richness(com::AbstractComMatrix{T}) where T<:Bool = vec(rowsum(com.occurrences))
richness(com::AbstractComMatrix{T}) where T<:Int = vec(mapslices(x->sum(i > 0 for i in x), com.occurrences, 2))

records(com::AbstractComMatrix{T}) where T<:Int = sum(i > 0 for i in com.occurrences)
records(com::AbstractComMatrix{T}) where T<:Bool = sum(com.occurrences)

size(com::AbstractComMatrix) = size(com.occurrences)
size(com::AbstractComMatrix, dims...) = size(com.occurrences, dims...)

summary(com::AbstractComMatrix) = "$(nsites(com))x$(nspecies(com)) $(typeof(com))"

function createsummaryline(vec::AbstractVector{T}) where T<:AbstractString
    linefunc(vec) = mapreduce(x->x*", ", *, vec[1:(end-1)])*vec[end]
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
