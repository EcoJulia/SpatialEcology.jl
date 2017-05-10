# the forward macro was copied in from Lazy.jl at the suggestion of @MikeInnes
macro forward(ex, fs)
    @capture(ex, T_.field_) || error("Syntax: @forward T.x f, g, h")
    T = esc(T)
    fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
    :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
        for f in fs]...);
    nothing)
end

@forward Assmbl.occ nspecies, nsites, specnames, sitenames, occupancy, richness, records
@forward AbstractOccFields.commatrix nspecies, nsites, specnames, sitenames, occupancy, richness, records

nspecies(com::AbstractComMatrix) = size(com.occurrences, 2)
nsites(com::AbstractComMatrix) = size(com.occurrences, 1)

nsites(sd::SpatialData) = size(coordinates(sd.site), 1)
nsites(sd::SiteFields) = DataFrames.nrow(sd.sitestats)

specnames(com::AbstractComMatrix) = NamedArrays.names(com.occurrences)[2]

sitenames(com::AbstractComMatrix) = NamedArrays.names(com.occurrences)[1]
sitenames(sd::SpatialData) = sitenames(sd.site)
sitenames(sd::SiteFields) = NamedArrays.names(coordinates(sd))[1]

occupancy{T<:Bool}(com::AbstractComMatrix{T}) = sum(com.occurrences, 1)[1,:]
occupancy{T<:Int}(com::AbstractComMatrix{T}) = mapslices(x->sum(i > 0 for i in x), com.occurrences, 1)[1,:]

richness{T<:Bool}(com::AbstractComMatrix{T}) = sum(com.occurrences, 2)[:,1]
richness{T<:Int}(com::AbstractComMatrix{T}) = mapslices(x->sum(i > 0 for i in x), com.occurrences, 2)[:,1]

records{T<:Int}(com::AbstractComMatrix{T}) = sum(i > 0 for i in com.occurrences)
records{T<:Bool}(com::AbstractComMatrix{T}) = sum(com.occurrences)

size(com::AbstractComMatrix) = size(com.occurrences)
size(com::AbstractComMatrix, dims...) = size(com.occurrences, dims...)

summary(com::AbstractComMatrix) = "$(nsites(com))x$(nspecies(com)) $(typeof(com))"

function createsummaryline{T<:AbstractString}(vec::Vector{T})
    linefunc(vec) = mapreduce(x->x*", ", *, vec[1:(end-1)])*vec[end]
    length(vec) < 6 && return linefunc(vec)
    linefunc(vec[1:3])*"..."*linefunc(vec[(end-1):end])
end

function show(io::IO, com::ComMatrix)
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "ComMatrix with $(records(com)) records of $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, com::Assemblage)
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "Assemblage with $(records(com)) records of $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, sd::SiteData)
    println(io, "Spatial data set with $(nsites(sd)) sites\n\nSite names:\n$(createsummarylines(sitenames(sd)))")
end

#TODO also create render functions for Juno

getindex(com::AbstractComMatrix, inds...) = ComMatrix(getindex(com.occurrences, inds...))

setindex!(com::AbstractComMatrix, X, inds...) = setindex!(com.occurrences, X, inds...)

function getindex{S <: SiteFields}(site::S, inds)
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
