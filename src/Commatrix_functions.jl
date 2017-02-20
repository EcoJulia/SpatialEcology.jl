
# the forward macro was copied in from Lazy.jl at the suggestion of @MikeInnes
macro forward(ex, fs)
  @capture(ex, T_.field_) || error("Syntax: @forward T.x f, g, h")
  T = esc(T)
  fs = isexpr(fs, :tuple) ? map(esc, fs.args) : [esc(fs)]
  :($([:($f(x::$T, args...) = (Base.@_inline_meta; $f(x.$field, args...)))
       for f in fs]...);
    nothing)
end

#TODO apply this everywhere below

nspecies(com::ComMatrix) = size(com.occurrences, 2)
nspecies(ocf::OccFields) = nspecies(ocf.commatrix)
nspecies(asm::Assmbl) = nspecies(asm.occ)



nsites(com::ComMatrix) = size(com.occurrences, 1)

nsites(ocf::OccFields) = nsites(ocf.commatrix)
nsites(asm::Assmbl) = nsites(asm.occ)
nsites(sd::SpatialData) = size(coordinates(sd.site), 1)
nsites(sd::SiteFields) = DataFrames.nrow(sd.sitestats)


specnames(com::ComMatrix) = NamedArrays.names(com.occurrences)[2]
specnames(ocf::OccFields) = specnames(ocf.commatrix)
specnames(asm::Assmbl) = specnames(asm.occ)


sitenames(com::ComMatrix) = NamedArrays.names(com.occurrences)[1]
sitenames(ocf::OccFields) = sitenames(ocf.commatrix)
sitenames(asm::Assmbl) = sitenames(asm.occ)
sitenames(sd::SpatialData) = sitenames(sd.site)
sitenames(sd::SiteFields) = NamedArrays.names(coordinates(sd))[1]

occupancy{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences, 1)[1,:]
occupancy{T<:Int}(com::ComMatrix{T}) = mapslices(x->sum(i > 0 for i in x), com.occurrences, 1)[1,:]
occupancy(ocf::OccFields) = occupancy(ocf.commatrix)
occupancy(asm::Assmbl) = occupancy(asm.occ)



richness{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences, 2)[:,1]
richness{T<:Int}(com::ComMatrix{T}) = mapslices(x->sum(i > 0 for i in x), com.occurrences, 2)[:,1]
richness(ocf::OccFields) = richness(ocf.commatrix)
richness(asm::Assmbl) = richness(asm.occ)



records{T<:Int}(com::ComMatrix{T}) = sum(i > 0 for i in com.occurrences)
records{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences)
records(ocf::OccFields) = records(ocf.commatrix)
records(asm::Assmbl) = records(asm.occ)



size(com::ComMatrix) = size(com.occurrences)
size(com::ComMatrix, dims...) = size(com.occurrences, dims...)



summary(com::ComMatrix) = "$(nsites(com))x$(nspecies(com)) $(typeof(com))"



function createsummaryline{T<:AbstractString}(vec::Vector{T})
    linefunc(vec) = mapreduce(x->x*", ", *, vec[1:(end-1)])*vec[end]
    length(vec) < 6 && return linefunc(vec)
    linefunc(vec[1:3])*"..."*linefunc(vec[(end-1):end])
end



function show(io::IO, com::ComMatrix)
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println("Community matrix with $(records(com)) records of $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, com::Assemblage)
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println("Assemblage with $(records(com)) records of $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, sd::SiteData)
    println("Spatial data set with $(nsites(sd)) sites\n\nSite names:\n$(createsummarylines(sitenames(sd)))")
end

# why did I do this?
#function getindex(com::ComMatrix, inds...)
#    com.occurrences = getindex(com.occurrences, inds...)
#    com
#end

getindex(com::ComMatrix, inds...) = getindex(com.occurrences, inds...)

setindex!(com::ComMatrix, X, inds...) = setindex!(com.occurrences, X, inds...)


function getindex{S <: SiteFields}(site::S, inds)
  S(coordinates(site)[inds,:], site.sitestats[inds,:])
end
