

nspecies(com::ComMatrix) = size(com.occurrences, 2)
nsites(com::ComMatrix) = size(com.occurrences, 1)
specnames(com::ComMatrix) = NamedArrays.allnames(com.occurrences)[2]
sitenames(com::ComMatrix) = NamedArrays.allnames(com.occurrences)[1]
occupancy{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences, 1)[1,:]
occupancy{T<:Int}(com::ComMatrix{T}) = sum(com.occurrences .> 0, 1)[1,:]
richness{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences, 2)[:,1]
richness{T<:Int}(com::ComMatrix{T}) = sum(com.occurrences .> 0, 2)[:,1]
records{T<:Int}(com::ComMatrix{T}) = sum(com.occurrences .> 0)
records{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences)

function getindex(com::ComMatrix, inds...)
    com.occurrences = getindex(com.occurrences, inds...)
    com
end

setindex!(com::ComMatrix, X, inds...) = setindex!(com.occurrences, X, inds...)

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

richness(ocf::OccFields) = richness(ocf.commatrix)
occupancy(ocf::OccFields) = occupancy(ocf.commatrix)
sitenames(ocf::OccFields) = sitenames(ocf.commatrix)
specnames(ocf::OccFields) = specnames(ocf.commatrix)
nsites(ocf::OccFields) = nsites(ocf.commatrix)
nspecies(ocf::OccFields) = nspecies(ocf.commatrix)
records(ocf::OccFields) = records(ocf.commatrix)

richness(asm::Assmbl) = richness(asm.occ)
occupancy(asm::Assmbl) = occupancy(asm.occ)
sitenames(asm::Assmbl) = sitenames(asm.occ)
specnames(asm::Assmbl) = specnames(asm.occ)
nsites(asm::Assmbl) = nsites(asm.occ)
nspecies(asm::Assmbl) = nspecies(asm.occ)
records(asm::Assmbl) = records(asm.occ)


function show(io::IO, com::Assemblage)
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println("Assemblage with $(records(com)) records of $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end


nsites(sd::SpatialData) = size(coordinates(sd.site), 1)
sitenames(sd::SpatialData) = NamedArrays.allnames(coordinates(sd.site))[1]

function show(io::IO, sd::SiteData)
    println("Spatial data set with $(nsites(sd)) sites\n\nSite names:\n$(createsummarylines(sitenames(sd)))")
end
