

nspecies(com::ComMatrix) = size(com.occurrences, 2)
nsites(com::ComMatrix) = size(com.occurrences, 1)
specnames(com::ComMatrix) = NamedArrays.allnames(com.occurrences)[2]
sitenames(com::ComMatrix) = NamedArrays.allnames(com.occurrences)[1]
occupancy{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences, 1)
occupancy{T<:Int}(com::ComMatrix{T}) = sum(com.occurrences .> 0, 1)
richness{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences, 2)
richness{T<:Int}(com::ComMatrix{T}) = sum(com.occurrences .> 0, 2)

function getindex(com::ComMatrix, inds...)
    com.occurrences = getindex(com.occurrences, inds...)
    com
end

setindex!(com::ComMatrix, X, inds...) = setindex!(com.occurrences, X, inds...)

size(com::ComMatrix) = size(com.occurrences)
size(com::ComMatrix, dims...) = size(com.occurrences, dims...)
