

nspecies(com::ComMatrix) = size(com.occurrences, 2)
nsites(com::ComMatrix) = size(com.occurrences, 1)
specnames(com::ComMatrix) = NamedArrays.allnames(com.occurrences)[2]
sitenames(com::ComMatrix) = NamedArrays.allnames(com.occurrences)[1]
occupancy{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences, 1)
occupancy{T<:Int}(com::ComMatrix{T}) = sum(com.occurrences .> 0, 1)
richness{T<:Bool}(com::ComMatrix{T}) = sum(com.occurrences, 2)
richness{T<:Int}(com::ComMatrix{T}) = sum(com.occurrences .> 0, 2)


#TODO specify getindex and setindex variations that just pass on to occurrences

getindex(com::ComMatrix, i) = getindex(com.occurrences, i)
getindex(com::ComMatrix, i, j) = getindex(com.occurrences, i, j)
setindex(com::ComMatrix, i) = setindex(com.occurrences, i)
setindex(com::ComMatrix, i, j) = setindex(com.occurrences, i, j)
