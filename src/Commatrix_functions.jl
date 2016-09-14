

nspecies(com::ComMatrix) = size(com.occurrences, 2)
nsites(com::ComMatrix) = size(com.occurrences, 1)
occupancy{T<:Bool}(com::ComMatrix{T}) = sum(com.occurences, 2)
occupancy{T<:Int}(com::ComMatrix{T}) = sum(com.occurences .> 0, 2)
richness{T<:Bool}(com::ComMatrix{T}) = sum(com.occurences, 1)
richness{T<:Int}(com::ComMatrix{T}) = sum(com.occurences .> 0, 1)
