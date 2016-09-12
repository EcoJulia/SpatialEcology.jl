
typealias OccMatrix{T} NamedArrays.NamedArray{T,2,SparseMatrixCSC{T,Int},Tuple{Dict{String,Int},Dict{String,Int}}}  #should not be using inheritance - instead use ownership and overload needed features
typealias AbundanceMatrix OccMatrix{Int} #Maybe have an inner constructor that prevents going under zero
typealias PAMatrix OccMatrix{Bool}

@enum coordstype grid points

abstract Assmbl

type Assemblage{T} <: Assmbl # A type to keep subtypes together, ensuring that they are all aligned at all times
    coords::Matrix{Float64}  #This should be spatialpoints - not yet implemented? The only one it MUST have
    cdtype::coordstype
    Occurrences::Nullable(OccMatrix{T})
    traits::Nullable(DataFrames.DataFrame)
    sitestats::Nullable(DataFrames.DataFrame)
end

type PhyloAssemblage{T} <: Assmbl # A type to keep subtypes together, ensuring that they are all aligned at all times
    coords::Matrix{Float64}  #This should be spatialpoints - not yet implemented? The only one it MUST have
    cdtype::coordstype
    Occurrences::Nullable(OccMatrix{T})
    traits::Nullable(DataFrames.DataFrame)
    sitestats::Nullable(DataFrames.DataFrame)
    phylo::Phylo.Phylogeny #This makes it a special type, because so many functions are only defined when there is a phylogeny
    nodespecies::Matrix{Bool}
end
