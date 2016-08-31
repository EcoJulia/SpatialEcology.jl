import DataFrames
import Bio.Phylo
import NamedArrays

typealias OccMatrix{T} NamedArrays.NamedArray{T,2,SparseMatrixCSC{T,Int},Tuple{Dict{String,Int},Dict{String,Int}}}  #should not be using inheritance - instead use ownership and overload needed features
typealias AbundanceMatrix OccMatrix{Int} #Maybe have an inner constructor that prevents going under zero
typealias PAMatrix OccMatrix{Bool}

function OccMatrix{T <: Union{Int, Bool}}(x::AbstractMatrix{T})
    println("success")
    NamedArrays.NamedArray(sparse(x))
end

function OccMatrix{T <: Union{Int, Bool}}(x::AbstractMatrix{T}, species::Vector{String}, sites::Vector{String})
    NamedArrays.NamedArray(sparse(x), (sites, species))
end

function tes{T <: Union{Int, Bool}}(x::AbstractMatrix{T})
    println("success")
end

tst = OccMatrix(zeros(Int, 2, 2))

# reconsidering this.
function OccMatrix{T <: Union{Int, Bool}}(x::Matrix{T}, species::Vector{String}, sites::Vector{String})
    NamedArrays.NamedArray(sparse(x), (sites, species))
end

#function PAMatrix(x::Matrix{Bool}; species = :auto, sites = :auto)
#    OccMatrix(x, species = species, sites = sites)
#end
#
#function AbundanceMatrix(x::Matrix{Int}; species = :auto, sites = :auto)
#    OccMatrix(x, species = species, sites = sites)
#end

@enum coordstype grid points

type Assemblage{T}  # A type to keep subtypes together, ensuring that they are all aligned at all times
    coords::Matrix{Float64}  #This should be spatialpoints - not yet implemented? The only one it MUST have
    cdtype::coordstype
    Occurrences::Nullable(OccMatrix{T})
    traits::Nullable(DataFrames.DataFrame)
    sitestats::Nullable(DataFrames.DataFrame)
    phylo::Nullable(Phylo.Phylogeny)
end
