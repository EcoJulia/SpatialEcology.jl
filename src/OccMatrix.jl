const AbstractSparseCSC = SparseArrays.AbstractSparseMatrixCSC

"""
    OccMatrix

The storage behind a `ComMatrix`: its species × sites occurrence matrix, held as a
`SparseMatrixCSC`, together with the transpose (sites × species). The matrix is
compressed by column, so reading a site is cheap but reading a species means
scanning every column; in the transpose a species is a column, so per-species
queries (`occupied(com, sp)`, `noccupied`, `getspecies`, Bool `occupancy`) read
it instead. The transpose is built when first asked for and kept until the
occurrences change.

`OccMatrix` is an `AbstractSparseMatrixCSC`, so SparseArrays' sparse methods work
on it directly, and `occurrences(com)` returns it.

Writing to it with `setindex!` (also through a `view` or a `ComMatrix`),
broadcasting assignment (`occ .= …`), `copyto!`, `fill!` or `dropzeros!` drops the
transpose and removes any stored zeros, so the transpose is never stale and a
written zero is never counted as a presence. Code that changes the stored values
in place by other means, e.g. through `nonzeros(occ)` or `rmul!`, must call
`SpatialEcology.invalidate!(occ)` afterwards.
"""
mutable struct OccMatrix{Tv} <: AbstractSparseCSC{Tv, Int}
    csc::SparseMatrixCSC{Tv, Int}
    t::Union{Nothing, SparseMatrixCSC{Tv, Int}}
    # Takes ownership of `csc` (no copy), removing stored zeros from it in place.
    OccMatrix{Tv}(csc::SparseMatrixCSC{Tv, Int}) where Tv = new{Tv}(dropzeros!(csc), nothing)
end
OccMatrix(csc::SparseMatrixCSC{Tv, Int}) where Tv = OccMatrix{Tv}(csc)
OccMatrix(m::AbstractMatrix{Tv}) where Tv = OccMatrix{Tv}(SparseMatrixCSC{Tv, Int}(sparse(m)))

# AbstractSparseMatrixCSC interface
Base.size(W::OccMatrix) = size(W.csc)
SparseArrays.getcolptr(W::OccMatrix) = SparseArrays.getcolptr(W.csc)
SparseArrays.rowvals(W::OccMatrix) = rowvals(W.csc)
SparseArrays.nonzeros(W::OccMatrix) = nonzeros(W.csc)
# The arrays are W.csc's, so broadcasting `W .= f.(W)` must see that they alias.
Base.dataids(W::OccMatrix) = Base.dataids(W.csc)
Base.copy(W::OccMatrix) = OccMatrix(copy(W.csc))
SparseArrays.sparse(W::OccMatrix) = copy(W.csc)

_csc(W::OccMatrix) = W.csc

"""
    transposed(occ::OccMatrix)

The sites × species transpose of `occ`, built on first use.
"""
transposed(W::OccMatrix) = (t = W.t; t === nothing ? (W.t = permutedims(W.csc)) : t)

"""
    invalidate!(occ::OccMatrix)

Drop the cached transpose, after the stored values were changed in place.
"""
invalidate!(W::OccMatrix) = (W.t = nothing; W)

# Install a new occurrence matrix (same size or not), taking ownership of it.
function _replacecsc!(W::OccMatrix{Tv}, csc::SparseMatrixCSC) where Tv
    W.csc = dropzeros!(csc isa SparseMatrixCSC{Tv, Int} ? csc : SparseMatrixCSC{Tv, Int}(csc))
    W.t = nothing
    W
end

# After any write: no stored zeros, no stale transpose.
_written!(W::OccMatrix) = (dropzeros!(W.csc); W.t = nothing; W)

# Each write method mirrors one of SparseArrays' methods for AbstractSparseMatrixCSC
# (with `OccMatrix` in the first position), so it is the more specific one and no
# call is ambiguous. The work is done on the `SparseMatrixCSC` itself.
const _CSCIndex = Union{Colon, Integer, AbstractVector{<:Integer}}

function Base.setindex!(W::OccMatrix, v, i::Integer, j::Integer)
    # Setting a stored entry to zero leaves it stored; remove it (writing a zero
    # where nothing is stored stores nothing, so needs no clean-up).
    storedzero = iszero(v) && Base.isstored(W.csc, i, j)
    W.csc[i, j] = v
    storedzero && dropzeros!(W.csc)
    W.t = nothing
    W
end
Base.setindex!(W::OccMatrix, v::AbstractVecOrMat, i::Integer, j::Integer) =
    (setindex!(W.csc, v, i, j); _written!(W))
Base.setindex!(W::OccMatrix, v, I::_CSCIndex, J::_CSCIndex) =
    (setindex!(W.csc, v, I, J); _written!(W))
Base.setindex!(W::OccMatrix, v::AbstractVecOrMat, I::_CSCIndex, J::_CSCIndex) =
    (setindex!(W.csc, v, I, J); _written!(W))
Base.setindex!(W::OccMatrix, v::AbstractArray, I::AbstractMatrix{Bool}) =
    (setindex!(W.csc, v, I); _written!(W))
Base.setindex!(W::OccMatrix, v::AbstractArray, I::AbstractVector{<:Integer}) =
    (setindex!(W.csc, v, I); _written!(W))

Base.copyto!(W::OccMatrix, src::AbstractSparseCSC) = (copyto!(W.csc, src); _written!(W))
Base.copyto!(W::OccMatrix, src::SparseArrays.AbstractCompressedVector) = (copyto!(W.csc, src); _written!(W))
Base.copyto!(W::OccMatrix, src::Union{Array, BitArray}) = (copyto!(W.csc, src); _written!(W))
Base.copyto!(W::OccMatrix, bc::Broadcast.Broadcasted) = (copyto!(W.csc, bc); _written!(W))
Base.copyto!(W::OccMatrix, bc::Broadcast.Broadcasted{Nothing}) = (copyto!(W.csc, bc); _written!(W))
Base.copyto!(W::OccMatrix, bc::Broadcast.Broadcasted{<:Broadcast.AbstractArrayStyle{0}}) =
    (copyto!(W.csc, bc); _written!(W))
Base.copyto!(W::OccMatrix, bc::Broadcast.Broadcasted{Nothing, <:Any, <:Any, Tuple{}}) =
    (copyto!(W.csc, bc); _written!(W))
Base.copyto!(W::OccMatrix, bc::Broadcast.Broadcasted{<:Union{SparseArrays.HigherOrderFns.SparseMatStyle, SparseArrays.HigherOrderFns.SparseVecStyle}}) =
    (copyto!(W.csc, bc); _written!(W))
Base.copyto!(W::OccMatrix, bc::Broadcast.Broadcasted{SparseArrays.HigherOrderFns.PromoteToSparse}) =
    (copyto!(W.csc, bc); _written!(W))
Base.fill!(W::OccMatrix, x) = (fill!(W.csc, x); _written!(W))
SparseArrays.dropzeros!(W::OccMatrix) = _written!(W)

# SparseArrays fills a view of a sparse matrix by writing into the parent's arrays
# directly, not through setindex! (and `v .= x` on a view ends up here), so a view
# of an OccMatrix needs the same treatment.
function Base.fill!(V::SubArray{Tv, <:Any, <:OccMatrix{Tv}, <:Tuple{Union{Integer, AbstractVector{<:Integer}}, Union{Integer, AbstractVector{<:Integer}}}}, x) where Tv
    W = parent(V)
    fill!(view(W.csc, parentindices(V)...), x)
    _written!(W)
    V
end
