import Random

"""
    CommunityRandomizer

The generator returned by `matrixrandomizer` for a Boolean `ComMatrix` or
`Assemblage`. It holds a private copy of the community and a RandomBooleanMatrices
generator over its occurrence matrix; `rand!` installs a new draw in that copy and
returns it, and `rand` returns a copy.
"""
struct CommunityRandomizer{T, G}
    community::T
    generator::G
end

Base.show(io::IO, r::CommunityRandomizer) =
    print(io, "matrixrandomizer for a $(typeof(r.community).name.name) with $(nspecies(r.community)) species and $(nsites(r.community)) sites")

matrixrandomizer(com::ComMatrix, rng = Xoroshiro128Plus();
                method::matrixrandomizations = curveball) = "Only defined for Boolean Assemblages"
matrixrandomizer(asm::SEAssemblage, rng = Xoroshiro128Plus();
                method::matrixrandomizations = curveball) = "Only defined for Boolean Assemblages"

"""
    matrixrandomizer(com [,rng]; method = curveball, trades)

Create a generator of random communities with the same occupancy and richness as
the Boolean `ComMatrix` or `Assemblage` `com`. `method` and `trades` are passed on
to RandomBooleanMatrices' `matrixrandomizer`, which describes the available
`matrixrandomizations`. Draw from the generator with `rand` or `rand!`.
"""
matrixrandomizer(com::ComMatrix{Bool}, rng = Xoroshiro128Plus(); kwargs...) =
    _communityrandomizer(copy(com), rng; kwargs...)
matrixrandomizer(asm::SEAssemblage{Bool}, rng = Xoroshiro128Plus(); kwargs...) =
    _communityrandomizer(copy(asm), rng; kwargs...)

# RandomBooleanMatrices' fast path takes a SparseMatrixCSC{Bool, Int} itself, so the
# generator gets its own copy of the CSC matrix inside the community's OccMatrix.
_communityrandomizer(com, rng; kwargs...) =
    CommunityRandomizer(com, matrixrandomizer(copy(_csc(occurrences(com))), rng; kwargs...))

Random.rand(r::CommunityRandomizer) = copy(Random.rand!(r))
function Random.rand!(r::CommunityRandomizer)
    _setoccurrences!(_commatrix(r.community), Random.rand!(r.generator))
    r.community
end

_commatrix(com::ComMatrix) = com
_commatrix(asm::SEAssemblage) = commatrix(asm)

# Install a new occurrence matrix; the storage drops its transpose, which is rebuilt
# on the next species query. The draw is copied so the community never shares
# storage with the generator's state.
function _setoccurrences!(com::ComMatrix, occ)
    _replacecsc!(occurrences(com), copy(occ))
    com
end
