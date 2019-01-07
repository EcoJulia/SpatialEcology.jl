import Random

matrixrandomizer(com::ComMatrix, rng = Xoroshiro128Plus();
                method::matrixrandomizations = curveball) = "Only defined for Boolean Assemblages"
matrixrandomizer(asm::SEAssemblage, rng = Xoroshiro128Plus();
                method::matrixrandomizations = curveball) = "Only defined for Boolean Assemblages"

function matrixrandomizer(asm::S, rng = Xoroshiro128Plus();
                        method::matrixrandomizations = curveball) where S <: SEAssemblage{Bool} where P
    as = copy(asm)
    ret = MatrixGenerator{typeof(rng), typeof(as)}(as, method, rng)
    dropzeros!(ret.m.occ.commatrix.occurrences)
    ret
end

function matrixrandomizer(com::C, rng = Xoroshiro128Plus();
                        method::matrixrandomizations = curveball) where C <: ComMatrix{Bool}
    cm = copy(com)
    ret = MatrixGenerator{typeof(rng), typeof(cm)}(as, method, rng)
    dropzeros!(ret.m.occurrences)
    ret
end

Random.rand(r::MatrixGenerator{R, <:Assemblage}) where R = copy(rand!(r))
function Random.rand!(r::MatrixGenerator{R, A}) where {R} where {A <: Assemblage}
    RandomBooleanMatrices._curveball!(r.m.occ.commatrix.occurrences)
    r.m
end

Random.rand(r::MatrixGenerator{<:ComMatrix}) where R = copy(rand!(r))
function Random.rand!(r::MatrixGenerator{C}) where {C <: ComMatrix}
    RandomBooleanMatrices._curveball!(r.m.occurrences)
    r.m
end
