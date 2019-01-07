import Random


matrixrandomizer(asm::SEAssemblage, rng = Xoroshiro128Plus();
                method::matrixrandomizations = curveball) = "Only defined for Boolean Assemblages"

function matrixrandomizer(asm::S, rng = Xoroshiro128Plus();
                        method::matrixrandomizations = curveball) where S <: SEAssemblage{Bool} where P
    as = copy(asm)
    ret = MatrixGenerator{typeof(rng), typeof(as)}(as, method, rng)
    dropzeros!(ret.m.occ.commatrix.occurrences)
    ret
end

Random.rand(r::MatrixGenerator{R, Assemblage}) where R = copy(rand!(r))
function Random.rand!(r::MatrixGenerator{R, A}) where {R} where {A <: Assemblage}
    RandomBooleanMatrices._curveball!(r.m.occ.commatrix.occurrences)
    r.m
end
