import Random

struct AssemblageRandomizer{P <: Locations}
    asm::Assemblage{Bool, P}
end

assemblage_randomizer(asm::SEAssemblage) = "Only defined for Boolean Assemblages"
function assemblage_randomizer(asm::SEAssemblage{Bool, P}) where P <: SELocations
    ret = AssemblageRandomizer{P}(copy(asm))
    dropzeros!(ret.asm.occ.commatrix.occurrences)
    ret
end

Random.rand(r::AssemblageRandomizer) = copy(rand!(r))
function Random.rand!(r::AssemblageRandomizer)
    RandomBooleanMatrices._curveball!(r.asm.occ.commatrix.occurrences)
    r.asm
end
