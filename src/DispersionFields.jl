
mutable struct DispersionField <: SESpatialData
    site::SELocations
    DFmat::Matrix{Int}

    function DispersionField(asm::Assemblage)
        df = asm.occ.occurrences' * asm.occ.occurrences
        new(asm.site, df)
    end
end
