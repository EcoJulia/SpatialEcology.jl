
mutable struct DispersionField{T} <: SESpatialData{T}
    site::SELocations
    DFmat::Matrix{Int}

    function DispersionField(asm::Assemblage)
        df = asm.occ.occurrences' * asm.occ.occurrences
        new{typeof(asm.site.coords)}(asm.site, df)
    end
end
