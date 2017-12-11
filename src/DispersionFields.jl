
mutable struct DispersionField <: SpatialData
    site::SiteFields
    DFmat::Matrix{Int}

    function DispersionField(asm::Assemblage)
        df = asm.occ.occurrences * asm.occ.occurrences'
        new(asm.site, df)
    end
end
