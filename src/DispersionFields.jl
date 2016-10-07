
type DispersionField
    site::SiteFields
    DFmat::NamedMatrix{Int}

    function DispersionField(asm::Assemblage)
        df = asm.occ.occurrences * asm.occ.occurrences'
        new(asm.site, df)
    end
end
