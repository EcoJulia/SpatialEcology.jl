coordinates(sd::SpatialData) = coordinates(sd.site)
coordinates(pd::AbstractPointData) = pd.coords

function coordinates(gd::AbstractGridData)
    index1 = xrange(gd.grid)[gd.indices[:,1]]
    index2 = yrange(gd.grid)[gd.indices[:,2]]
    NamedArrays.NamedArray(hcat(index1, index2), NamedArrays.names(gd.indices), dimnames(gd.indices))
end

traits(occ::AbstractOccFields) = occ.traits
traits(asm::Assmbl) = traits(asm.occ)

sitestats(asm::Assmbl) = asm.site.sitestats

traitnames(asm::Assmbl) = names(traits(asm))
sitestatnames(asm::Assmbl) = names(sitestats(asm))
