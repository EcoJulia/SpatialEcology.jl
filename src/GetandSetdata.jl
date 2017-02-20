coordinates(sd::SpatialData) = coordinates(sd.site)
coordinates(pd::PointData) = pd.coords

function coordinates(gd::GridData)
    index1 = xrange(gd.grid)[gd.indices[:,1]]
    index2 = yrange(gd.grid)[gd.indices[:,2]]
    NamedArrays.NamedArray(hcat(index1, index2), NamedArrays.names(gd.indices), dimnames(gd.indices))
end

traits(occ::OccFields) = occ.traits
traits(asm::Assemblage) = traits(asm.occ)
