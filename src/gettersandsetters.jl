coordinates(sd::SpatialData) = coordinates(sd.site)
coordinates(pd::PointData) = pd.coords
function coordinates(gd::GridData)
    index1 = xrange(gd.grid)[gd.indices[:,1]]
    index2 = yrange(gd.grid)[gd.indices[:,2]]
    NamedArrays.NamedArray(hcat(index1, index2), allnames(gd.indices), dimnames(gd.indices))
end

addshape!(as::Union{Assmbl, SiteData}, shape::Shapefile.Handle) = (as.site.shape = shape)
deleteshape!(as::Union{Assmbl, SiteData}) = (as.site.shape = Nullable{Shapefile.Handle}())

traits(occ::OccFields) = occ.traits
traits(asm::Assemblage) = traits(asm.occ)
