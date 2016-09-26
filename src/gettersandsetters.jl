coordinates(sd::SpatialData) = sd.site.coords

addshape!(as::Union{Assmbl, SiteData}, shape::Shapefile.Handle) = (as.site.shape = shape)
deleteshape!(as::Union{Assmbl, SiteData}) = (as.site.shape = Nullable{Shapefile.Handle}())

traits(occ::OccFields) = occ.traits
traits(asm::Assemblage) = traits(asm.occ)
