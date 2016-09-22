coords(sd::SpatialData) = sd.site.coords
coordtype(as::Union{Assmbl, SiteData}) = as.site.cdtype

addshape!(as::Union{Assmbl, SiteData}, shape::Shapefile.Handle) = (as.site.shape = shape)
deleteshape!(as::Union{Assmbl, SiteData}) = (as.site.shape = Nullable{Shapefile.Handle}())
