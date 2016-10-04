
function getRobject(name::String)
    R"""
    obj = get($name)
    obj$shape <- NULL
    #if(!is.null(obj$shape)){
    #    if(inherits(obj$shape, "RasterLayer"))
    #        obj$shape <- as.matrix(obj$shape)
    #    } # the idea is to later support shape types - not now
    obj$sitestat = obj$coords@data
    obj$coords = sp::coordinates(dat$coords)
    """
    rcopy(R"obj")
end

function Assemblage{T <: Symbol, S}(rdict::Dict{T, S})
    cd_type = rdict[:type] == "grid" ? griddata : pointdata #this code only works as long as there are only those two types
    Assemblage(rdict[:comm], rdict[:coords] )#,
        #cdtype = cd_type, traits = rdict[:species_stats],
        #sitestats = rdict[:sitestat])
end
