
function getRobject(file::String, name::String)
    R"""
    library(nodiv)
    load($file)
    """
    getRobject(name)
end

function getRobject(name::String)
    R"""
    obj = get($name)
    if(!inherits(obj, "distrib_data"))
        stop(paste(name, "is not a distrib_data object!"))

    obj$shape <- NULL
    obj$sitestat = obj$coords@data
    obj$coords = sp::coordinates(obj$coords)
    """
    rcopy(R"obj")
end

function Assemblage{T <: Symbol, S}(rdict::Dict{T, S})

    cd_type = rdict[:type] == "grid" ? griddata : pointdata #this code only works as long as there are only those two types
    Assemblage(rdict[:comm], rdict[:coords] ,
        cdtype = cd_type, traits = rdict[:species_stats],
        sitestats = rdict[:sitestat])
end
