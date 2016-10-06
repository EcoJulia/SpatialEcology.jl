
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
    obj$tempshapefile <- FALSE
    if(!inherits(obj, "distrib_data"))
        stop(paste(name, "is not a distrib_data object!"))
    if(!is.null(obj$shape)){
        print(class(obj$shape))
        print(inherits(obj$shape, "SpatialPolygons"))
        if(inherits(obj$shape, "SpatialPolygons")){
            if(requireNamespace("rgdal")){
                obj$tempshapefile <- TRUE
                rgdal::writeOGR(obj$shape, ".", "tempshapefile", driver = "ESRI Shapefile", overwrite_layer = T)
            } else {
                if(requireNamespace("maptools")){
                    maptools::writePolyShape(obj$shape, "tempshapefile")
                    obj$tempshapefile <- TRUE
                } else print("Shapefile dropped as rgdal or maptools were not installed")
            }
        } else warning("shape element dropped as only polygon shapefiles are currently supported")

        obj$shape <- NULL
    }
    obj$sitestat = obj$coords@data
    obj$coords = sp::coordinates(obj$coords)
    """
    rcopy(R"obj")
end

function Assemblage{T <: Symbol, S}(rdict::Dict{T, S})
    if rdict[:tempshapefile]
        shp = open("tempshapefile.shp") do fd
            read(fd, Shapefile.Handle)
        end
        shp = Nullable(shp)   # I really don't see how it works with these silly nullables
    else
        shp = Nullable{Shapefile.Handle}()
    end

    cd_type = rdict[:type] == "grid" ? griddata : pointdata #this code only works as long as there are only those two types
    Assemblage(rdict[:comm], rdict[:coords] ,
        cdtype = cd_type, traits = rdict[:species_stats],
        sitestats = rdict[:sitestat], shape = shp)
end
