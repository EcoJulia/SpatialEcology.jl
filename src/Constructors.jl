
Assemblage(assm::Assmbl) = Assemblage(assm.site, assm.occ) # Not a copy constructor - I think those are automatic? Just a function that will reduce derived types to the base type

# I think all of the constructors should be able to take traits, sites, dropemptys etc.

# a constructor that takes occ and coords as one single DataFrame format and separates them
function Assemblage(occ::DataFrames.DataFrame; kwargs...)
  occ, coords = parsesingleDataFrame(occ)
  Assemblage(occ, coords; kwargs...)
end


# a constructor that takes occ as a DataFrame #should this just be kwargs...?
function Assemblage(occ::DataFrames.DataFrame, coords::Union{AbstractMatrix, DataFrames.DataFrame}; kwargs...)
  occ = parseDataFrame(occ)
  Assemblage(occ, coords; kwargs...)
end

# a constructor that takes occ as a normal matrix
#function Assemblage(occ::AbstractMatrix, coords::Union{AbstractMatrix, DataFrames.DataFrame},
#      sites::Vector{String}, species::Vector{String}; cdtype::coordstype = auto,
#      dropemptyspecies::Bool = true, dropemptysites::Bool = true, match_to_coords = true,
#      shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

function Assemblage(occ::AbstractMatrix, coords::Union{AbstractMatrix, DataFrames.DataFrame}, sites::Vector{String}, species::Vector{String}; kwargs...)
  occ = NamedArrays.NamedArray(occ, (sites, species))
  Assemblage(occ, coords; kwargs...)
end

# a constructor that takes coords as a data.frame
function Assemblage(occ::NamedArrays.NamedArray, coords::DataFrames.DataFrame; kwargs...)

  if DataFrames.ncol(coords) == 2 && all(map(x -> x<:Number, eltypes(coords)))
    coords = dataFrametoNamedMatrix(coords, sparsematrix = false, dimnames = ("sites", "coordinates"))
elseif DataFrames.ncol(coords) == 3
    xind, yind = guess_xycols(coords)
    siteind = setdiff(1:3, [xind, yind])[1]
    coords = dataFrametoNamedMatrix(coords[[xind, yind]], coords[siteind], sparsematrix = false, dimnames = ("sites", "coordinates"))
  else
    error("coords must be a DataFrame with a column for sites and two columns for coordinates")
  end

  Assemblage(occ, coords; kwargs...)
end

Assemblage(occ::NamedArrays.NamedArray, coords::AbstractMatrix; kwargs...) = Assemblage(ComMatrix(occ), coords; kwargs...)


#Assemblage(occ::ComMatrix, sitedata::SiteData; kwargs...) =
#          Assemblage(occ, sitedata.site.coords, cdtype = sitedata.site.cdtype,
#          sitestats = sitedata.site.sitestats, shape = sitedata.site.shape,
#          kwargs...)

Assemblage(occ::ComMatrix, coords::AbstractMatrix; kwargs...) = Assemblage(occ, NamedArrays.NamedArray(coords); kwargs...)

function Assemblage(occ::ComMatrix, coords::NamedArrays.NamedArray;
      dropemptyspecies::Bool = true, dropemptysites::Bool = true, match_to_coords = true,
      traits = DataFrames.DataFrame(name = specnames(occ)), sitestats = DataFrames.DataFrame(sites = sitenames(occ)),
      cdtype::coordstype = auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

    if match_to_coords
        occ, coords, sitestats = match_commat_coords(occ, coords, sitestats)
    end

    Assemblage(createSiteFields(coords, cdtype, sitestats, shape), OccFields(occ, traits))
  end

function Assemblage{T <: Union{Bool, Int}, S <: SiteFields}(site::S, occ::OccFields{T};
    dropemptyspecies::Bool = true, dropemptysites::Bool = true)

    if dropemptyspecies
        dropspecies!(occ)
    end
    if dropemptysites
        dropsites!(occ, site)
    end
    Assemblage{S}{T}(site, occ)
end

function createSiteFields(coords::AbstractMatrix, cdtype::coordstype = auto,  #by design, this is not type stable, but maybe that is OK for type constructors
        sitestats = DataFrames.DataFrame(sites = sitenames(occ)),
        shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

    cdtype == pointdata && return PointData(coords, sitestats, shape)
    cdtype == griddata && return GridData(coords, sitestats, shape)
    if cdtype == auto
        try
            return GridData(coords, sitestats, shape)
        catch
            return PointData(coords, sitestats, shape)
        end
    end
end


OccFields{T <: Union{Bool, Int}}(commatrix::ComMatrix{T}, traits::DataFrames.DataFrame) = OccFields{T}(commatrix, traits)
OccFields(com::ComMatrix) = OccFields(com, DataFrames.DataFrame(id = specnames(commatrix)))

function GridData(coords::NamedArrays.NamedMatrix{Float64},
        sitestats::DataFrames.DataFrame = DataFrames.DataFrame(id = 1:size(coords,1)),
        shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())
    grid = creategrid(coords)
    indices = getindices(coords, grid)
    GridData(indices, grid, sitestats, shape)
end
