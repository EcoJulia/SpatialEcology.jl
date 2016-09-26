
Assemblage(assm::Assmbl) = Assemblage(assm.site, assm.occ) # Not a copy constructor - I think those are automatic? Just a function that will reduce derived types to the base type

# I think all of the constructors should be able to take traits, sites, dropemptys etc.

# a constructor that takes occ and coords as one single DataFrame format and separates them
function Assemblage(occ::DataFrames.DataFrame; dropemptyspecies::Bool = true,
            dropemptysites::Bool = true, match_to_coords = true,
            traits = DataFrames.DataFrame(name = specnames(occ)), sitestats = DataFrames.DataFrame(sites = sitenames(occ)),
            cdtype::coordstype = auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  occ, coords = parsesingleDataFrame(occ)
  Assemblage(occ, coords, cdtype = cdtype, shape = shape)
end

# a constructor that takes occ as a DataFrame
function Assemblage(occ::DataFrames.DataFrame, coords::Union{AbstractMatrix, DataFrames.DataFrame};
      dropemptyspecies::Bool = true, dropemptysites::Bool = true, match_to_coords = true,
      cdtype::coordstype = auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  occ = parseDataFrame(occ)
  Assemblage(occ, coords, cdtype = cdtype, shape = shape)
end

# a constructor that takes occ as a normal matrix
function Assemblage(occ::AbstractMatrix, coords::Union{AbstractMatrix, DataFrames.DataFrame},
      sites::Vector{String}, species::Vector{String}; cdtype::coordstype = auto,
      dropemptyspecies::Bool = true, dropemptysites::Bool = true, match_to_coords = true,
      shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  occ = NamedArrays.NamedArray(occ, (sites, species))
  Assemblage(occ, coords, cdtype = cdtype, shape = shape)
end

# a constructor that takes coords as a data.frame
function Assemblage(occ::NamedArrays.NamedArray, coords::DataFrames.DataFrame; dropemptyspecies::Bool = true,
      dropemptysites::Bool = true, match_to_coords = true,
      cdtype::coordstype = auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  if ncol(coords) == 2 && all(map(x -> x<:Number, eltypes(coords)))
    coords = dataFrametoNamedMatrix(coords, sparsematrix = false)
  elseif ncol(coords) == 3
    xind, yind = guess_xycols(coords)
    siteind = setdiff(1:3, [xind, yind])[1]
    coords = dataFrametoNamedMatrix(coords[[xind, yind]], coords[siteind], sparsematrix = false)
  else
    error("coords must be a DataFrame with a column for sites and two columns for coordinates")
  end

  Assemblage(occ, coords, cdtype = cdtype, shape = shape)
end

function Assemblage(occ::NamedArrays.NamedArray, coords::AbstractMatrix;
      dropemptyspecies::Bool = true, dropemptysites::Bool = true, match_to_coords = true,
      cdtype::coordstype = auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  Assemblage(ComMatrix(occ), coords, cdtype = cdtype, shape = shape)
end

Assemblage(occ::ComMatrix, sitedata::SiteData; dropemptyspecies::Bool = true,
      dropemptysites::Bool = true, match_to_coords = true) =
          Assemblage(occ, sitedata.site.coords, cdtype = sitedata.site.cdtype,
          sitestats = sitedata.site.sitestats, shape = sitedata.site.shape,
          dropemptyspecies = dropemptyspecies, dropemptysites = dropemptysites, match_to_coords = match_to_coords)

function Assemblage(occ::ComMatrix, coords::AbstractMatrix;
      dropemptyspecies::Bool = true, dropemptysites::Bool = true, match_to_coords = true,
      traits = DataFrames.DataFrame(name = specnames(occ)), sitestats = DataFrames.DataFrame(sites = sitenames(occ)),
      cdtype::coordstype = auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

    if match_to_coords
        occ, coords, sitestats = match_commat_coords!(occ, coords, sitestats)
    end

    if dropemptyspecies
        occ, traits = dropspecies!(occ, traits)
    end

    if dropemptysites
        occ, coords, sitestats = dropsites!(occ, coords, sitestats)
    end

    Assemblage(createSiteFields(coords, cdtype, sitestats, shape), OccFields(occ, traits))
  end

function Assemblage{T <: Union{Bool, Int}}(site::SiteFields, occ::OccFields{T};
    dropemptyspecies::Bool = true, dropemptysites::Bool = true)

    if dropemptyspecies
        occ.commatrix, occ.traits = dropspecies!(occ.commatrix, occ.traits)
    end

    if dropemptysites
        occ.commatrix, site.coords, site.sitestats = dropsites!(occ.commatrix, site.coords, site.sitestats)
    end

    Assemblage{T}(site, occ)
end

function createSiteFields(coords::AbstractMatrix, cdtype::coordstype = auto,  #by design, this is not type stable, but maybe that is OK for type constructors
        sitestats = DataFrames.DataFrame(sites = sitenames(occ)),
        shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

    cdtype == points && return PointData(coords, sitestats, shape)
    cdtype == grid && return GridData(coords, sitestats, shape)
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
    grid, coords = creategrid(coords)
    GridData(coords, grid, sitestats, shape)
end
