
Assemblage(assm::Assmbl) = Assemblage(assm.site, assm.occ) # Not a copy constructor - I think those are automatic? Just a function that will reduce derived types to the base type

# a constructor that takes occ and coords as one single DataFrame format and separates them
function Assemblage(occ::DataFrames.DataFrame; cdtype::coordstype = auto,
      shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

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
          Assemblage(occ, sitedata.SiteFields.coords, cdtype = sitedata.SiteFields.cdtype,
          sitestats = sitedata.SiteFields.sitestats, shape = sitedata.SiteFields.shape,
          dropemptyspecies = dropemptyspecies, dropemptysites = dropemptysites, match_to_coords = match_to_coords)

function Assemblage(occ::ComMatrix, coords::AbstractMatrix; dropemptyspecies::Bool = true,
      dropemptysites::Bool = true, match_to_coords = true,
      traits = DataFrames.DataFrame(name = specnames(occ)), sitestats = DataFrames.DataFrame(sites = sitenames(occ)),
      cdtype::coordstype = auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

    match_to_coords && match_commat_coords!(occ, coords, sitestats)
    dropemptyspecies && dropspecies!(occ, traits)
    dropemptysites && dropsites!(occ, coords, sitestats)

    Assemblage(SiteFields(coords, cdtype, sitestats, shape), OccFields(occ, traits))
  end

Assemblage{T <: Union{Bool, Int}}(site::SiteFields, occ::OccFields{T}) = Assemblage{T}(site, occ)

OccFields{T <: Union{Bool, Int}}(commatrix::ComMatrix{T}, traits::DataFrames.DataFrame) = OccFields{T}(commatrix, traits)
OccFields(com::ComMatrix) = OccFields(com, DataFrames.DataFrame(id = specnames(commatrix)))
