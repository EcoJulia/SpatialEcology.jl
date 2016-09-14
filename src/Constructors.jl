
include("Constructor_helperfunctions.jl")

function PhyloAssemblage(site::SiteFields, occ::OccFields, phylo::Phylogeny)
  nodespec = createNodeBySpeciesMatrix(phylo)
  PhyloAssemblage(site, occ, PhyloFields(phylo, nodespec))
end

PhyloAssemblage(assm::Assemblage, phylo::Phylogeny) = PhyloAssemblage(assm.site, assm.occ, phylo)

Assemblage(assm::Assmbl) = Assemblage(assm.site, assm.occ) # Not a copy constructor - I think those are automatic? Just a function that will reduce derived types to the base type

# a constructor that takes occ and coords as one single DataFrame format and separates them
function Assemblage(occ::DataFrames.DataFrame; cdtype::coordstype = :auto,
      shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  occ, coords = parsesingleDataFrame(occ)
  Assemblage(occ, coords, cdtype, shape)
end

# a constructor that takes occ as a DataFrame
function Assemblage(occ::DataFrames.DataFrame, coords::Union{AbstractMatrix, DataFrames.DataFrame};
      dropemptyspecies::Bool = true, dropemptysites::Bool = true, match_to_coords = true,
      cdtype::coordstype = :auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  occ = parseDataFrame(occ)
  Assemblage(occ, coords, cdtype, shape)
end

# a constructor that takes occ as a normal matrix
function Assemblage(occ::AbstractMatrix, coords::Union{AbstractMatrix, DataFrames.DataFrame},
      sites::Vector{String}, species::Vector{String}; cdtype::coordstype = :auto,
      dropemptyspecies::Bool = true, dropemptysites::Bool = true, match_to_coords = true,
      shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  occ = NamedArrays.NamedArray(occ, (sites, species))
  Assemblage(occ, coords, cdtype, shape)
end

# a constructor that takes coords as a data.frame
function Assemblage(occ::NamedArrays.NamedArray, coords::DataFrames.DataFrame; dropemptyspecies::Bool = true,
      dropemptysites::Bool = true, match_to_coords = true,
      cdtype::coordstype = :auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  if ncol(coords) == 2 && all(map(x -> x<:Number, eltypes(dd)))
    coords = dataFrametoNamedMatrix(coords)
  elseif ncol(coords) == 3 && all(map(x -> x<:Number, eltypes(dd)[2:3]))
    coords = dataFrametoNamedMatrix(coords[2:3], coords[1])
  else
    error("coords must be a DataFrame with a column for sites and two columns for coordinates")
  end

  Assemblage(occ, coords, cdtype, shape)
end

function Assemblage(occ::NamedArrays.NamedArray, coords::AbstractMatrix;
      sitestats::DataFrames.DataFrame = DataFrames.DataFrame,
      dropemptyspecies::Bool = true, dropemptysites::Bool = true, match_to_coords = true,
      cdtype::coordstype = :auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

  Assemblage(ComMatrix(occ), coords, cdtype, shape)
end

Assemblage(occ::ComMatrix, sitedata::SiteData; dropemptyspecies::Bool = true,
      dropemptysites::Bool = true, match_to_coords = true) =
          Assemblage(occ, sitedata.SiteFields.coords, cdtype = sitedata.SiteFields.cdtype,
          sitestats = sitedata.SiteFields.sitestats, shape = sitedata.SiteFields.shape,
          dropemptyspecies = dropemptyspecies, dropemptysites = dropemptysites, match_to_coords = match_to_coords)

function Assemblage(occ::ComMatrix, coords::AbstractMatrix; dropemptyspecies::Bool = true,
      dropemptysites::Bool = true, match_to_coords = true,
      cdtype::coordstype = :auto, shape::Nullable{Shapefile.Handle} = Nullable{Shapefile.Handle}())

    match_to_coords && match_commat_coords!(occ, coords, sitestats)
    dropemptyspecies && dropspecies!(occ, coords, sitestats)
    dropemptysites && dropsites!(occ, traits)

    Assemblage(SiteFields(coords, cdtype, sitestats, shape), OccFields(occ, traits))
  end
