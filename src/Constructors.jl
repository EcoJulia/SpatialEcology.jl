
include("Constructor_helperfunctions.jl")

function PhyloAssemblage(site::SiteFields, occ::OccFields, phylo::Phylogeny)
  nodespec = createNodeBySpeciesMatrix(phylo)
  PhyloAssemblage(site, occ, PhyloFields(phylo, nodespec))
end

PhyloAssemblage(assm::Assemblage, phylo::Phylogeny) = PhyloAssemblage(assm.site, assm.occ, phylo)

Assemblage(assm::Assmbl) = Assemblage(assm.site, assm.occ) # Not a copy constructor - I think those are automatic?

function Assemblage(occ::DataFrames.DataFrame; cdtype::coordstype = :auto,
      shape::Nullable{ShapeFiles.ShapeFile} = Nullable{ShapeFiles.ShapeFile}())

  occ, coords = parsesingleDataFrame(occ)
  Assemblage(occ, coords, cdtype, shape)
end

function Assemblage(occ::DataFrames.DataFrame, coords::Union{AbstractMatrix, DataFrames.DataFrame};
      cdtype::coordstype = :auto, shape::Nullable{ShapeFiles.ShapeFile} = Nullable{ShapeFiles.ShapeFile}())

  occ = parseDataFrame(occ)
  Assemblage(occ, coords, cdtype, shape)
end

function Assemblage(occ::AbstractMatrix, coords::Union{AbstractMatrix, DataFrames.DataFrame},
      sites::Vector{String}, species::Vector{String}; cdtype::coordstype = :auto,
      shape::Nullable{ShapeFiles.ShapeFile} = Nullable{ShapeFiles.ShapeFile}())

  occ = NamedArrays.NamedArray(occ, (sites, species))
  Assemblage(occ, coords, cdtype, shape)
end

function Assemblage(occ::NamedArrays.NamedArray, coords::DataFrames.DataFrame:
      cdtype::coordstype = :auto, shape::Nullable{ShapeFiles.ShapeFile} = Nullable{ShapeFiles.ShapeFile}())

  ## Here is a constructor that does work on coords
end

function Assemblage(occ::NamedArrays.NamedArray, coords::AbstractMatrix:
      cdtype::coordstype = :auto, shape::Nullable{ShapeFiles.ShapeFile} = Nullable{ShapeFiles.ShapeFile}())

  ## Here is the constructor that actually constructs the object and does the work from nodiv

end
