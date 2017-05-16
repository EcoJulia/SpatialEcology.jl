
# Big TODO we are still missing the functionality that does the aligning in the constructors

Assemblage(assm::Assmbl) = Assemblage(assm.site, assm.occ) # Not a copy constructor - Just a function that will reduce derived types to the base type

# a constructor that takes occ and coords as one single DataFrame format and separates them
function Assemblage(occ::DataFrames.DataFrame; kwargs...)
  occ, coords = parsesingleDataFrame(occ)
  Assemblage(occ, coords; kwargs...)
end

# a constructor that takes occ as a DataFrame
Assemblage(occ::DataFrames.DataFrame, coords::Union{AbstractMatrix, DataFrames.DataFrame}; kwargs...) = Assemblage(ComMatrix(occ), coords; kwargs...)

Assemblage(occ::AbstractMatrix, coords::Union{AbstractMatrix, DataFrames.DataFrame}, sites::Vector{String}, species::Vector{String}; kwargs...) = Assemblage(ComMatrix(occ, species, sites), coords; kwargs...)

# a constructor that takes coords as a data.frame
function Assemblage(occ::ComMatrix, coords::DataFrames.DataFrame; kwargs...)

    if DataFrames.ncol(coords) == 2 && all(map(x -> x<:Number, eltypes(coords)))
        coords = convert(Array, coords)
    elseif DataFrames.ncol(coords) == 3
        xind, yind = guess_xycols(coords)
        siteind = setdiff(1:3, [xind, yind])[1]
        coords = convert(Array, coords[[xind, yind]])
    else
        error("coords must be a DataFrame with a column for sites and two columns for coordinates")
    end

    Assemblage(occ, coords; kwargs...)
end

function Assemblage(occ::ComMatrix, coords::AbstractMatrix;
      dropemptyspecies::Bool = false, dropemptysites::Bool = false, match_to_coords = true,
      traits = DataFrames.DataFrame(name = specnames(occ)), sitestats = DataFrames.DataFrame(sites = sitenames(occ)),
      cdtype::coordstype = auto)

    if match_to_coords
        occ, coords, sitestats = match_commat_coords(occ, coords, sitestats)
    end

    Assemblage(createSiteFields(coords, cdtype, sitestats), OccFields(occ, traits))
  end

function Assemblage{T <: Union{Bool, Int}, S <: SiteFields}(site::S, occ::OccFields{T};
    dropemptyspecies::Bool = false, dropemptysites::Bool = false)

    if dropemptyspecies
        dropspecies!(occ)
    end
    if dropemptysites
        dropsites!(occ, site)
    end
    Assemblage{S}{T}(site, occ)
end

function createSiteFields(coords::AbstractMatrix, cdtype::coordstype = auto,  #by design, this is not type stable, but maybe that is OK for type constructors
        sitestats = DataFrames.DataFrame(sites = sitenames(occ)))

    cdtype == pointdata && return PointData(coords, sitestats)
    cdtype == griddata && return GridData(coords, sitestats)
    if cdtype == auto
        try
            return GridData(coords, sitestats)
        catch
            return PointData(coords, sitestats)
        end
    end
end


OccFields{T <: Union{Bool, Int}}(commatrix::ComMatrix{T}, traits::DataFrames.DataFrame) = OccFields{T}(commatrix, traits)
OccFields(com::ComMatrix) = OccFields(com, DataFrames.DataFrame(id = specnames(commatrix)))

function GridData(coords::Matrix{Float64},
        sitestats::DataFrames.DataFrame = DataFrames.DataFrame(id = 1:size(coords,1)))
    grid = creategrid(coords)
    indices = getindices(coords, grid)
    GridData(indices, grid, sitestats)
end

function ComMatrix(occ::DataFrames.DataFrame)
    if DataFrames.ncol(occ) == 3 && eltypes(occ)[3] <: String
        println("Data format identified as Phylocom")
        sites = unique(a[1])
        species = unique(a[3])
        is = indexin(a[1], sites)
        js = indexin(a[3], species
        occ = maximum(a[2]) == 1 ? sparse(is, js, true) : sparse(is, js, a[2])
        return ComMatrix(occ, string.(species), string.(sites))
    end

    sites = collect(string.(occ[1]))
    species = collect(string.(names(occ)))[2:end]

    try
        occ = dataFrametoSparseMatrix(occ[2:end], Bool)
        println("Matrix data assumed to be presence-absence")
    catch
        occ = dataFrametoSparseMatrix(occ[2:end], Int) #TODO This line means that this code is not completely type stable.
        println("Matrix data assumed to be abundances, minimum $(minimum(occ)), maximum $(maximum(occ))")
    end

    ComMatrix(occ, species, sites)
end
