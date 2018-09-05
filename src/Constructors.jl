
# Big TODO we are still missing the functionality that does the aligning in the constructors

Assemblage(assm::SEAssemblage) = Assemblage(assm.site, assm.occ) # Not a copy constructor - Just a function that will reduce derived types to the base type


# a constructor that takes occ and coords as one single DataFrame format and separates them
function Assemblage(occ::DataFrames.DataFrame; kwargs...)
  occ, coords = parsesingleDataFrame(occ)
  Assemblage(occ, coords; kwargs...)
end

# a constructor that takes occ as a DataFrame
Assemblage(occ::DataFrames.DataFrame, coords::Union{AbstractMatrix, DataFrames.DataFrame}; sitecolumns = true, kwargs...) = Assemblage(ComMatrix(occ; sitecolumns = sitecolumns), coords; kwargs...)

Assemblage(occ::AbstractMatrix, coords::Union{AbstractMatrix, DataFrames.DataFrame}, sites::Vector{String}, species::Vector{String}; sitecolumns = true, kwargs...) = Assemblage(ComMatrix(occ, species, sites; sitecolumns = sitecolumns), coords; kwargs...)

# a constructor that takes coords as a data.frame
function Assemblage(occ::ComMatrix, coords::DataFrames.DataFrame; kwargs...)

    if DataFrames.ncol(coords) == 2 && all(map(x -> x<:Number, eltypest(coords)))
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
      traits = DataFrames.DataFrame(name = specnames(occ)), sitestat = DataFrames.DataFrame(sites = sitenames(occ)),
      cdtype::coordstype = auto)

    if match_to_coords
        occ, coords, sitestat = match_commat_coords(occ, coords, sitestat)
    end

    Assemblage(createSiteFields(coords, cdtype, sitestat), OccFields(occ, traits))
  end

function Assemblage(site::S, occ::OccFields{T};
    dropemptyspecies::Bool = false, dropemptysites::Bool = false) where {T <: OccTypes, S <: SiteFields}

    if dropemptyspecies
        dropspecies!(occ)
    end
    if dropemptysites
        dropsites!(occ, site)
    end
    Assemblage{S}{T}(site, occ)
end

function createSiteFields(coords::AbstractMatrix, cdtype::coordstype = auto,  #by design, this is not type stable, but maybe that is OK for type constructors
        sitestat = DataFrames.DataFrame(sites = sitenames(occ)))

    cdtype == pointdata && return PointData(coords, sitestat)
    cdtype == griddata && return GridData(coords, sitestat)
    if cdtype == auto
        try
            return GridData(coords, sitestat)
        catch
            return PointData(coords, sitestat)
        end
    end
end


OccFields(commatrix::ComMatrix{T}, traits::DataFrames.DataFrame) where T <: OccTypes = OccFields{T}(commatrix, traits)
OccFields(com::ComMatrix) = OccFields(com, DataFrames.DataFrame(id = specnames(commatrix)))


function ComMatrix(occ::DataFrames.DataFrame; sitecolumns = true)
    if DataFrames.ncol(occ) == 3 && eltypest(occ)[3] <: AbstractString
        println("Data format identified as Phylocom")
        sites = unique(occ[1])
        species = unique(occ[3])
        js = indexin(occ[1], sites)
        is = indexin(occ[3], species)
        occ = maximum(occ[2]) == 1 ? sparse(is, js, true) : sparse(is, js, occ[2])
        return ComMatrix(occ, string.(collect(species)), string.(collect(sites)))
    end

    if eltypet(occ[1]) <: AbstractString
        species = string.(collect(occ[1]))
        occ = occ[2:end]
        sites = string.(collect(names(occ)))
    else
        species = string.(1:DataFrames.nrow(occ))
        sites = string.(collect(names(occ)))
    end

    try
        occ = dataFrametoSparseMatrix(occ, Bool)
        println("Matrix data assumed to be presence-absence")
    catch
        try
            occ = dataFrametoSparseMatrix(occ, Int) #TODO These lines mean that this code is not completely type stable.
            println("Matrix data assumed to be abundances, minimum $(minimum(occ)), maximum $(maximum(occ))")
        catch
            occ = dataFrametoSparseMatrix(occ, Float64)
            println("Matrix data assumed to be relative abundances, minimum $(minimum(occ)), maximum $(maximum(occ))")
            (minimum(occ) < 0 || maximum(occ) > 1) && info("Values don't fall in the 0,1 range")
        end
    end

    if !sitecolumns
        return ComMatrix(collect(occ'), sites, species)
    end

    ComMatrix(occ, species, sites)
end

ComMatrix(occurrences::AbstractMatrix, specnames, sitenames; sitecolumns = true) =
    ComMatrix(occurrences; specnames = specnames, sitenames = sitenames, sitecolumns = sitecolumns)


function ComMatrix(occurrences; specnames = :auto, sitenames = :auto, sitecolumns = true)
    if !sitecolumns
        occurrences = occurrences'
    end
    if sitenames == :auto
        sitenames = ["site$i" for i in 1:size(occurrences, 2)]
    end
    if specnames == :auto
        specnames = ["species$i" for i in 1:size(occurrences, 1)]
    end

    length(specnames) == size(occurrences, 1) || throw(ArgumentError("length of specnames ($(length(specnames))) different from number of species in occurrences ($(size(occurrences, 1)))"))
    length(sitenames) == size(occurrences, 2) || throw(ArgumentError("length of sitenames ($(length(sitenames))) different from number of sites in occurrences ($(size(occurrences, 2)))"))
    ComMatrix{eltype(occurrences)}(sparse(occurrences), string.(specnames), string.(sitenames))
end

function GridData(coords::AbstractMatrix{<:Union{AbstractFloat, Missings.Missing}},
        sitestats::DataFrames.DataFrame = DataFrames.DataFrame(id = 1:size(coords,1)))
    grid = creategrid(coords)
    indices = getindices(coords, grid)
    GridData(indices, grid, sitestats)
end
