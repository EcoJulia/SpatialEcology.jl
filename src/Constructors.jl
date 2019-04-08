
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
        coords = convert(Matrix, coords)
    elseif DataFrames.ncol(coords) == 3
        xind, yind = guess_xycols(coords)
        siteind = setdiff(1:3, [xind, yind])[1]
        coords = convert(Matrix, coords[[xind, yind]])
    else
        error("coords must be a DataFrame with a column for sites and two columns for coordinates")
    end

    Assemblage(occ, coords; kwargs...)
end

function Assemblage(occ::ComMatrix, coords::AbstractMatrix;
      dropemptyspecies::Bool = false, dropemptysites::Bool = false, match_to_coords = true,
      traits = DataFrames.DataFrame(name = speciesnames(occ)), sitestat = DataFrames.DataFrame(sites = sitenames(occ)),
      cdtype::coordstype = auto)

    if match_to_coords
        occ, coords, sitestat = match_commat_coords(occ, coords, sitestat)
    end

    Assemblage(createLocations(coords, cdtype, sitestat), SpeciesData(occ, traits))
  end

function Assemblage(site::P, occ::SpeciesData{D};
    dropemptyspecies::Bool = false, dropemptysites::Bool = false) where {D <: Real, P <: SELocations}

    if dropemptyspecies
        dropspecies!(occ)
    end
    if dropemptysites
        dropsites!(occ, site)
    end
    Assemblage{D, P}(site, occ)
end

function createLocations(coords::AbstractMatrix, cdtype::coordstype = auto,  #by design, this is not type stable, but maybe that is OK for type constructors
        sitestat = DataFrames.DataFrame(sites = 1:size(coords, 1)))

    cdtype == pointdata && return Locations{PointData}(PointData(coords), sitestat)
    cdtype == griddata && return Locations{GridData}(GridData(coords), sitestat)
    if cdtype == auto
        try
            return Locations{GridData}(GridData(coords), sitestat)
        catch
            return Locations{PointData}(PointData(coords), sitestat)
        end
    end
end

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
            (minimum(occ) < 0 || maximum(occ) > 1) && @info("Values don't fall in the 0,1 range")
        end
    end

    if !sitecolumns
        return ComMatrix(collect(occ'), sites, species)
    end

    ComMatrix(occ, species, sites)
end

function ComMatrix(; species::AbstractVector = error("keyword `species` must be specified"), sites::AbstractVector = error("keyword `sites` must be specified"), abundances = 0)
    usites = unique(sites)
    uspecies = unique(species)
    js = indexin(sites, usites)
    is = indexin(species, uspecies)
    occ = abundances isa AbstractVector && length(unique(abundances)) > 1 ? sparse(is, js, abundances) : sparse(is, js, true)
    return ComMatrix(occ, string.(collect(uspecies)), string.(collect(usites)))
end

ComMatrix(occurrences::AbstractMatrix, speciesnames, sitenames; sitecolumns = true) =
    ComMatrix(occurrences; speciesnames = speciesnames, sitenames = sitenames, sitecolumns = sitecolumns)


function ComMatrix(occs; speciesnames = :auto, sitenames = :auto, sitecolumns = true)
    occurrences = sitecolumns ? occs : occs'

    if sitenames == :auto
        sitenames = ["site$i" for i in 1:size(occurrences, 2)]
    end
    if speciesnames == :auto
        speciesnames = ["species$i" for i in 1:size(occurrences, 1)]
    end

    length(speciesnames) == size(occurrences, 1) || throw(ArgumentError("length of speciesnames ($(length(speciesnames))) different from number of species in occurrences ($(size(occurrences, 1)))"))
    length(sitenames) == size(occurrences, 2) || throw(ArgumentError("length of sitenames ($(length(sitenames))) different from number of sites in occurrences ($(size(occurrences, 2)))"))
    ComMatrix{eltype(occurrences)}(sparse(occurrences), string.(speciesnames), string.(sitenames))
end

function GridData(coords::AbstractMatrix{<:Union{AbstractFloat, Missing}})
    grid = creategrid(coords)
    indices = getindices(coords, grid)
    GridData(indices, grid)
end

GridTopology(min_x, max_x, cellsize_x, min_y, max_y, cellsize_y) =
    GridTopology(range(min_x, stop = max_x, step = cellsize_x), range(min_y, stop = max_y, step = cellsize_y))
