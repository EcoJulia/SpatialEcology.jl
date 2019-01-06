## Functions for subsetting data objects

#SubDataTypes

# Definition is the same, but importantly this keeps a Subarray
mutable struct SubComMatrix{D <: Real} <: AbstractComMatrix{D}
    occurrences::SubArray{D,2}
    speciesnames::SubArray{String,1}
    sitenames::SubArray{String,1}
end

mutable struct SubSpeciesData{D <: Real} <: SEThings{D}
    commatrix::SubComMatrix{D}
    traits::DataFrames.SubDataFrame
end

mutable struct SubGridData <: SEGrid
    indices::SubArray{Int,2}
    grid::GridTopology
end

mutable struct SubPointData
    coords::SubArray{Float64,2}
end

mutable struct SubLocations{T<:Union{SubGridData, SubPointData}} <: SELocations
    coords::T
    sitestats::DataFrames.SubDataFrame
end

mutable struct SubAssemblage{D <: Real, P <: SubLocations} <: SEAssemblage{D, SubSpeciesData{D}, P} # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::P
    occ::SubSpeciesData{D}
end

# TODO delete
mutable struct SubSiteData{S} <: SESpatialData where S <: SubLocations
    site::S
end

# TODO not sure this is necessary anymore - perhaps remove, or update with a string (for names)

asindices(x::AbstractArray{T}) where T <: Integer = x
asindices(x::AbstractArray{T}) where T <: Bool = findall(x)
asindices(x, y) = asindices(x)
asindices(x::AbstractArray{T}, y::AbstractArray{T}) where T <: AbstractString = indexin(x, y)
# creating views
view(occ::SEThings; species = 1:nspecies(occ), sites = 1:nsites(occ)) = SubSpeciesData(view(occ.commatrix, sites = sites, species = species), view(occ.traits,species, :))
# The SELocations things are missing as of yet - need to go by the dropbyindex functionality
function view(com::AbstractComMatrix; species = 1:nspecies(com), sites = 1:nsites(com))
    sit = asindices(sites, sitenames(com))
    spec = asindices(species, speciesnames(com))
    SubComMatrix(view(com.occurrences, spec, sit), view(com.speciesnames, spec), view(com.sitenames, sit)) #TODO change the order of these in the object to fit the array index order
end

view(pd::SEPointData, sites) = SubPointData(view(pd.coords, sites, :))
view(gd::SEGrid, sites) = SubGridData(view(gd.indices, sites, :), gd.grid)
view(lo::SELocations, sites) = SubLocations{SubGridData}(view(lo.coords, sites), view(lo.sitestats, sites, :))
view(sp::SESpatialData, sites = 1:nsites(sp)) = SubSiteData(view(sp.site, sites))

 function view(asm::SEAssemblage{D, P}; species = 1:nspecies(asm), sites = 1:nsites(asm), dropsites = false, dropspecies = false, dropempty = false) where D where P
    occ = view(asm.occ, species = species, sites = sites)
    site = view(asm.site, sites)

    if dropsites || dropempty
        hasspecies = occupied(occ)
        occ = view(occ, sites = hasspecies)
        site = view(site, hasspecies)
    end

    if dropspecies || dropempty
        occ = view(occ, species = occurring(occ))
    end

    SubAssemblage{D, typeof(site)}(site, occ)
end

Assemblage(assm::SubAssemblage) = copy(assm)

copy(asm::SEAssemblage) = Assemblage(copy(asm.site), copy(asm.occ))
copy(sp::SESpatialData) = SiteData(copy(sp.site))
copy(pd::SEPointData) = PointData(copy(pd.coords), copy(pd.sitestats))
copy(pd::AbstractComMatrix) = ComMatrix(copy(pd.occurrences), copy(pd.speciesnames), copy(pd.sitenames))
copy(occ::SEThings) = SpeciesData(copy(occ.commatrix), my_dataframe_copy(occ.traits))

function copy(gd::SEGrid)
    indices = copy(gd.indices)
    grid = subsetgrid(indices, gd.grid)
    x_shift = Int.((xmin(grid) - xmin(gd.grid)) / xcellsize(gd.grid))
    y_shift = Int.((ymin(grid) - ymin(gd.grid)) / ycellsize(gd.grid))
    indices[:,1] .-= x_shift
    indices[:,2] .-= y_shift
    GridData(indices, grid, my_dataframe_copy(gd.sitestats))
end

# because I cannot define a new copy method for DataFrames
function my_dataframe_copy(sdf::AbstractDataFrame)
    ret = DataFrame()
    for n in names(sdf)
        ret[n] = sdf[n]
    end
    ret
end
# Helper functions

function subsetgrid(indices, grid)
   xmin = xrange(grid)[minimum(indices[:,1])]
   ymin = yrange(grid)[minimum(indices[:,2])]
   xcells = maxrange(indices[:,1]) + 1
   ycells = maxrange(indices[:,2]) + 1
   GridTopology(xmin, grid.xcellsize, xcells, ymin, grid.ycellsize, ycells)
end

function show(io::IO, com::SubComMatrix)
    sp = createsummaryline(speciesnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "SubComMatrix with $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, com::SubAssemblage)
    sp = createsummaryline(speciesnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "SubAssemblage with $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end
