## Functions for subsetting data objects

#SubDataTypes

# Definition is the same, but importantly this keeps a Subarray
type SubComMatrix{T <: Union{Bool, Int}} <: AbstractComMatrix{T}
    occurrences::SubArray{T,2}
    specnames::SubArray{String,1}
    sitenames::SubArray{String,1}
end

type SubOccFields{T <: Union{Bool, Int}} <: AbstractOccFields{T}
    commatrix::SubComMatrix{T}
    traits::DataFrames.SubDataFrame
end

type SubGridData <: AbstractGridData
    indices::SubArray{Int,2}
    grid::GridTopology
    sitestats::DataFrames.SubDataFrame
end

type SubPointData <: AbstractPointData
    coords::SubArray{Float64,2}
    sitestats::DataFrames.SubDataFrame
end

type SubAssemblage{S <: Union{SubGridData, SubPointData}, T <: Union{Bool, Int}} <: AbstractAssemblage # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::S
    occ::SubOccFields{T}
end

type SubSiteData{S <: Union{SubGridData, SubPointData}} <: AbstractSiteData
    site::S
end

# TODO not sure this is necessary anymore - perhaps remove, or update with a string (for names)
asindices{T <: Integer}(x::AbstractArray{T}) = x
asindices{T <: Bool}(x::AbstractArray{T}) = find(x)
# creating views
view(occ::AbstractOccFields; species = 1:nspecies(occ), sites = 1:nsites(occ)) = SubOccFields(view(occ.commatrix, sites = sites, species = species), view(occ.traits,species))
# The SiteFields things are missing as of yet - need to go by the dropbyindex functionality
function view(com::AbstractComMatrix; species = 1:nspecies(com), sites = 1:nsites(com))
    sit = asindices(sites)
    spec = asindices(species)
    SubComMatrix(view(com.occurrences, sit, spec), view(com.specnames, spec), view(com.sitenames, sit)) #TODO change the order of these in the object to fit the array index order
end

view(pd::AbstractPointData, sites) = SubPointData(view(pd.coords, sites), view(pd.sitestats, sites))

view(gd::AbstractGridData, sites) = SubGridData(view(gd.indices, sites, :), gd.grid, view(gd.sitestats, sites))
view(sp::AbstractSiteData, sites = 1:nsites(sp)) = SubSiteData(view(sp.site, sites))

function view(asm::AbstractAssemblage; species = 1:nspecies(asm), sites = 1:nsites(asm), dropsites = false, dropspecies = false, dropempty = false)
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

    SubAssemblage(site, occ)
end

copy(asm::AbstractAssemblage) = Assemblage(copy(asm.site), copy(asm.occ))
copy(sp::AbstractSiteData) = SiteData(copy(sp.site))
copy(pd::AbstractPointData) = PointData(copy(pd.coords), copy(pd.sitestats))
copy(pd::AbstractComMatrix) = ComMatrix(copy(pd.occurrences), copy(pd.specnames), copy(pd.sitenames))
copy(occ::AbstractOccFields) = OccFields(copy(occ.commatrix), my_dataframe_copy(occ.traits))

function copy(gd::AbstractGridData)
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
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "SubComMatrix with $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, com::SubAssemblage)
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println(io, "SubAssemblage with $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end
