## Functions for subsetting data objects

#SubDataTypes

# Definition is the same, but importantly this keeps a Named Subarray
type SubComMatrix{T <: Union{Bool, Int}} <: AbstractComMatrix{T}
    occurrences::NamedArrays.NamedArray{T, 2} #Not sure how to specify this is a Named Subarray
end

type SubOccFields{T <: Union{Bool, Int}} <: AbstractOccFields{T}
    commatrix::SubComMatrix{T}
    traits::DataFrames.SubDataFrame
end

type SubGridData <: AbstractGridData
    indices::NamedArrays.NamedMatrix{Int} #SubArray
    grid::GridTopology
    sitestats::DataFrames.SubDataFrame
end

type SubPointData <: AbstractPointData
    coords::NamedArrays.NamedMatrix{Float64} #SubArray
    sitestats::DataFrames.SubDataFrame
end

type SubAssemblage{S <: Union{SubGridData, SubPointData}, T <: Union{Bool, Int}} <: AbstractAssemblage # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::S
    occ::SubOccFields{T}
end

type SubSiteData{S <: Union{SubGridData, SubPointData}} <: AbstractSiteData
    site::S
end

# this is because there is issues in NamedArrays when passing a Boolean PooledDataArray, which is useful for subsetting though
asindices{T <: Integer}(x::AbstractArray{T}) = x
asindices{T <: Bool}(x::AbstractArray{T}) = find(x)
# creating views
view(occ::AbstractOccFields; species = 1:nspecies(occ), sites = 1:nsites(occ)) = SubOccFields(view(occ.commatrix, sites = sites, species = species), view(occ.traits,species))
# The SiteFields things are missing as of yet - need to go by the dropbyindex functionality
view(com::AbstractComMatrix; species = 1:nspecies(com), sites = 1:nsites(com)) = SubComMatrix(view(com.occurrences, asindices(sites), asindices(species)))
view(pd::AbstractPointData, sites) = SubPointData(view(pd.coords, sites), view(pd.sitestats, sites))


# In actual fact, only the sitestats DataFrames is subsetted, the rest creates a new object in memory
function view(gd::AbstractGridData, sites)
    indices = gd.indices[sites, :]
    grid = subsetgrid(indices, gd.grid)
    x_shift::Int = (xmin(grid) - xmin(gd.grid)) / xcellsize(gd.grid)
    y_shift::Int = (ymin(grid) - ymin(gd.grid)) / ycellsize(gd.grid)
    indices[:,1] .-= x_shift
    indices[:,2] .-= y_shift
    SubGridData(indices, grid, view(gd.sitestats, sites))
end

view(sp::AbstractSiteData, sites = 1:nsites(sp)) = SubSiteData(view(sp.site, sites))

# the alternative to :occupied and :occurring is :all in both cases
function view(asm::AbstractAssemblage; species = 1:nspecies(asm), sites = 1:nsites(asm), keepsites = :occupied, keepspecies = :occurring)
    occ = view(asm.occ, species = species, sites = sites)
    site = view(asm.site, sites)

    if keepsites == :occupied
        hasspecies = find(richness(occ) .> 0)
        occ = view(occ, sites = hasspecies)
        site = view(site, hasspecies)
    else
        if ! keepsites == :all
            error("keepsites must be :occupied or :all")
        end
    end

    if keepspecies == :occurring
        occ = view(occ, species = find(occupancy(occ) .> 0))
    else
        if ! keepspecies == :all
            error("keepspecies must be :occupied or :all")
        end
    end

    SubAssemblage(site, occ)
end

copy(asm::AbstractAssemblage) = Assemblage(copy(asm.site), copy(asm.occ))
copy(sp::AbstractSiteData) = SiteData(copy(sp.site))
copy(pd::AbstractPointData) = PointData(copy(pd.coords), copy(pd.sitestats))
copy(gd::AbstractGridData) = GridData(copy(gd.indices), gd.grid, my_dataframe_copy(gd.sitestats))
copy(pd::AbstractComMatrix) = ComMatrix(copy(pd.occurrences))
copy(occ::AbstractOccFields) = OccFields(copy(occ.commatrix), my_dataframe_copy(occ.traits))

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
