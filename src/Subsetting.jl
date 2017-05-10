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



# creating views
view(occ::AbstractOccFields; species = 1:nspecies(occ), sites = 1:nsites(occ)) = SubOccFields(view(occ.commatrix, sites = sites, species = species), sub(occ.traits,species))
# The SiteFields things are missing as of yet - need to go by the dropbyindex functionality

view(com::AbstractComMatrix; species = 1:nspecies(com), sites = 1:nsites(com)) = SubComMatrix(view(com.occurrences, sites, species))

view(pd::AbstractPointData, sites) = SubPointData(view(pd.coords, sites), view(pd.sitestats, sites))


#I need proper show functions for the views to make it nice.

# TODO Make sure that indices are 1-based! check this with a subset! #NOTE have tried to fix that now by not altering grid - remember that for the copy function!
function view(gd::AbstractGridData, sites)
    indices = view(gd.indices, sites, :)
    #grid = subsetgrid(indices, gd.grid)
    SubGridData(indices, gd.grid, view(gd.sitestats, sites))
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


## TODO Need the copy functions before this is truly useful



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
    println("SubComMatrix with $(records(com)) records of $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end

function show(io::IO, com::SubAssemblage)
    sp = createsummaryline(specnames(com))
    si = createsummaryline(sitenames(com))
    println("SubAssemblage with $(records(com)) records of $(nspecies(com)) species in $(nsites(com)) sites\n\nSpecies names:\n$(sp)\n\nSite names:\n$(si)")
end
