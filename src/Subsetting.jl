## Functions for subsetting data objects

#SubDataTypes

# Definition is the same, but importantly this keeps a Named Subarray
type SubComMatrix{T <: Union{Bool, Int}} <: AbstractComMatrix
    occurrences::NamedArrays.NamedArray{T, 2} #Not sure how to specify this is a Named Subarray
end

type SubOccFields{T <: Union{Bool, Int}} <: AbstractOccFields
    commatrix::SubComMatrix{T}
    traits::DataFrames.SubDataFrame
end

type SubAssemblage{S <: SubSiteFields, T <: Union{Bool, Int}} <: AbstractAssemblage # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::S
    occ::SubOccFields{T}
end

type SubGridData <: SubSiteFields
    indices::NamedArrays.NamedMatrix{Int} #SubArray
    grid::GridTopology
    sitestats::DataFrames.SubDataFrame
end

type SubPointData <: SubSiteFields
    coords::NamedArrays.NamedMatrix{Float64} #SubArray
    sitestats::DataFrames.SubDataFrame
end



# creating views
view(occ:OccFields, species = 1:nspecies(occ), sites = 1:nsites(occ)) = subOccFields(view(occ.commatrix, sites, species), sub(occ.traits,species))
# The SiteFields things are missing as of yet - need to go by the dropbyindex functionality

view(com::ComMatrix, sites, species) = SubComMatrix(view(com.occurrences, sites, species))

view(grd::GridData, sites) = SubGridData(view(grd.indices, sites), )




######## OLD #####################
subset!(occ::OccFields, species = 1:nspecies(occ), sites = 1:nsites(occ)) = OccFields(occ.commatrix[sites, species], occ.traits[species,:])
function subset(occ::OccFields, species = 1:nspecies(occ), sites = 1:nsites(occ))
    ret = deepcopy(occ)
    subset!(ret, species, sites)
    ret
end

subset!(site::SiteFields, sites = 1:nsites(occ)) = dropbyindex!(site, sites) #TODO maybe replace dropbyindex with subset! to begin with
function subset(site::SiteFields, sites = 1:nsites(occ))
    ret = deepcopy(site)
    subset!(ret, sites)
    ret
end

subset(asm::Assemblage; species = 1:nspecies(asm), sites = 1:nsites(asm), dropemptyspecies = true, dropemptysites = true) = Assemblage(subset(asm.site, sites), subset(asm.occ, species, sites), dropemptysites = dropemptysites, dropemptyspecies = dropemptyspecies)
subset!(asm::Assemblage; species = 1:nspecies(asm), sites = 1:nsites(asm), dropemptyspecies = true, dropemptysites = true) = Assemblage(subset!(asm.site, sites), subset!(asm.occ, species, sites), dropemptysites = dropemptysites, dropemptyspecies = dropemptyspecies)

subset!(sp::SiteData, sites = 1:nsites(sp)) = SiteData(subset!(sp.site, sites))
subset(sp::SiteData, sites = 1:nsites(sp)) = SiteData(subset(sp.site, sites))






# @enum keepsite allsites occupied
# @enum keepspecies allspecies present
#
# function subset!(occ::OccFields, species = present, sites = occupied)
#     keepsites = false
#     sites == allsites && keepsites = true
#     (keepsites || sites == occupied) && sites = 1:nsites(occ)
#
#     keepspecies = false
#     species == allspecies && keepspecies = true
#     (keepspecies || species == present) && species = 1:nspecies(occ)
#
#     occ.commatrix = occ.commatrix[sites, species]  #There are three big allocations of the same array here - this code should be rewritten to use views instead so most can be done inplace
#     keepspecies || occ.commatrix = occ.commatrix[:, occupancy(occ.commatrix) .> 0]
#     keepsites || occ.commatrix = occ.commatrix[richness(occ.commatrix) .> 0, :]
#
#     #it'd be smarter to do this with indices, I think, as I also need to subset the other bits.
# end
#
#
#
# function subset!(site::SiteFields, sites)
#     site.coords = site.coords[sites, :]
#     site.sitestats = site.sitestats[sites, :]
# end
#
# subset(site::SiteFields, sites) = subset!(deepcopy(site), sites)
#
# function subset!(asm::Assemblage; species = present, sites = occupied)
#     subset!(asm.occ, species, sites)
#     subset!(asm.site, findin(sitenames(asm.site, sitenames(asm.occ))))
# end
