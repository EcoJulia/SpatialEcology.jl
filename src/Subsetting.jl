## Functions for subsetting data objects

subset!(occ::OccFields, species = 1:nspecies(occ), sites = 1:nsites(occ)) = OccFields(occ.commatrix[sites, species], occ.traits[species,:])
subset(occ::OccFields, species = 1:nspecies(occ), sites = 1:nsites(occ)) = subset!(deepcopy(occ), species, sites)

subset!(site::SiteFields, sites = 1:nsites(occ)) = SiteFields(site.coords[sites, :], site.cdtype, site.sitestats[sites, :], site.shape)
subset(site::SiteFields, sites = 1:nsites(occ)) = subset!(deepcopy(site), sites)

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
