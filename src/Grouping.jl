#TODO move this to EcoBase

groupsites(a::EcoBase.AbstractAssemblage, s::Symbol; kwargs...) = groupsites(a, a[s]; kwargs...)
groupsites(a::EcoBase.AbstractAssemblage, s::AbstractVector; kwargs...) = [view(a, sites = [(!ismissing(y) && y == x) for y in s]; kwargs...) for x in sort(unique(s))]
groupspecies(a::EcoBase.AbstractAssemblage, s::Symbol; kwargs...) = groupspecies(a, a[s]; kwargs...)
groupspecies(a::EcoBase.AbstractAssemblage, s::AbstractVector; kwargs...) = [view(a, species = [(!ismissing(y) && y == x) for y in s]; kwargs...) for x in sort(unique(s))]
const GroupedAssemblage = Vector{T} where T<:EcoBase.AbstractAssemblage

function show(io::IO, com::GroupedAssemblage)
    println(io, "GroupedAssemblage with $(length(com)) assemblages")
end
