#TODO move this to EcoBase

groupsites(a::EcoBase.AbstractAssemblage, s::Symbol) = groupsites(a, a[s])
groupsites(a::EcoBase.AbstractAssemblage, s::AbstractVector) = [view(a, sites = [(!ismissing(y) && y == x) for y in s]) for x in sort(unique(s))]
groupspecies(a::EcoBase.AbstractAssemblage, s::Symbol) = groupspecies(a, a[s])
groupspecies(a::EcoBase.AbstractAssemblage, s::AbstractVector) = [view(a, species = [(!ismissing(y) && y == x) for y in s]) for x in sort(unique(s))]
const GroupedAssemblage = Vector{T} where T<:EcoBase.AbstractAssemblage

function show(io::IO, com::GroupedAssemblage)
    println(io, "GroupedAssemblage with $(length(com)) assemblages")
end
