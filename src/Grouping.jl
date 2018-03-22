groupsites(a::AbstractAssemblage, s::Symbol) = groupsites(a, a[s])
groupsites(a::AbstractAssemblage, s::AbstractVector) = [view(a, sites = [(!ismissing(y) && y == x) for y in s]) for x in sort(unique(s))]
groupspecies(a::AbstractAssemblage, s::Symbol) = groupspecies(a, a[s])
groupspecies(a::AbstractAssemblage, s::AbstractVector) = [view(a, species = [(!ismissing(y) && y == x) for y in s]) for x in sort(unique(s))]
const GroupedAssemblage = Vector{<:SubAssemblage}

function show(io::IO, com::GroupedAssemblage)
    println(io, "GroupedAssemblage with $(length(com)) SubAssemblages")
end
