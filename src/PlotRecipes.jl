
function convert_to_image(var::AbstractVector, grd::AbstractGridData)
    x = Matrix{Float64}(reverse(cells(grd))...)
    fill!(x, NaN)
    xind, yind =  grd.indices[:,1], grd.indices[:,2] #since matrices are probably drawn from upper left corner
    [x[yind[i], xind[i]] = val for (i, val) in enumerate(var)]
    x
end

RecipesBase.@recipe function f(var::AbstractVector, grd::AbstractGridData)
    seriestype := :heatmap
    aspect_ratio --> :equal
    grid --> false
    convert_to_image(var, grd)
end

RecipesBase.@recipe function f(sit::SiteFields)
    ones(nsites(sit)), sit
end

RecipesBase.@recipe function f(var::AbstractVector, pnt::PointData)
    seriestype := :scatter
    aspect_ratio --> :equal
    grid --> false
    marker_z := var
    legend --> false
    cd = coordinates(pnt)
    cd[:,1], cd[:,2]
end

RecipesBase.@recipe function f(asm::AbstractAssemblage)
    richness(asm), asm.site
end


RecipesBase.@recipe function f(var::AbstractVector, asm::AbstractAssemblage)
    var, asm.site
end
