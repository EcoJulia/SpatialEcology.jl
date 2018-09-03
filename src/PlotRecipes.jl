
function convert_to_image(var::AbstractVector, grd::AbstractGridData)
    x = Matrix{Float64}(undef, reverse(cells(grd))...)
    fill!(x, NaN)
    xind, yind =  grd.indices[:,1], grd.indices[:,2] #since matrices are drawn from upper left corner
    [x[yind[i], xind[i]] = val for (i, val) in enumerate(var)]
    x
end

RecipesBase.@recipe function f(var::AbstractVector, grd::AbstractGridData)
    seriestype := :heatmap
    aspect_ratio --> :equal
    grid --> false
    framestyle --> :none
    xrange(grd), yrange(grd), convert_to_image(var, grd)
end

RecipesBase.@recipe function f(sit::SiteFields)
    ones(nsites(sit)), sit
end

RecipesBase.@recipe function f(var::AbstractVector, pnt::AbstractPointData)
    seriestype := :scatter
    aspect_ratio --> :equal
    grid --> false
    marker_z := var
    legend --> false
    framestyle --> :none
    cd = coordinates(pnt)
    cd[:,1], cd[:,2]
end

RecipesBase.@recipe function f(asm::AbstractAssemblage; occupied = true)
    var = richness(asm)
    if occupied
        var = [Float64(v) for v in var]
        (var[var.==0] .= NaN)
    end
    var, asm.site
end


RecipesBase.@recipe function f(var::AbstractVector, asm::AbstractAssemblage)
    var, asm.site
end

RecipesBase.@recipe function f(var::Symbol, asm::AbstractAssemblage)
    recode(asm.site.sitestats[var], Missings.missing => NaN), asm.site
end

RecipesBase.@recipe function f(g::Function, asm::AbstractAssemblage)
    g(asm), asm.site
end
