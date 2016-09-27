import RecipesBase




# Make this a type recipe

convert_to_image{G <: GridData, T <: Any}(var::AbstractVector, asm::Assemblage{G, T}) = convert_to_image(var, asm.site)

function convert_to_image(var::AbstractVector, grd::GridData)
    x = Matrix{Float64}(reverse(cells(grd))...)
    fill!(x, NaN)
    xind, yind =  grd.indices[:,1], grd.indices[:,2] #since matrices are probably drawn from upper left corner
    [x[yind[i], xind[i]] = val for (i, val) in enumerate(var)]
    x
end

x = convert_to_image(richness(mam), mam.site)
heatmap(x, aspect_ratio = :equal, grid = false) #so far only works in plotlyjs()
