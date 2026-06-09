using RecipesBase

RecipesBase.@recipe function f(var::Symbol, asm::SEAssemblage)
    replace(asm.site.sitestats[!,var], missing => NaN), getcoords(places(asm))
end

RecipesBase.@recipe function f(var::AbstractVector{<:Union{Number, Missing}}, asm::SEAssemblage)
    replace(var, missing => NaN), getcoords(places(asm))
end
