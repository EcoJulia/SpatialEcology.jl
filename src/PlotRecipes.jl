using RecipesBase

RecipesBase.@recipe function f(var::Symbol, asm::SEAssemblage)
    replace(asm.site.sitestats[var], missing => NaN), getcoords(places(asm))
end
