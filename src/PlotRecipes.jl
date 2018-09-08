using RecipesBase

RecipesBase.@recipe function f(var::Symbol, asm::SEAssemblage)
    asm.site.sitestats[var], getcoords(places(asm))
end
