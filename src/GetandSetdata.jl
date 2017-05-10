coordinates(sd::SpatialData) = coordinates(sd.site)
coordinates(pd::AbstractPointData) = pd.coords

function coordinates(gd::AbstractGridData)
    index1 = xrange(gd.grid)[gd.indices[:,1]]
    index2 = yrange(gd.grid)[gd.indices[:,2]]
    NamedArrays.NamedArray(hcat(index1, index2), NamedArrays.names(gd.indices), dimnames(gd.indices))
end

traits(occ::AbstractOccFields) = occ.traits
traits(asm::Assmbl) = traits(asm.occ)

sitestats(asm::Assmbl) = asm.site.sitestats

traitnames(asm::Assmbl) = names(traits(asm))
sitestatnames(asm::Assmbl) = names(sitestats(asm))

function addtraits!(asm::Assemblage, newtraits::DataFrames.DataFrame, species::Symbol; validate = true, tolerance = 0.5)
    if validate
        dif, left, right = setdiff(newtraits[species], specnames(asm)) , size(newtraits,1), nspecies(asm)
        max(dif/left, dif/right) == 0 && error("No match between species names, aborting join")
        println("$dif matching species names,\n",
                "\t$(100*dif/left)% of $left species in the Assemblage\n",
                "\t$(100*dif/right)% of $right species in the new traits data\n")
        max(dif/left, dif/right) < tolerance && error("Aborting join, as fit was smaller than the tolerance of $tolerance . To perform the join decrease the tolerance value")
    end
    #ugly workaround for the join
    nam = names(newtraits)
    nam[findin(nam, species)] = :species
    asm.occ.traits = join(asm.occ.traits, newtraits, kind = :left, on = :species)
end

function addsitestats!(asm::Assemblage, newsites::DataFrames.DataFrame, sites::Symbol; validate = true, tolerance = 0.5)
    if validate
        dif, left, right = setdiff(newsites[sites], sitenames(asm)) , size(newsites,1), nsites(asm)
        max(dif/left, dif/right) == 0 && error("No match between site names, aborting join")
        println("$dif matching site names,\n",
                "\t$(100*dif/left)% of $left sites in the Assemblage\n",
                "\t$(100*dif/right)% of $right sites in the new sitestats data\n")
        max(dif/left, dif/right) < tolerance && error("Aborting join, as fit was smaller than the tolerance of $tolerance . To perform the join decrease the tolerance value")
    end
    #ugly workaround for the join
    nam = names(newsites)
    nam[findin(nam, sites)] = :sites
    asm.occ.traits = join(asm.occ.traits, newtraits, kind = :left, on = :sites)
end
