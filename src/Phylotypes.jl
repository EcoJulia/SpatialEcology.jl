


type PhyloFields
    phylo::Phylo.Phylogeny
    nodespecies::Matrix{Bool}

    function PhyloFields(phylo, nodespecies)
        (Ntip(phylo) == size(nodespecies, 1) && Nnode(phylo) == size(nodespecies, 2))
            || throw(DimensionMismatch("Dimension mismatch between nodespecies matrix and phylogeny"))
        new(phylo, nodespecies)
    end
end


type PhyloAssemblage{T} <: SEAssemblage # A type to keep subtypes together, ensuring that they are all aligned at all times
    site::SELocations
    occ::SpeciesData{T}
    phy::PhyloFields

    function PhyloAssemblage(site, occ, phy)
        size(occ.commatrix.occurrences, 2) == size(coordinates(site), 1) || throw(DimensionMismatch("Length mismatch between occurrence matrix and coordinates"))
        Ntip(phy.Phylo) == size(occ.commatrix.occurrences, 1) || throw(DimensionMismatch("Occurrence matrix and phylogeny do not match in species numbers"))
        new(site, occ, phy)
    end
end



function PhyloAssemblage(site::SELocations, occ::SpeciesData, phylo::Phylogeny)
  nodespec = createNodeBySpeciesMatrix(phylo)
  PhyloAssemblage(site, occ, PhyloFields(phylo, nodespec))
end

PhyloAssemblage(assm::Assemblage, phylo::Phylogeny) = PhyloAssemblage(assm.site, assm.occ, phylo)



function createNodeBySpeciesMatrix(tree::Phylogenetics.Phylogeny)
   colname = tree.tipLabel
   rowname = ["$num" for num in nodeNumbers(tree)]
   nodespecies = zeros(Int, nNode(tree), nTip(tree))

   function loc(tree::Phylogenetics.Phylogeny, node::Int)
     if node <= nTip(tree)
       return node
     end
     ret = [loc(tree, descendants(node, tree)[1])..., loc(tree, descendants(node, tree)[2])...]
     nodespecies[node-nTip(tree), ret] = 1
     ret
   end
   loc(tree, basalNode(tree))

   nodespecies
end
