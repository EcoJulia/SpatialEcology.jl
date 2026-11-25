using SparseArrays
using Random
using SpatialEcology
using Distances
using StableRNGs
using Test

@testset "ComMatrix" begin
    rng = StableRNG(1337)

    datb = sprand(rng, Bool, 12,8,0.9)

    dati = sprand(rng, 8,13,0.9).*100
    dati = round.(Int, dati)

    datf = sprand(rng, 11,9,0.9)

    cmb = ComMatrix(datb)
    @test cmb isa ComMatrix{Bool}
    cmi = ComMatrix(dati,
        speciesnames = ["sp1", "sp2", "sp3", "sp4", "sp5", "sp6", "sp7", "sp8"],
        sitenames = [:a, :b, :c, :d, :e, :f, :g, :h, :i, :j, :k, :l, :m])
    @test cmi isa ComMatrix{<:Union{Int64, Int32}}
    cmf = ComMatrix(datf, sitecolumns = false)

    @test occurrences(ComMatrix(Matrix(datf))) == occurrences(ComMatrix(datf))

    @test occupancy(cmb) == [6, 6, 8, 7, 5, 8, 7, 8, 6, 7, 7, 8]
    @test occupancy(cmi) == [11, 13, 13, 10, 12, 12, 11, 11]
    @test occupancy(cmf) == [9, 10, 11, 11, 9, 9, 11, 9, 11]

    @test richness(cmb) == [9, 10, 12, 11, 11, 9, 11, 10]
    @test richness(cmi) == [8, 8, 7, 7, 7, 7, 6, 7, 8, 6, 8, 7, 7]
    @test richness(cmf) == [8, 9, 9, 8, 8, 7, 9, 8, 9, 8, 7]

    @test nsites(cmb) == 8
    @test nsites(cmi) == 13
    @test nsites(cmf) == 11

    @test nspecies(cmb) == 12
    @test nspecies(cmi) == 8
    @test nspecies(cmf) == 9

    @test occurring(cmb) == 1:12
    @test occurring(cmi) == 1:8
    @test occurring(cmf) == 1:9

    @test occupied(cmb) == 1:8
    @test occupied(cmi) == 1:13
    @test occupied(cmf) == 1:11

    @test noccurring(cmb) == 12
    @test noccurring(cmi) == 8
    @test noccurring(cmf) == 9

    @test noccupied(cmb) == 8
    @test noccupied(cmi) == 13
    @test noccupied(cmf) == 11

    @test getspecies(cmb, 3) == [true, true, true, true, true, true, true, true]
    @test getspecies(cmb, "species2") == [false, true, true, false, true, true, true, true]
    @test getspecies(cmi, 2) == [10, 61, 48, 6, 11, 69, 11, 36, 14, 91, 28, 93, 58]
    gcm = getspecies(cmf, 6)
    @test gcm[3] == 0.9426149574692801
    @test gcm isa SparseVector     # a copy, not a view into the stored transpose

    @test speciesnames(cmi)[8] == "sp8"
    @test speciesnames(cmf)[3] == "species3"

    @test sitenames(cmi)[4] == "d"
    @test sitenames(cmb)[5] == "site5"

    @test length(sitenames(cmf)) == nsites(cmf)
    @test length(sitenames(cmi)) == nsites(cmi)
    @test length(sitenames(cmb)) == nsites(cmb)

    @test length(speciesnames(cmf)) == nspecies(cmf)
    @test length(speciesnames(cmi)) == nspecies(cmi)
    @test length(speciesnames(cmb)) == nspecies(cmb)

    @test sitetotals(cmb) == [9, 10, 12, 11, 11, 9, 11, 10]
    @test sitetotals(cmb) == richness(cmb)
    @test sitetotals(cmi) == [336, 520, 259, 301, 264, 243, 285, 450, 400, 280, 431, 470, 424]
    @test length(sitetotals(cmf)) == nsites(cmf)
    @test sitetotals(cmf)[5] ≈ 3.985356922989832

    @test speciestotals(cmb) == occupancy(cmb)
    @test speciestotals(cmi) == [572, 536, 703, 537, 729, 549, 514, 523]
    @test length(speciestotals(cmf)) == nspecies(cmf)
    @test speciestotals(cmf)[8] ≈ 5.171424049687359

    @test size(cmb) == (12, 8)
    @test size(cmi, 1) == 8

    @test cooccurring(cmb, 1, 3) == [true, false, true, true, true, false, true, true]
    @test cooccurring(cmf, [8, 3]) == [true, true, true, true, true, false, true, true, true, false, true]

    dist = pairwise(BrayCurtis(), view(cmi, sites=1:2))
    @test round.(dist, digits=1) == [0.0 0.3; 0.3 0.0]

    @testset "indexing (DimensionalData)" begin
        SE = SpatialEcology
        DD = SpatialEcology.DD
        # internally consistent: the stored transpose matches, no stored zeros,
        # and one name per row and per column
        consistent(c::ComMatrix) =
            SE.transposed(occurrences(c)) == permutedims(SE._csc(occurrences(c))) &&
            nnz(occurrences(c)) == count(!iszero, Matrix(occurrences(c))) &&
            length(speciesnames(c)) == nspecies(c) && length(sitenames(c)) == nsites(c)
        rowsites(c, i) = findall(!iszero, Matrix(occurrences(c))[i, :])
        rowcounts(c) = vec(count(!iszero, Matrix(occurrences(c)), dims = 2))

        m = Bool[1 0 1; 0 1 1; 1 1 0]
        com = ComMatrix(m, ["a", "b", "c"], ["s1", "s2", "s3"])
        icom = ComMatrix([2 0 5; 0 1 0; 3 0 0], ["a", "b", "c"], ["s1", "s2", "s3"])

        @testset "getindex" begin
            @test com isa AbstractMatrix{Bool}
            @test com[1, 3] === true
            @test com[2, 1] === false
            @test icom[1, 3] === 5
            @test com[At("a"), At("s3")] == com[1, 3]
            @test com[Species = At("c"), Sites = At("s1")] == com[3, 1]

            # a single site drops its dimension, as for an Array, but keeps its name
            col = com[1:2, 1]
            @test col isa DD.AbstractDimVector && !(col isa ComMatrix)
            @test parent(DD.lookup(col, Species)) == ["a", "b"]
            @test parent(DD.lookup(only(DD.refdims(col)))) == ["s1"]
            @test size(com[1:2, 1:1]) == (2, 1)    # a range keeps it

            row1 = com[1:1, :]
            @test row1 isa ComMatrix && speciesnames(row1) == ["a"]
            sub = com[[3, 1], 2:3]
            @test sub isa ComMatrix
            @test speciesnames(sub) == ["c", "a"] && sitenames(sub) == ["s2", "s3"]
            @test Matrix(occurrences(sub)) == m[[3, 1], 2:3]
            @test com[At(["c", "a"]), At(["s2", "s3"])] == sub
            @test speciesnames(com[end:end, :]) == ["c"]
            @test speciesnames(com[[true, false, true], :]) == ["a", "c"]
            @test sitenames(com[:, [false, true, true]]) == ["s2", "s3"]
            whole = com[:, :]
            @test whole == com && speciesnames(whole) == speciesnames(com)
            @test all(consistent, (row1, sub, whole, copy(com), com .> 0))
            @test_throws DD.Lookups.SelectorError com[At("zzz"), :]

            sub[1, 1] = false                   # the result is a copy
            @test com[3, 2] == true

            v = view(com, species = [3, 1], sites = 2:3)
            @test v isa SpatialEcology.SubComMatrix
            @test v[1, 1] == com[3, 2]
            vsub = v[:, 2:2]
            @test vsub isa ComMatrix && !(vsub isa SpatialEcology.SubComMatrix)
            @test speciesnames(vsub) == ["c", "a"] && sitenames(vsub) == ["s3"]
            @test Matrix(occurrences(vsub)) == m[[3, 1], 3:3]
            @test consistent(vsub)
            # keyword views keep both dimensions
            @test size(view(com, species = "a")) == (1, 3)
            @test speciesnames(view(com, species = "a")) == ["a"]
        end

        @testset "setindex!" begin
            b = copy(com)
            b[1, 2] = true                      # inserts a new entry
            @test occupied(b, 1) == rowsites(b, 1) == [1, 2, 3]
            @test noccupied(b, 1) == 3
            @test collect(getspecies(b, 1)) == Matrix(occurrences(b))[1, :]
            @test occupancy(b) == rowcounts(b)
            @test consistent(b)

            n = nrecords(b)
            b[1, 1] = false                     # a zero is removed, not stored
            @test nrecords(b) == n - 1
            @test occupied(b, 1) == rowsites(b, 1) == [2, 3]
            @test occupancy(b) == rowcounts(b)
            @test consistent(b)

            b[At("b"), At("s1")] = true
            @test occupied(b, "b") == [1, 2, 3]
            b[1:2, 1:2] = falses(2, 2)
            @test Matrix(occurrences(b))[1:2, 1:2] == falses(2, 2)
            @test occupancy(b) == rowcounts(b) && consistent(b)

            b .= false                          # broadcasting assignment
            @test nrecords(b) == 0 && occupancy(b) == [0, 0, 0] && consistent(b)
            fill!(occurrences(b), true)
            @test occupancy(b) == [3, 3, 3] && consistent(b)
            occurrences(b) .= false
            @test nrecords(b) == 0 && consistent(b)

            ib = copy(icom)
            ib[1, 1] = 0
            @test occupied(ib, 1) == [3] && consistent(ib)
            ib[2, 3] = 7
            @test occupied(ib, 2) == [2, 3] && speciestotals(ib) == [5, 8, 3]
            @test consistent(ib)
        end

        @testset "every write path of the storage keeps it consistent" begin
            occ(c) = occurrences(c)
            check(c, expected) = Matrix(occ(c)) == expected && consistent(c)
            notm = .!m

            # an array into one entry fails as it does for a SparseMatrixCSC
            c = copy(com)
            @test_throws MethodError (occ(c)[1, 1] = fill(false, 1, 1))
            @test check(c, m)
            c = copy(com); occ(c)[1:2, :] .= false              # one value into a block
            e = copy(m); e[1:2, :] .= false
            @test check(c, e)
            c = copy(com); occ(c)[m] = falses(count(m))         # a logical matrix index
            @test check(c, falses(3, 3))
            c = copy(com); occ(c)[[1, 5]] = [false, false]      # linear indices
            e = copy(m); e[[1, 5]] .= false
            @test check(c, e)

            c = copy(com); copyto!(occ(c), sparse(notm));         @test check(c, notm)
            c = copy(com); copyto!(occ(c), sparsevec(vec(notm))); @test check(c, notm)
            c = copy(com); copyto!(occ(c), notm);                 @test check(c, notm)

            # broadcasting assignment, reading the storage it writes into
            c = copy(com); occ(c) .= occ(c) .& sparse(notm)     # all sparse
            @test check(c, falses(3, 3))
            c = copy(com); occ(c) .= occ(c) .| notm             # sparse and dense
            @test check(c, trues(3, 3))
            # these copyto! methods keep dispatch unambiguous; call them as Base would
            c = copy(com)
            copyto!(occ(c), Broadcast.instantiate(Broadcast.Broadcasted(identity, (notm,))))
            @test check(c, notm)
            c = copy(com)
            copyto!(occ(c), Broadcast.Broadcasted{Nothing}(identity, (notm,), axes(occ(c))))
            @test check(c, notm)
            c = copy(com)
            copyto!(occ(c), Broadcast.Broadcasted{Nothing}(() -> false, (), axes(occ(c))))
            @test check(c, falses(3, 3))

            # changing stored values by other means needs invalidate! (and dropzeros!)
            c = copy(com)
            nonzeros(occ(c)) .= false
            SE.invalidate!(occ(c))
            dropzeros!(occ(c))
            @test check(c, falses(3, 3))

            @test SE.OccMatrix(Matrix(m)) == occ(com)
            @test copy(occ(com)) isa SE.OccMatrix && copy(occ(com)) == occ(com)
            @test sparse(occ(com)) isa SparseMatrixCSC && sparse(occ(com)) == occ(com)
        end

        @testset "axes and colon views" begin
            @test axes(com) == (Base.OneTo(3), Base.OneTo(3))
            @test axes(com, 2) == Base.OneTo(3)
            @test speciesnames(com[begin:begin, :]) == ["a"]
            @test size(view(com, species = :)) == (3, 3)
        end

        @testset "getspecies is a copy" begin
            b = copy(com)
            g = getspecies(b, 1)
            g[2] = true                         # does not write into b
            @test b[1, 2] == false && consistent(b)
            b[1, 2] = true
            @test getspecies(b, 1)[2] == true
        end

        @testset "writing through views" begin
            b = copy(com)
            v = view(b)
            v[2, 1] = true                      # reaches the parent's storage
            @test b[2, 1] == true
            @test occupied(b, 2) == rowsites(b, 2) == [1, 2, 3] && consistent(b)

            w = @view b[2:3, :]
            @test w isa SpatialEcology.SubComMatrix
            @test speciesnames(w) == ["b", "c"] && sitenames(w) == ["s1", "s2", "s3"]
            w2 = view(w, 2:2, :)                # a view of a view: species c
            @test size(w2) == (1, 3) && speciesnames(w2) == ["c"]
            w2[1, 3] = true
            @test b[3, 3] == true && occupied(b, 3) == [1, 2, 3] && consistent(b)
            w2[1, 1] = false
            @test occupied(b, 3) == [2, 3] && occupancy(b) == rowcounts(b) && consistent(b)
            w2 .= false                         # broadcasting through a view
            @test isempty(occupied(b, 3)) && consistent(b)
            fill!(occurrences(view(b, sites = 1:2)), true)
            @test Matrix(occurrences(b))[:, 1:2] == trues(3, 2)
            @test occupancy(b) == rowcounts(b) && consistent(b)
            @test (@view b[1, 2]) isa DD.AbstractDimArray{Bool, 0}
        end

        @test isempty(Test.detect_ambiguities(SpatialEcology; recursive = true))
    end

    @testset "fast occupied / noccupied / thingoccurrences via occurrences_t" begin
        # occupied(com, integer) — fast path uses occ_t column, must agree with slow fallback
        for i in 1:nspecies(cmb)
            fast = occupied(cmb, i)
            slow = findall(!iszero, Matrix(occurrences(cmb))[i, :])
            @test sort(fast) == slow
        end

        # occupied(com, name) — routes through asindices then the fast integer path
        @test sort(occupied(cmb, "species3")) == sort(occupied(cmb, 3))

        # occupied(com, Vector{Int}) — multi-species union (uses nzrows on occ_t view)
        union2 = occupied(cmb, [1, 2])
        @test all(j -> any(!iszero, occurrences(cmb)[[1,2], j]), union2)

        # noccupied fast path — just a length(nzrange), no sites allocation
        for i in 1:nspecies(cmb)
            @test noccupied(cmb, i) == length(occupied(cmb, i))
        end
        @test noccupied(cmb, "species1") == noccupied(cmb, 1)

        # getspecies fast path — column view of occ_t, same values as row of occ
        for i in 1:nspecies(cmb)
            @test collect(getspecies(cmb, i)) == collect(occurrences(cmb)[i, :])
        end

        # occupancy(ComMatrix{Bool}) uses fast colptr diff — must match old rowsum path
        @test occupancy(cmb) == [6, 6, 8, 7, 5, 8, 7, 8, 6, 7, 7, 8]
        @test speciestotals(cmb) == occupancy(cmb)
    end

end
