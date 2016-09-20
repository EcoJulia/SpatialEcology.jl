

using DataFrames
mam = readtable("src/mam.csv")
for i in 4:10
    mam[i] = round(rand(10))
end


mam[3]
mam[3] = map(x->"$x", mam[3])

BenHoltMatrix(mam)

occ, coords = parsesingleDataFrame(mam)

occ = parseDataFrame(occ)

occ = ComMatrix(occ)

specnames(occ)

occurring = find(occupancy(occ) .> 0)

occ = occ[:, occurring]

traits

nspecies(occ)
size(occ)

guess_xycols(coords)

Assemblage(mam)

specnames(occ)

methods(OccFields)
OccFields(occ)
