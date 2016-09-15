

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

guess_xycols(coords)

Assemblage(mam)

methods(Assemblage)

Assemblage(occ, coords)

dat[3]



dat = coords

string.(names(mam))
