
@testset "Views" begin
    datf = sprand(11,9,0.6)
    datf./=SpatialEcology.colsum(datf)
    com = ComMatrix(datf)

    amphdat = CSV.read("../data/amph_Europe.csv")
    amph = Assemblage(amphdat[4:end], amphdat[1:3], sitecolumns = false)

    vcom = view(com, species = [1,2,3,5,7,8], sites = [2,3,6,8,9])
    comc = copy(vcom)






end
