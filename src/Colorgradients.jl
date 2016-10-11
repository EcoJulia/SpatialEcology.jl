
#registers some new color gradients for use in Plots
function registercolors()
    PlotUtils.register_gradient_colors(:redblue, [
                colorant"#67001f", colorant"#b2182b", colorant"#d6604d", colorant"#f4a582", colorant"#fddbc7",
                colorant"#f7f7f7", colorant"#d1e5f0", colorant"#92c5de", colorant"#4393c3", colorant"#2166ac",
                colorant"#053061"])

    PlotUtils.register_gradient_colors(:moreland, [
                colorant"#3B4CC0", colorant"#445ACC", colorant"#4D68D7", colorant"#5775E1", colorant"#6282EA",
                colorant"#6C8EF1", colorant"#779AF7", colorant"#82A5FB", colorant"#8DB0FE", colorant"#98B9FF",
                colorant"#A3C2FF", colorant"#AEC9FD", colorant"#B8D0F9", colorant"#C2D5F4", colorant"#CCD9EE",
                colorant"#D5DBE6", colorant"#DDDDDD", colorant"#E5D8D1", colorant"#ECD3C5", colorant"#F1CCB9",
                colorant"#F5C4AD", colorant"#F7BBA0", colorant"#F7B194", colorant"#F7A687", colorant"#F49A7B",
                colorant"#F18D6F", colorant"#EC7F63", colorant"#E57058", colorant"#DE604D", colorant"#D55042",
                colorant"#CB3E38", colorant"#C0282F", colorant"#B40426"])

    PlotUtils.register_gradient_colors(:parula, PlotUtils.sample_evenly([
                colorant"#352A87", colorant"#363093", colorant"#3637A0", colorant"#353DAD", colorant"#3243BA",
                colorant"#2C4AC7", colorant"#2053D4", colorant"#0F5CDD", colorant"#0363E1", colorant"#0268E1",
                colorant"#046DE0", colorant"#0871DE", colorant"#0D75DC", colorant"#1079DA", colorant"#127DD8",
                colorant"#1481D6", colorant"#1485D4", colorant"#1389D3", colorant"#108ED2", colorant"#0C93D2",
                colorant"#0998D1", colorant"#079CCF", colorant"#06A0CD", colorant"#06A4CA", colorant"#06A7C6",
                colorant"#07A9C2", colorant"#0AACBE", colorant"#0FAEB9", colorant"#15B1B4", colorant"#1DB3AF",
                colorant"#25B5A9", colorant"#2EB7A4", colorant"#38B99E", colorant"#42BB98", colorant"#4DBC92",
                colorant"#59BD8C", colorant"#65BE86", colorant"#71BF80", colorant"#7CBF7B", colorant"#87BF77",
                colorant"#92BF73", colorant"#9CBF6F", colorant"#A5BE6B", colorant"#AEBE67", colorant"#B7BD64",
                colorant"#C0BC60", colorant"#C8BC5D", colorant"#D1BB59", colorant"#D9BA56", colorant"#E1B952",
                colorant"#E9B94E", colorant"#F1B94A", colorant"#F8BB44", colorant"#FDBE3D", colorant"#FFC337",
                colorant"#FEC832", colorant"#FCCE2E", colorant"#FAD32A", colorant"#F7D826", colorant"#F5DE21",
                colorant"#F5E41D", colorant"#F5EB18", colorant"#F6F313", colorant"#F9FB0E"], 30))

    PlotUtils.register_gradient_colors(:blueyellow, PlotUtils.sample_evenly([
                colorant"#0707FE", colorant"#1717FC", colorant"#1E1EFA", colorant"#2424F8", colorant"#2828F7",
                colorant"#2C2CF5", colorant"#2F2FF3", colorant"#3232F2", colorant"#3434F0", colorant"#3737EF",
                colorant"#3939EE", colorant"#3B3BEC", colorant"#3D3DEB", colorant"#3F3FEA", colorant"#4141E9",
                colorant"#4242E7", colorant"#4444E6", colorant"#4545E5", colorant"#4747E4", colorant"#4848E3",
                colorant"#4A4AE2", colorant"#4B4BE1", colorant"#4C4CE1", colorant"#4E4EE0", colorant"#4F4FDF",
                colorant"#5050DE", colorant"#5151DD", colorant"#5252DD", colorant"#5454DC", colorant"#5555DB",
                colorant"#5656DA", colorant"#5757DA", colorant"#5858D9", colorant"#5959D8", colorant"#5A5AD8",
                colorant"#5B5BD7", colorant"#5C5CD6", colorant"#5D5DD6", colorant"#5E5ED5", colorant"#5F5FD5",
                colorant"#6060D4", colorant"#6161D4", colorant"#6262D3", colorant"#6262D2", colorant"#6363D2",
                colorant"#6464D1", colorant"#6565D1", colorant"#6666D0", colorant"#6767D0", colorant"#6868D0",
                colorant"#6969CF", colorant"#6969CF", colorant"#6A6ACE", colorant"#6B6BCE", colorant"#6C6CCD",
                colorant"#6D6DCD", colorant"#6E6ECC", colorant"#6E6ECC", colorant"#6F6FCC", colorant"#7070CB",
                colorant"#7171CB", colorant"#7272CA", colorant"#7272CA", colorant"#7373CA", colorant"#7474C9",
                colorant"#7575C9", colorant"#7676C8", colorant"#7676C8", colorant"#7777C8", colorant"#7878C7",
                colorant"#7979C7", colorant"#7979C7", colorant"#7A7AC6", colorant"#7B7BC6", colorant"#7C7CC6",
                colorant"#7C7CC5", colorant"#7D7DC5", colorant"#7E7EC5", colorant"#7F7FC4", colorant"#8080C4",
                colorant"#8080C3", colorant"#8181C3", colorant"#8282C3", colorant"#8282C2", colorant"#8383C2",
                colorant"#8484C2", colorant"#8585C1", colorant"#8585C1", colorant"#8686C1", colorant"#8787C0",
                colorant"#8888C0", colorant"#8888C0", colorant"#8989BF", colorant"#8A8ABF", colorant"#8B8BBF",
                colorant"#8B8BBE", colorant"#8C8CBE", colorant"#8D8DBE", colorant"#8E8EBD", colorant"#8E8EBD",
                colorant"#8F8FBD", colorant"#9090BC", colorant"#9090BC", colorant"#9191BC", colorant"#9292BB",
                colorant"#9393BB", colorant"#9393BB", colorant"#9494BA", colorant"#9595BA", colorant"#9595BA",
                colorant"#9696B9", colorant"#9797B9", colorant"#9898B9", colorant"#9898B8", colorant"#9999B8",
                colorant"#9A9AB8", colorant"#9A9AB7", colorant"#9B9BB7", colorant"#9C9CB6", colorant"#9D9DB6",
                colorant"#9D9DB6", colorant"#9E9EB5", colorant"#9F9FB5", colorant"#9F9FB5", colorant"#A0A0B4",
                colorant"#A1A1B4", colorant"#A2A2B4", colorant"#A2A2B3", colorant"#A3A3B3", colorant"#A4A4B2",
                colorant"#A4A4B2", colorant"#A5A5B2", colorant"#A6A6B1", colorant"#A7A7B1", colorant"#A7A7B0",
                colorant"#A8A8B0", colorant"#A9A9B0", colorant"#A9A9AF", colorant"#AAAAAF", colorant"#ABABAE",
                colorant"#ACACAE", colorant"#ACACAD", colorant"#ADADAD", colorant"#AEAEAD", colorant"#AEAEAC",
                colorant"#AFAFAC", colorant"#B0B0AB", colorant"#B1B1AB", colorant"#B1B1AA", colorant"#B2B2AA",
                colorant"#B3B3A9", colorant"#B3B3A9", colorant"#B4B4A8", colorant"#B5B5A8", colorant"#B5B5A7",
                colorant"#B6B6A7", colorant"#B7B7A6", colorant"#B8B8A6", colorant"#B8B8A5", colorant"#B9B9A5",
                colorant"#BABAA4", colorant"#BABAA4", colorant"#BBBBA3", colorant"#BCBCA3", colorant"#BDBDA2",
                colorant"#BDBDA2", colorant"#BEBEA1", colorant"#BFBFA1", colorant"#BFBFA0", colorant"#C0C09F",
                colorant"#C1C19F", colorant"#C2C29E", colorant"#C2C29E", colorant"#C3C39D", colorant"#C4C49D",
                colorant"#C4C49C", colorant"#C5C59B", colorant"#C6C69B", colorant"#C7C79A", colorant"#C7C799",
                colorant"#C8C899", colorant"#C9C998", colorant"#C9C997", colorant"#CACA97", colorant"#CBCB96",
                colorant"#CCCC95", colorant"#CCCC95", colorant"#CDCD94", colorant"#CECE93", colorant"#CECE92",
                colorant"#CFCF92", colorant"#D0D091", colorant"#D1D190", colorant"#D1D18F", colorant"#D2D28F",
                colorant"#D3D38E", colorant"#D3D38D", colorant"#D4D48C", colorant"#D5D58B", colorant"#D6D68A",
                colorant"#D6D68A", colorant"#D7D789", colorant"#D8D888", colorant"#D8D887", colorant"#D9D986",
                colorant"#DADA85", colorant"#DBDB84", colorant"#DBDB83", colorant"#DCDC82", colorant"#DDDD81",
                colorant"#DDDD80", colorant"#DEDE7F", colorant"#DFDF7E", colorant"#E0E07D", colorant"#E0E07C",
                colorant"#E1E17B", colorant"#E2E27A", colorant"#E2E279", colorant"#E3E377", colorant"#E4E476",
                colorant"#E5E575", colorant"#E5E574", colorant"#E6E672", colorant"#E7E771", colorant"#E8E870",
                colorant"#E8E86E", colorant"#E9E96D", colorant"#EAEA6B", colorant"#EAEA6A", colorant"#EBEB68",
                colorant"#ECEC67", colorant"#EDED65", colorant"#EDED64", colorant"#EEEE62", colorant"#EFEF60",
                colorant"#EFEF5E", colorant"#F0F05C", colorant"#F1F15B", colorant"#F2F259", colorant"#F2F256",
                colorant"#F3F354", colorant"#F4F452", colorant"#F5F550", colorant"#F5F54D", colorant"#F6F64A",
                colorant"#F7F748", colorant"#F7F745", colorant"#F8F841", colorant"#F9F93E", colorant"#FAFA3A",
                colorant"#FAFA36", colorant"#FBFB31", colorant"#FCFC2C", colorant"#FDFD25", colorant"#FDFD1C",
                colorant"#FEFE0D"], 30))


    PlotUtils.register_gradient_colors(:jet, PlotUtils.sample_evenly([
                colorant"#00008F", colorant"#00009F", colorant"#0000AF", colorant"#0000BF", colorant"#0000CF",
                colorant"#0000DF", colorant"#0000EF", colorant"#0000FF", colorant"#0010FF", colorant"#0020FF",
                colorant"#0030FF", colorant"#0040FF", colorant"#0050FF", colorant"#0060FF", colorant"#0070FF",
                colorant"#0080FF", colorant"#008FFF", colorant"#009FFF", colorant"#00AFFF", colorant"#00BFFF",
                colorant"#00CFFF", colorant"#00DFFF", colorant"#00EFFF", colorant"#00FFFF", colorant"#10FFEF",
                colorant"#20FFDF", colorant"#30FFCF", colorant"#40FFBF", colorant"#50FFAF", colorant"#60FF9F",
                colorant"#70FF8F", colorant"#80FF80", colorant"#8FFF70", colorant"#9FFF60", colorant"#AFFF50",
                colorant"#BFFF40", colorant"#CFFF30", colorant"#DFFF20", colorant"#EFFF10", colorant"#FFFF00",
                colorant"#FFEF00", colorant"#FFDF00", colorant"#FFCF00", colorant"#FFBF00", colorant"#FFAF00",
                colorant"#FF9F00", colorant"#FF8F00", colorant"#FF8000", colorant"#FF7000", colorant"#FF6000",
                colorant"#FF5000", colorant"#FF4000", colorant"#FF3000", colorant"#FF2000", colorant"#FF1000",
                colorant"#FF0000", colorant"#EF0000", colorant"#DF0000", colorant"#CF0000", colorant"#BF0000",
                colorant"#AF0000", colorant"#9F0000", colorant"#8F0000", colorant"#800000"], 30))


    PlotUtils.register_gradient_colors(:cube, PlotUtils.sample_evenly([
                colorant"#740081", colorant"#760085", colorant"#770088", colorant"#78008B", colorant"#79008E",
                colorant"#7A0091", colorant"#7B0094", colorant"#7C0097", colorant"#7D009A", colorant"#7E009D",
                colorant"#7F01A0", colorant"#8002A3", colorant"#8004A6", colorant"#8106A9", colorant"#8208AC",
                colorant"#830BB0", colorant"#830DB3", colorant"#8410B6", colorant"#8512B9", colorant"#8514BC",
                colorant"#8517BF", colorant"#851AC2", colorant"#851DC5", colorant"#8520C8", colorant"#8523CB",
                colorant"#8526CD", colorant"#8429D0", colorant"#832CD2", colorant"#832FD4", colorant"#8232D6",
                colorant"#8135D8", colorant"#8038DA", colorant"#803BDC", colorant"#7F3EDE", colorant"#7E40E0",
                colorant"#7E42E2", colorant"#7D45E4", colorant"#7C47E7", colorant"#7C49E9", colorant"#7B4BEB",
                colorant"#794DED", colorant"#784FF0", colorant"#7751F2", colorant"#7653F4", colorant"#7555F6",
                colorant"#7457F7", colorant"#7359F9", colorant"#725BFA", colorant"#715DFC", colorant"#705EFC",
                colorant"#6F60FD", colorant"#6E62FD", colorant"#6D64FD", colorant"#6C66FD", colorant"#6B68FD",
                colorant"#6A6BFC", colorant"#696DFC", colorant"#686FFB", colorant"#6771FB", colorant"#6673FA",
                colorant"#6675F9", colorant"#6577F8", colorant"#6479F7", colorant"#637BF7", colorant"#627DF6",
                colorant"#617EF5", colorant"#6080F4", colorant"#5F82F3", colorant"#5E84F2", colorant"#5D86F1",
                colorant"#5C87F0", colorant"#5B89EF", colorant"#5A8AEE", colorant"#598CEC", colorant"#588EEB",
                colorant"#578FEA", colorant"#5691E8", colorant"#5592E6", colorant"#5394E5", colorant"#5295E3",
                colorant"#5197E2", colorant"#5098E0", colorant"#4F99DE", colorant"#4E9BDD", colorant"#4D9CDB",
                colorant"#4C9ED9", colorant"#4A9FD7", colorant"#49A1D6", colorant"#48A2D4", colorant"#46A4D2",
                colorant"#45A5D0", colorant"#43A7CE", colorant"#42A8CC", colorant"#40AACA", colorant"#3FABC9",
                colorant"#3DADC7", colorant"#3CAEC5", colorant"#3BAFC2", colorant"#3AB1C0", colorant"#39B2BE",
                colorant"#39B3BC", colorant"#38B5BA", colorant"#38B6B8", colorant"#38B7B6", colorant"#38B8B4",
                colorant"#39B9B2", colorant"#39BBB0", colorant"#3ABCAE", colorant"#3ABDAB", colorant"#3BBEA9",
                colorant"#3CBFA7", colorant"#3DC0A5", colorant"#3EC1A3", colorant"#3EC2A0", colorant"#3FC39E",
                colorant"#40C49C", colorant"#41C59A", colorant"#42C697", colorant"#43C795", colorant"#43C893",
                colorant"#44C991", colorant"#45CA8E", colorant"#45CB8C", colorant"#46CB8A", colorant"#47CC87",
                colorant"#47CD85", colorant"#48CE83", colorant"#49CF81", colorant"#49D07E", colorant"#4AD07C",
                colorant"#4BD17A", colorant"#4BD277", colorant"#4CD375", colorant"#4DD373", colorant"#4DD471",
                colorant"#4ED56F", colorant"#4FD66C", colorant"#50D76A", colorant"#50D868", colorant"#51D865",
                colorant"#52D962", colorant"#52DA60", colorant"#53DB5D", colorant"#54DC5A", colorant"#54DD57",
                colorant"#55DE55", colorant"#56DE52", colorant"#57DF50", colorant"#57E04E", colorant"#58E14C",
                colorant"#59E14B", colorant"#5AE24A", colorant"#5BE349", colorant"#5CE349", colorant"#5EE449",
                colorant"#5FE549", colorant"#61E549", colorant"#63E64A", colorant"#65E64A", colorant"#68E74A",
                colorant"#6AE74B", colorant"#6DE84B", colorant"#6FE84C", colorant"#72E94C", colorant"#75E94D",
                colorant"#78EA4E", colorant"#7AEA4E", colorant"#7DEA4F", colorant"#80EB4F", colorant"#82EB50",
                colorant"#85EB50", colorant"#87EB50", colorant"#89EB51", colorant"#8CEB51", colorant"#8EEB52",
                colorant"#91EB52", colorant"#93EB52", colorant"#96EC53", colorant"#98EC53", colorant"#9BEC54",
                colorant"#9DEC54", colorant"#A0EC54", colorant"#A2EC55", colorant"#A5EC55", colorant"#A7EC55",
                colorant"#A9EC56", colorant"#ABEC56", colorant"#ADEC56", colorant"#AFEC57", colorant"#B1EC57",
                colorant"#B4EC57", colorant"#B6EC57", colorant"#B8EC58", colorant"#B9EC58", colorant"#BBEC58",
                colorant"#BDEC58", colorant"#BFEC59", colorant"#C1EC59", colorant"#C3EC59", colorant"#C4EC59",
                colorant"#C6EC59", colorant"#C8EC59", colorant"#C9EC5A", colorant"#CBEC5A", colorant"#CCEC5A",
                colorant"#CDEC5A", colorant"#CFEC5A", colorant"#D0EB5A", colorant"#D1EA5B", colorant"#D2EA5B",
                colorant"#D3E95B", colorant"#D4E85B", colorant"#D5E65B", colorant"#D6E55B", colorant"#D7E45B",
                colorant"#D8E25B", colorant"#D9E15B", colorant"#DAE05C", colorant"#DBDE5C", colorant"#DCDD5C",
                colorant"#DDDB5C", colorant"#DEDA5C", colorant"#DFD95C", colorant"#E0D75C", colorant"#E2D65C",
                colorant"#E3D55D", colorant"#E5D35D", colorant"#E6D25D", colorant"#E7D05D", colorant"#E9CE5D",
                colorant"#EACD5D", colorant"#ECCB5D", colorant"#EDC95E", colorant"#EEC85E", colorant"#EFC65E",
                colorant"#F0C45E", colorant"#F1C25E", colorant"#F2C05E", colorant"#F3BE5E", colorant"#F3BC5E",
                colorant"#F4BA5E", colorant"#F4B85E", colorant"#F5B65E", colorant"#F5B45E", colorant"#F6B25E",
                colorant"#F6B05D", colorant"#F7AD5D", colorant"#F7AB5D", colorant"#F8A85D", colorant"#F8A65D",
                colorant"#F8A35C", colorant"#F9A15C", colorant"#F99E5C", colorant"#F99C5C", colorant"#F9995B",
                colorant"#F9965B"], 30))


    PlotUtils.register_gradient_colors(:blackbody, PlotUtils.sample_evenly([
                colorant"#000000", colorant"#230000", colorant"#340000", colorant"#3C0000", colorant"#3F0100",
                colorant"#400200", colorant"#440500", colorant"#450600", colorant"#480800", colorant"#4A0A00",
                colorant"#4D0C00", colorant"#4E0E00", colorant"#511000", colorant"#531100", colorant"#551300",
                colorant"#561400", colorant"#591600", colorant"#5B1800", colorant"#5C1900", colorant"#5E1A00",
                colorant"#5F1C00", colorant"#621E00", colorant"#641F00", colorant"#662100", colorant"#672200",
                colorant"#692300", colorant"#6A2400", colorant"#6C2600", colorant"#6D2700", colorant"#6F2800",
                colorant"#702A00", colorant"#722B00", colorant"#732C00", colorant"#752D00", colorant"#772F00",
                colorant"#772F00", colorant"#783000", colorant"#7A3100", colorant"#7B3300", colorant"#7D3400",
                colorant"#7D3400", colorant"#7E3500", colorant"#803600", colorant"#813800", colorant"#813800",
                colorant"#833900", colorant"#843A00", colorant"#863B00", colorant"#863B00", colorant"#883D00",
                colorant"#893E00", colorant"#893E00", colorant"#8B3F00", colorant"#8B3F00", colorant"#8C4100",
                colorant"#8E4200", colorant"#8E4200", colorant"#8F4300", colorant"#8F4300", colorant"#914400",
                colorant"#914400", colorant"#924600", colorant"#924600", colorant"#944700", colorant"#944700",
                colorant"#954800", colorant"#954800", colorant"#974900", colorant"#974900", colorant"#994B00",
                colorant"#994B00", colorant"#9A4C00", colorant"#9A4C00", colorant"#9A4C00", colorant"#9C4D00",
                colorant"#9C4D00", colorant"#9D4F00", colorant"#9D4F00", colorant"#9F5000", colorant"#9F5000",
                colorant"#9F5000", colorant"#A05100", colorant"#A05100", colorant"#A25200", colorant"#A25200",
                colorant"#A35400", colorant"#A35400", colorant"#A55500", colorant"#A55500", colorant"#A65600",
                colorant"#A65600", colorant"#A65600", colorant"#A85700", colorant"#A85700", colorant"#AA5900",
                colorant"#AA5900", colorant"#AB5A00", colorant"#AB5A00", colorant"#AD5B00", colorant"#AD5B00",
                colorant"#AE5D00", colorant"#AE5D00", colorant"#B05E00", colorant"#B05E00", colorant"#B15F00",
                colorant"#B15F00", colorant"#B36000", colorant"#B36000", colorant"#B46200", colorant"#B66300",
                colorant"#B66300", colorant"#B76400", colorant"#B76400", colorant"#B96600", colorant"#B96600",
                colorant"#BB6700", colorant"#BB6700", colorant"#BC6800", colorant"#BC6800", colorant"#BE6900",
                colorant"#BF6B00", colorant"#BF6B00", colorant"#C16C00", colorant"#C16C00", colorant"#C26D00",
                colorant"#C46E00", colorant"#C46E00", colorant"#C57000", colorant"#C57000", colorant"#C77100",
                colorant"#C87200", colorant"#C87200", colorant"#CA7400", colorant"#CA7400", colorant"#CC7500",
                colorant"#CD7600", colorant"#CD7600", colorant"#CF7700", colorant"#D07900", colorant"#D07900",
                colorant"#D27A00", colorant"#D37B00", colorant"#D37B00", colorant"#D57C00", colorant"#D67E00",
                colorant"#D67E00", colorant"#D87F00", colorant"#D98000", colorant"#D98000", colorant"#DB8200",
                colorant"#DD8300", colorant"#DD8300", colorant"#DE8400", colorant"#E08500", colorant"#E08500",
                colorant"#E18700", colorant"#E38800", colorant"#E38800", colorant"#E48900", colorant"#E68A00",
                colorant"#E68A00", colorant"#E78C00", colorant"#E98D00", colorant"#E98D00", colorant"#EA8E00",
                colorant"#EC9000", colorant"#EC9000", colorant"#EE9100", colorant"#EF9200", colorant"#F19300",
                colorant"#F19300", colorant"#F29500", colorant"#F49600", colorant"#F49600", colorant"#F59700",
                colorant"#F79900", colorant"#F79900", colorant"#F89A00", colorant"#FA9B00", colorant"#FB9C00",
                colorant"#FB9C00", colorant"#FD9E00", colorant"#FF9F00", colorant"#FF9F00", colorant"#FFA000",
                colorant"#FFA100", colorant"#FFA300", colorant"#FFA300", colorant"#FFA400", colorant"#FFA500",
                colorant"#FFA700", colorant"#FFA700", colorant"#FFA800", colorant"#FFA900", colorant"#FFA900",
                colorant"#FFAA00", colorant"#FFAC00", colorant"#FFAD00", colorant"#FFAD00", colorant"#FFAE00",
                colorant"#FFAF00", colorant"#FFB100", colorant"#FFB200", colorant"#FFB300", colorant"#FFB500",
                colorant"#FFB500", colorant"#FFB600", colorant"#FFB700", colorant"#FFB800", colorant"#FFBB07",
                colorant"#FFBC0A", colorant"#FFBD0E", colorant"#FFBF12", colorant"#FFC015", colorant"#FFC119",
                colorant"#FFC31D", colorant"#FFC524", colorant"#FFC628", colorant"#FFC82B", colorant"#FFCA33",
                colorant"#FFCC36", colorant"#FFCE3D", colorant"#FFCF41", colorant"#FFD248", colorant"#FFD34C",
                colorant"#FFD653", colorant"#FFD85B", colorant"#FFDB62", colorant"#FFDD69", colorant"#FFDF6D",
                colorant"#FFE174", colorant"#FFE47B", colorant"#FFE886", colorant"#FFEA8E", colorant"#FFED95",
                colorant"#FFEF9C", colorant"#FFF0A0", colorant"#FFF3A7", colorant"#FFF6AE", colorant"#FFF8B6",
                colorant"#FFF9B9", colorant"#FFFCC1", colorant"#FFFDC4", colorant"#FFFFCC", colorant"#FFFFCF",
                colorant"#FFFFD3", colorant"#FFFFDA", colorant"#FFFFDE", colorant"#FFFFE1", colorant"#FFFFE5",
                colorant"#FFFFE9", colorant"#FFFFEC", colorant"#FFFFF0", colorant"#FFFFF4", colorant"#FFFFF7",
                colorant"#FFFFFF"], 30))


    PlotUtils.register_gradient_colors(:HMrainbow, PlotUtils.sample_evenly([
                colorant"#000000", colorant"#2D0024", colorant"#38002E", colorant"#3C0031", colorant"#430036",
                colorant"#46003B", colorant"#47003D", colorant"#4B0044", colorant"#4A0049", colorant"#4A004D",
                colorant"#490051", colorant"#470057", colorant"#45015A", colorant"#44025E", colorant"#420361",
                colorant"#3F0666", colorant"#3D076A", colorant"#3A0A6D", colorant"#380C71", colorant"#350F74",
                colorant"#301277", colorant"#2F1479", colorant"#2C177C", colorant"#291B80", colorant"#281C81",
                colorant"#252084", colorant"#222486", colorant"#1D2B89", colorant"#19348A", colorant"#18398B",
                colorant"#183E8D", colorant"#18408E", colorant"#17418E", colorant"#17458F", colorant"#17478E",
                colorant"#17478E", colorant"#17498E", colorant"#174B8E", colorant"#174B8E", colorant"#174E8E",
                colorant"#17508E", colorant"#17508E", colorant"#17528D", colorant"#17558D", colorant"#17558D",
                colorant"#17578C", colorant"#17578C", colorant"#185A8C", colorant"#185A8C", colorant"#185D8B",
                colorant"#185D8B", colorant"#185D8B", colorant"#185D8B", colorant"#18618B", colorant"#18618B",
                colorant"#19658A", colorant"#19658A", colorant"#196889", colorant"#196889", colorant"#196889",
                colorant"#1A6C89", colorant"#1A6C89", colorant"#1B6F88", colorant"#1B6F88", colorant"#1B6F88",
                colorant"#1B7387", colorant"#1B7387", colorant"#1C7686", colorant"#1C7686", colorant"#1D7A85",
                colorant"#1D7A85", colorant"#1D7A85", colorant"#1D7A85", colorant"#1D7D84", colorant"#1D7D84",
                colorant"#1E8083", colorant"#1E8083", colorant"#1F8382", colorant"#1F8382", colorant"#1F8382",
                colorant"#208680", colorant"#208680", colorant"#21897F", colorant"#21897F", colorant"#21897F",
                colorant"#228C7D", colorant"#228C7D", colorant"#238E7B", colorant"#238E7B", colorant"#249179",
                colorant"#249179", colorant"#249179", colorant"#259376", colorant"#259376", colorant"#269674",
                colorant"#269674", colorant"#289871", colorant"#289871", colorant"#299A6F", colorant"#299A6F",
                colorant"#2A9C6C", colorant"#2A9C6C", colorant"#2B9E6A", colorant"#2B9E6A", colorant"#2B9E6A",
                colorant"#2DA068", colorant"#2DA068", colorant"#2EA265", colorant"#2EA265", colorant"#30A463",
                colorant"#30A463", colorant"#32A661", colorant"#32A661", colorant"#33A85F", colorant"#35AA5D",
                colorant"#35AA5D", colorant"#35AA5D", colorant"#37AC5B", colorant"#37AC5B", colorant"#39AE58",
                colorant"#39AE58", colorant"#3BAF56", colorant"#3EB154", colorant"#40B252", colorant"#40B252",
                colorant"#43B450", colorant"#43B450", colorant"#45B54F", colorant"#48B74D", colorant"#48B74D",
                colorant"#48B74D", colorant"#4BB84C", colorant"#4DBA4A", colorant"#50BB49", colorant"#53BD48",
                colorant"#57BE48", colorant"#5BBF47", colorant"#5FC046", colorant"#63C146", colorant"#67C246",
                colorant"#6BC346", colorant"#6FC446", colorant"#6FC446", colorant"#73C446", colorant"#77C546",
                colorant"#7BC546", colorant"#82C647", colorant"#85C747", colorant"#89C748", colorant"#8CC748",
                colorant"#8FC749", colorant"#8FC749", colorant"#93C749", colorant"#96C74A", colorant"#99C74A",
                colorant"#9CC74B", colorant"#A0C84C", colorant"#A7C84E", colorant"#AAC84F", colorant"#ADC84F",
                colorant"#ADC84F", colorant"#B1C850", colorant"#B4C851", colorant"#B7C752", colorant"#BAC752",
                colorant"#BEC753", colorant"#C4C755", colorant"#C7C655", colorant"#C7C655", colorant"#CBC656",
                colorant"#CEC557", colorant"#D4C559", colorant"#D7C45A", colorant"#DAC35B", colorant"#E0C25E",
                colorant"#E0C25E", colorant"#E6C160", colorant"#E9C062", colorant"#ECBE64", colorant"#EEBD68",
                colorant"#F0BC6A", colorant"#F0BC6A", colorant"#F2BB6E", colorant"#F4B972", colorant"#F5B874",
                colorant"#F7B778", colorant"#F8B67B", colorant"#F8B67B", colorant"#FAB57D", colorant"#FBB480",
                colorant"#FCB482", colorant"#FDB485", colorant"#FDB485", colorant"#FEB486", colorant"#FEB38A",
                colorant"#FFB38E", colorant"#FFB391", colorant"#FFB391", colorant"#FFB398", colorant"#FFB4A1",
                colorant"#FFB4A4", colorant"#FFB4A7", colorant"#FFB4A7", colorant"#FFB5A9", colorant"#FFB5AA",
                colorant"#FFB6AD", colorant"#FFB7B0", colorant"#FFB7B0", colorant"#FFB8B3", colorant"#FFB9B3",
                colorant"#FFB9B6", colorant"#FFBAB6", colorant"#FFBAB6", colorant"#FFBBB9", colorant"#FFBCB9",
                colorant"#FFBDBC", colorant"#FFBDBC", colorant"#FFBEBC", colorant"#FFBFBF", colorant"#FFC0BF",
                colorant"#FFC2C2", colorant"#FFC2C2", colorant"#FFC5C5", colorant"#FFC6C6", colorant"#FFC8C8",
                colorant"#FFC9C9", colorant"#FFC9C9", colorant"#FFCACA", colorant"#FFCBCB", colorant"#FFCDCD",
                colorant"#FFCECE", colorant"#FFCECE", colorant"#FFD0D0", colorant"#FFD1D1", colorant"#FFD3D3",
                colorant"#FFD7D7", colorant"#FFD8D8", colorant"#FFD8D8", colorant"#FFDADA", colorant"#FFDBDB",
                colorant"#FFDDDD", colorant"#FFDFDF", colorant"#FFE2E2", colorant"#FFE4E4", colorant"#FFE6E6",
                colorant"#FFE6E6", colorant"#FFE8E8", colorant"#FFEBEB", colorant"#FFEDED", colorant"#FFF0F0",
                colorant"#FFF3F3", colorant"#FFF6F6", colorant"#FFF9F9", colorant"#FFFBFB", colorant"#FFFDFD",
                colorant"#FFFFFF"], 30))
end