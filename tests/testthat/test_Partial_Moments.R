# FROM NNS-Python
x <- c(0.6964691855978616, 0.28613933495037946, 0.2268514535642031, 0.5513147690828912, 0.7194689697855631, 0.42310646012446096, 0.9807641983846155, 0.6848297385848633, 0.48093190148436094, 0.3921175181941505, 0.3431780161508694, 0.7290497073840416, 0.4385722446796244, 0.05967789660956835, 0.3980442553304314, 0.7379954057320357, 0.18249173045349998, 0.17545175614749253, 0.5315513738418384, 0.5318275870968661, 0.6344009585513211, 0.8494317940777896, 0.7244553248606352, 0.6110235106775829, 0.7224433825702216, 0.3229589138531782, 0.3617886556223141, 0.22826323087895561, 0.29371404638882936, 0.6309761238544878, 0.09210493994507518, 0.43370117267952824, 0.4308627633296438, 0.4936850976503062, 0.425830290295828, 0.3122612229724653, 0.4263513069628082, 0.8933891631171348, 0.9441600182038796, 0.5018366758843366, 0.6239529517921112, 0.11561839507929572, 0.3172854818203209, 0.4148262119536318, 0.8663091578833659, 0.2504553653965067, 0.48303426426270435, 0.985559785610705, 0.5194851192598093, 0.6128945257629677, 0.12062866599032374, 0.8263408005068332, 0.6030601284109274, 0.5450680064664649, 0.3427638337743084, 0.3041207890271841, 0.4170222110247016, 0.6813007657927966, 0.8754568417951749, 0.5104223374780111, 0.6693137829622723, 0.5859365525622129, 0.6249035020955999, 0.6746890509878248, 0.8423424376202573, 0.08319498833243877, 0.7636828414433382, 0.243666374536874, 0.19422296057877086, 0.5724569574914731, 0.09571251661238711, 0.8853268262751396, 0.6272489720512687, 0.7234163581899548, 0.01612920669501683, 0.5944318794450425, 0.5567851923942887, 0.15895964414472274, 0.1530705151247731, 0.6955295287709109, 0.31876642638187636, 0.6919702955318197, 0.5543832497177721, 0.3889505741231446, 0.9251324896139861, 0.8416699969127163, 0.35739756668317624, 0.04359146379904055, 0.30476807341109746, 0.398185681917981, 0.7049588304513622, 0.9953584820340174, 0.35591486571745956, 0.7625478137854338, 0.5931769165622212, 0.6917017987001771, 0.15112745234808023, 0.39887629272615654, 0.24085589772362448, 0.34345601404832493)
y <- c(0.9290953494701337, 0.3001447577944899, 0.20646816984143224, 0.7712467017344186, 0.179207683251417, 0.7203696347073341, 0.2978651188274144, 0.6843301478774432, 0.6020774780838681, 0.8762070150459621, 0.7616916032270227, 0.6492402854114879, 0.3486146126960078, 0.5308900543442001, 0.31884300700035195, 0.6911215594221642, 0.7845248814489976, 0.8626202294885787, 0.4135895282244193, 0.8672153808700541, 0.8063467153755893, 0.7473209976914339, 0.08726848196743031, 0.023957638562143946, 0.050611236457549946, 0.4663642370285497, 0.4223981453920743, 0.474489623129292, 0.534186315014437, 0.7809131772951494, 0.8198754325768683, 0.7111791151322316, 0.49975889646204175, 0.5018097125708618, 0.7991356578408818, 0.03560152015693441, 0.921601798248779, 0.2733414160633679, 0.7824828518318679, 0.395582605302746, 0.48270235978971854, 0.5931259692926043, 0.2731798106977692, 0.8570159493264954, 0.5319561444631024, 0.1455315278392807, 0.6755524321238062, 0.27625359167650576, 0.2723010177649897, 0.6810977486565571, 0.9493047259244862, 0.807623816061548, 0.9451528088524095, 0.6402025296719795, 0.8258783277528565, 0.6300644920352498, 0.3893090155420259, 0.24163970305689175, 0.18402759570852467, 0.6031603131688895, 0.6566703304734626, 0.21177484928830181, 0.4359435889362071, 0.22965129132316398, 0.13087653733774363, 0.5989734941782344, 0.6688357426448118, 0.8093723729154483, 0.36209409565006223, 0.8513351315065957, 0.6551606487241549, 0.8554790691017261, 0.13596214615618918, 0.10883347378170816, 0.5448015917555307, 0.8728114143337533, 0.6621652225678912, 0.8701363950944805, 0.8453249339337617, 0.6283199211390311, 0.20690841095962864, 0.5176511518958, 0.6448515562981659, 0.42666354124364536, 0.9718610781333566, 0.24973274985042482, 0.05193778223157797, 0.6469719787522865, 0.3698392148054457, 0.8167218997483684, 0.710280810455504, 0.260673487453131, 0.4218711567383805, 0.793490082297006, 0.9398115107412777, 0.7625379749026492, 0.039750173274282985, 0.040137387046519146, 0.16805410857991787, 0.78433600580123)
z <- c(0.19999193561416084,0.6010279101158327,0.9788327513669298,0.8608964619298911,0.7601684508905298,0.12397506746787612,0.5394401401912896,0.8969279890952392,0.3839893553453263,0.5974293052436022,0.06516937735345008,0.15292545930437007,0.533669687225804,0.5430715864428796,0.8676197246411066,0.9298956526581725,0.6460088459791522,0.006548180072424414,0.6025139026895475,0.36841377074834125,0.44801794989436194,0.5048619249681798,0.4000809850582463,0.763740516980946,0.34083865579228434,0.5424284677884146,0.9587984735763967,0.5859672618993342,0.8422555318312421,0.5153219248350965,0.8358609378832195,0.787997995901579,0.2741451405223151,0.6444057500854898,0.02596405447571548,0.2797463018215405,0.10295252828980817,0.4354164588706081,0.26211152577662666,0.6998708543101617,0.37283691796585705,0.3227717548199931,0.1370286323274963,0.8070990185408966,0.7360223497043797,0.34991170542178995,0.9307716779643572,0.8134995545754865,0.32999762541477007,0.7009778150431946,0.9592132203954723,0.285109164298465,0.005404210183425628,0.7840965908154933,0.6534845192821737,0.22306404635944888,0.5599264352651063,0.9126415066887666,0.20749150526588522,0.769668024293192,0.7563728166813091,0.07231316109809582,0.44492578689736473,0.7211553193518122,0.8758657804680099,0.01890807847890197,0.11581293306751883,0.17126277092356368,0.8602241279326432,0.1371855605933343,0.5539492279716964,0.7663649743593801,0.19398868259207802,0.9569799507956978,0.24749785606958874,0.7610819645861326,0.567591973275089,0.7770410669374613,0.0733167994187951,0.845138899921509,0.867602249399254,0.32704688986389774,0.6298085331238098,0.019754547108759235,0.39450735124570824,0.5754821972966637,0.9506549185034494,0.6165089490060033,0.7456130158491189,0.8764042203221318,0.520223244392622,0.8123527374664891,0.8251058874981864,0.6842790562674221,0.4753605948189793,0.7491417107396956,0.4062763059892013,0.5738846393238041,0.32205678990789743,0.5765251949731963)
x_df <- as.data.frame(x)
y_df <- as.data.frame(y)
z_df <- as.data.frame(z)

test_that(
	"LPM", {
	  expect_equal(LPM(0, NULL, x), 0.49, tolerance=1e-5)
	  expect_equal(LPM(0, mean(x), x), 0.49, tolerance=1e-5)
	  expect_equal(LPM(1, mean(x), x), 0.1032933, tolerance=1e-5)
	  expect_equal(LPM(2, mean(x), x), 0.02993767, tolerance=1e-5)
	
	  expect_equal(LPM(0, NULL, x_df), 0.49, tolerance=1e-5)
	  expect_equal(LPM(0, colMeans(x_df), x_df), 0.49, tolerance=1e-5)
	  expect_equal(LPM(1, colMeans(x_df), x_df), 0.1032933, tolerance=1e-5)
	  expect_equal(LPM(2, colMeans(x_df), x_df), 0.02993767, tolerance=1e-5)
	}
)

test_that(
	"UPM", {
	  expect_equal(UPM(0, NULL, x), 0.51, tolerance=1e-5)
	  expect_equal(UPM(0, mean(x), x), 0.51, tolerance=1e-5)
	  expect_equal(UPM(1, mean(x), x), 0.1032933, tolerance=1e-5)
	  expect_equal(UPM(2, mean(x), x), 0.03027411, tolerance=1e-5)
	
	  expect_equal(UPM(0, NULL, x_df), 0.51, tolerance=1e-5)
	  expect_equal(UPM(0, colMeans(x_df), x_df), 0.51, tolerance=1e-5)
	  expect_equal(UPM(1, colMeans(x_df), x_df), 0.1032933, tolerance=1e-5)
	  expect_equal(UPM(2, colMeans(x_df), x_df), 0.03027411, tolerance=1e-5)
	}
)

test_that(
	"Co.UPM", {
		expect_equal(Co.UPM(0, x, y, NULL, NULL), 0.28, tolerance=1e-5)
        expect_equal(Co.UPM(0, x, y, mean(x), mean(y)), 0.28, tolerance=1e-5)
        expect_equal(Co.UPM(1, x, y, mean(x), mean(y)), 0.01204606, tolerance=1e-5)
        expect_equal(Co.UPM(2, x, y, mean(x), mean(y)), 0.0009799173, tolerance=1e-5)

		expect_equal(Co.UPM(0, x_df, y_df, NULL, NULL), 0.28, tolerance=1e-5)
        expect_equal(Co.UPM(0, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.28, tolerance=1e-5)
        expect_equal(Co.UPM(1, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.01204606, tolerance=1e-5)
        expect_equal(Co.UPM(2, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.0009799173, tolerance=1e-5)
	}
)

test_that(
	"Co.LPM", {
        expect_equal(Co.LPM(0, x, y, NULL, NULL), 0.24, tolerance=1e-5)
        expect_equal(Co.LPM(0, x, y, mean(x), mean(y)), 0.24, tolerance=1e-5)
        expect_equal(Co.LPM(1, x, y, mean(x), mean(y)), 0.01058035, tolerance=1e-5)
        expect_equal(Co.LPM(2, x, y, mean(x), mean(y)), 0.0008940764, tolerance=1e-5)

        expect_equal(Co.LPM(0, x_df, y_df, NULL, NULL), 0.24, tolerance=1e-5)
        expect_equal(Co.LPM(0, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.24, tolerance=1e-5)
        expect_equal(Co.LPM(1, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.01058035, tolerance=1e-5)
        expect_equal(Co.LPM(2, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.0008940764, tolerance=1e-5)
	}
)

test_that(
	"D.LPM", {
        expect_equal(D.LPM(0, 0, x, y, NULL, NULL), 0.23, tolerance=1e-5)
        expect_equal(D.LPM(0, 0, x, y, mean(x), mean(y)), 0.23, tolerance=1e-5)
        expect_equal(D.LPM(1, 0, x, y, mean(x), mean(y)), 0.06404049, tolerance=1e-5)
        expect_equal(D.LPM(0, 1, x, y, mean(x), mean(y)), 0.05311669, tolerance=1e-5)
        expect_equal(D.LPM(1, 1, x, y, mean(x), mean(y)), 0.01513793, tolerance=1e-5)
        expect_equal(D.LPM(2, 0, x, y, mean(x), mean(y)), 0.02248309, tolerance=1e-5)
        expect_equal(D.LPM(0, 2, x, y, mean(x), mean(y)), 0.01727327, tolerance=1e-5)
        expect_equal(D.LPM(2, 2, x, y, mean(x), mean(y)), 0.001554909, tolerance=1e-5)

        expect_equal(D.LPM(0, 0, x_df, y_df, NULL, NULL), 0.23, tolerance=1e-5)
        expect_equal(D.LPM(0, 0, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.23, tolerance=1e-5)
        expect_equal(D.LPM(1, 0, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.06404049, tolerance=1e-5)
        expect_equal(D.LPM(0, 1, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.05311669, tolerance=1e-5)
        expect_equal(D.LPM(1, 1, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.01513793, tolerance=1e-5)
        expect_equal(D.LPM(2, 0, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.02248309, tolerance=1e-5)
        expect_equal(D.LPM(0, 2, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.01727327, tolerance=1e-5)
        expect_equal(D.LPM(2, 2, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.001554909, tolerance=1e-5)
	}
)

test_that(
	"D.UPM", {
        expect_equal(D.UPM(0, 0, x, y, NULL, NULL), 0.25, tolerance=1e-5)
        expect_equal(D.UPM(0, 0, x, y, mean(x), mean(y)), 0.25, tolerance=1e-5)
        expect_equal(D.UPM(0, 1, x, y, mean(x), mean(y)), 0.05488706, tolerance=1e-5)
        expect_equal(D.UPM(1, 0, x, y, mean(x), mean(y)), 0.05843498, tolerance=1e-5)
        expect_equal(D.UPM(1, 1, x, y, mean(x), mean(y)), 0.01199175, tolerance=1e-5)
        expect_equal(D.UPM(0, 2, x, y, mean(x), mean(y)), 0.01512857, tolerance=1e-5)
        expect_equal(D.UPM(2, 0, x, y, mean(x), mean(y)), 0.01926167, tolerance=1e-5)
        expect_equal(D.UPM(2, 2, x, y, mean(x), mean(y)), 0.0009941733, tolerance=1e-5)

        expect_equal(D.UPM(0, 0, x_df, y_df, NULL, NULL), 0.25, tolerance=1e-5)
        expect_equal(D.UPM(0, 0, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.25, tolerance=1e-5)
        expect_equal(D.UPM(0, 1, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.05488706, tolerance=1e-5)
        expect_equal(D.UPM(1, 0, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.05843498, tolerance=1e-5)
        expect_equal(D.UPM(1, 1, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.01199175, tolerance=1e-5)
        expect_equal(D.UPM(0, 2, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.01512857, tolerance=1e-5)
        expect_equal(D.UPM(2, 0, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.01926167, tolerance=1e-5)
        expect_equal(D.UPM(2, 2, x_df, y_df, colMeans(x_df), colMeans(y_df)), 0.0009941733, tolerance=1e-5)
	}
)

test_that(
	"LPM.ratio", {
        expect_equal(LPM.ratio(degree=0, target=mean(x), variable=x), 0.49, tolerance=1e-5)
        expect_equal(LPM.ratio(degree=1, target=mean(x), variable=x), 0.5000000000000002, tolerance=1e-5)
        expect_equal(LPM.ratio(degree=2, target=mean(x), variable=x), 0.49720627, tolerance=1e-5)

        expect_equal(LPM.ratio(degree=0, target=colMeans(x_df), variable=x_df), 0.49, tolerance=1e-5)
        expect_equal(LPM.ratio(degree=1, target=colMeans(x_df), variable=x_df), 0.5000000000000002, tolerance=1e-5)
        expect_equal(LPM.ratio(degree=2, target=colMeans(x_df), variable=x_df), 0.49720627, tolerance=1e-5)
	}
)

test_that(
	"UPM.ratio", {
        expect_equal(UPM.ratio(degree=0, target=mean(x), variable=x), 0.51, tolerance=1e-5)
        expect_equal(UPM.ratio(degree=1, target=mean(x), variable=x), 0.4999999999999999, tolerance=1e-5)
        expect_equal(UPM.ratio(degree=2, target=mean(x), variable=x), 0.5027937984146681, tolerance=1e-5)

        expect_equal(UPM.ratio(degree=0, target=colMeans(x_df), variable=x_df), 0.51, tolerance=1e-5)
        expect_equal(UPM.ratio(degree=1, target=colMeans(x_df), variable=x_df), 0.4999999999999999, tolerance=1e-5)
        expect_equal(UPM.ratio(degree=2, target=colMeans(x_df), variable=x_df), 0.5027937984146681, tolerance=1e-5)
	}
)

############################################################################
A <- matrix(c(1,1,3,2,2,3), ncol = 2)
T1 <- matrix(c(1.3333333, 0.6666667, 0.6666667, 0.3333333), ncol=2)
T2 <- matrix(c(0.8888889, 0.4444444, 0.4444444, 0.2222222), ncol=2)
T1_n <- T1
T2_n <- T2
rownames(T1_n) <- c("V1", "V2")
colnames(T1_n) <- c("V1", "V2")
rownames(T2_n) <- c("V1", "V2")
colnames(T2_n) <- c("V1", "V2")

R1 <- NNS::PM.matrix(1,1,colMeans(A), A, pop_adj = TRUE)$cov.matrix
R2 <- NNS::PM.matrix(1,1,colMeans(A), A, pop_adj = FALSE)$cov.matrix
test_that(
  "NNS::PM.matrix - Mean Target", {
    expect_equal(T1, cov(A), tolerance=1e-5)
    expect_equal(R1, T1, tolerance=1e-5)
    expect_equal(R2, T2, tolerance=1e-5)
  }
)

R1 <- NNS::PM.matrix(1,1,NULL, A, pop_adj = TRUE)$cov.matrix
R2 <- NNS::PM.matrix(1,1,NULL, A, pop_adj = FALSE)$cov.matrix
test_that(
  "NNS::PM.matrix - NULL Target", {
    expect_equal(T1, cov(A), tolerance=1e-5)
    expect_equal(R1, T1, tolerance=1e-5)
    expect_equal(R2, T2, tolerance=1e-5)
  }
)

A <- as.data.frame(A)
R1 <- NNS::PM.matrix(1,1,colMeans(A), A, pop_adj = TRUE)$cov.matrix
R2 <- NNS::PM.matrix(1,1,colMeans(A), A, pop_adj = FALSE)$cov.matrix
test_that(
  "NNS::PM.matrix - Mean Target - DataFrame", {
    expect_equal(R1, T1_n, tolerance=1e-5)
    expect_equal(R2, T2_n, tolerance=1e-5)
  }
)

R1 <- NNS::PM.matrix(1,1,NULL, A, pop_adj = TRUE)$cov.matrix
R2 <- NNS::PM.matrix(1,1,NULL, A, pop_adj = FALSE)$cov.matrix
test_that(
  "NNS::PM.matrix - NULL Target - DataFrame", {
    expect_equal(R1, T1_n, tolerance=1e-5)
    expect_equal(R2, T2_n, tolerance=1e-5)
  }
)

#########################################################################
# CDF

# SURVIVAL
A<-c(1,1,2,2,3,3,4,4,5,5,2.5)
T1<-data.table::data.table(matrix(
  c(
    1.0, 1.0, 2.0, 2.0, 2.5, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0, 
    0.8181818, 0.8181818, 0.6363636, 0.6363636, 0.5454545, 0.3636364, 0.3636364, 0.1818182, 0.1818182,0.0000000,0.0000000
  ), 
  ncol=2
))
colnames(T1) <- c("A", "S(x)")
B<-NNS.CDF(A, type="survival")
test_that(
  "NNS.CDF", {
    expect_equal(B$Function, T1, tolerance=1e-5)
    expect_equal(B$target.value, numeric(0), tolerance=1e-5)
  }
)
