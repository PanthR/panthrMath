var main = require('..');
var chai = require('chai');
var expect = chai.expect;
var precision = 1e-10;
var chisq = main.chisq;
var utils = require('../panthrMath/utils');

/* Rcode generating the tests:
   options(digits = 20)
   x = c(exp(-5:5))             # x >= 0
   a = c(2:10, 20, 80, 200)     # df
   as = rep(a, each = length(x))
   xs = rep(x, length(a))
   xs = signif(xs, digits=7)
   as = signif(as, digits=7)
   chisqs = dchisq(xs, as, log=TRUE)
   ps = pchisq(xs, as, lower.tail = TRUE, log.p=TRUE)
   qs = pchisq(xs, as, lower.tail = FALSE, log.p=TRUE)
   xs = formatC(xs, digits=7, format="g")
   as = formatC(as, digits=7, format="g")
   ps = formatC(ps, digits=15, format="g")
   qs = formatC(qs, digits=15, format="g")
   s = paste(xs, as, chisqs, sep=", ", collapse="],\n[")
   s = paste("[", s, "]", sep="")
   s2 = paste(xs, as, ps, qs, sep=", ", collapse="],\n[")
   s2 = paste("[", s, "]", sep="")
 */

describe('chi squared', function() {
   it('dchisq', function() {
      [
      [0.006737947,        2, -0.696516154059945],
      [0.01831564,        2, -0.702305000559945],
      [0.04978707,        2, -0.718040715559945],
      [0.1353353,        2, -0.760814830559945],
      [0.3678794,        2, -0.877086880559945],
      [       1,        2, -1.19314718055995],
      [2.718282,        2, -2.05228818055995],
      [7.389056,        2, -4.38767518055995],
      [20.08554,        2, -10.7359171805599],
      [54.59815,        2, -27.9922221805599],
      [148.4132,        2, -74.8997471805599],
      [0.006737947,        3, -3.42230750663681],
      [0.01831564,        3, -2.92809632286814],
      [0.04978707,        3, -2.44383205181351],
      [0.1353353,        3, -1.98660612127187],
      [0.3678794,        3, -1.60287828916247],
      [       1,        3, -1.41893853320467],
      [2.718282,        3, -1.77807950165148],
      [7.389056,        3, -3.61346653989908],
      [20.08554,        3, -9.46170845661195],
      [54.59815,        3, -26.2180135335082],
      [148.4132,        3, -72.6255383954223],
      [0.006737947,        4, -6.38966333448416],
      [0.01831564,        4, -5.39545212044684],
      [0.04978707,        4, -4.41118786333756],
      [0.1353353,        4, -3.45396188725429],
      [0.3678794,        4, -2.57023417303548],
      [       1,        4, -1.88629436111989],
      [2.718282,        4, -1.7454352980135],
      [7.389056,        4, -3.0808223745087],
      [20.08554,        4, -8.42906420793444],
      [54.59815,        4, -24.6853693617269],
      [148.4132,        4, -70.5928940855553],
      [0.006737947,        5, -9.52091979516919],
      [0.01831564,        5, -8.0267085508632],
      [0.04978707,        5, -6.54244430769929],
      [0.1353353,        5, -5.08521828607438],
      [0.3678794,        5, -3.70149068974617],
      [       1,        5, -2.51755082187278],
      [2.718282,        5, -1.8766917272132],
      [7.389056,        5, -2.71207884195599],
      [20.08554,        5, -7.5603205920946],
      [54.59815,        5, -23.3166258227834],
      [148.4132,        5, -68.7241504085258],
      [0.006737947,        6, -12.7759576954683],
      [0.01831564,        6, -10.7817464208937],
      [0.04978707,        6, -8.79748219167512],
      [0.1353353,        6, -6.84025612450858],
      [0.3678794,        6, -4.95652864607096],
      [       1,        6, -3.27258872223978],
      [2.718282,        6, -2.131729596027],
      [7.389056,        6, -2.4671167490174],
      [20.08554,        6, -6.81535841586887],
      [54.59815,        6, -22.0716637234539],
      [148.4132,        6, -66.9791881711105],
      [0.006737947,        7, -16.1303577074676],
      [0.01831564,        7, -13.6361464026242],
      [0.04978707,        7, -11.1518821873511],
      [0.1353353,        7, -8.69465607464288],
      [0.3678794,        7, -6.31092871409586],
      [       1,        7, -4.12698873430688],
      [2.718282,        7, -2.48612957654091],
      [7.389056,        7, -2.3215167677789],
      [20.08554,        7, -6.16975835134325],
      [54.59815,        7, -20.9260637358245],
      [148.4132,        7, -65.3335880453953],
      [0.006737947,        8, -19.5677171645606],
      [0.01831564,        8, -16.5735058294487],
      [0.04978707,        8, -13.5892416281209],
      [0.1353353,        8, -10.632015469871],
      [0.3678794,        8, -7.74828822721461],
      [       1,        8, -5.06434819146784],
      [2.718282,        8, -2.92348900214867],
      [7.389056,        8, -2.25887623163426],
      [20.08554,        8, -5.60711773191147],
      [54.59815,        8, -19.863423193289],
      [148.4132,        8, -63.7709473647739],
      [0.006737947,        9, -23.0762678563871],
      [0.01831564,        9, -19.5820564910065],
      [0.04978707,        9, -16.097792303624],
      [0.1353353,        9, -12.6405660998326],
      [0.3678794,        9, -9.25683897506676],
      [       1,        9, -6.0728988833622],
      [2.718282,        9, -3.43203966248984],
      [7.389056,        9, -2.26742693022302],
      [20.08554,        9, -5.11566834721311],
      [54.59815,        9, -18.8719738854869],
      [148.4132,        9, -62.279497918886],
      [0.006737947,       10, -26.6471587061048],
      [0.01831564,       10, -22.6529473104554],
      [0.04978707,       10, -18.6686831370184],
      [0.1353353,       10, -14.7114568876853],
      [0.3678794,       10, -10.82772988081],
      [       1,       10, -7.14378973314767],
      [2.718282,       10, -4.00293048072212],
      [7.389056,       10, -2.3383177867029],
      [20.08554,       10, -4.68655912040586],
      [54.59815,       10, -17.9428647355759],
      [148.4132,       10, -60.8503886308891],
      [0.006737947,       20, -64.7366682579594],
      [0.01831564,       20, -55.7424565596234],
      [0.04978707,       20, -46.75819252564],
      [0.1353353,       20, -37.8009658208905],
      [0.3678794,       20, -28.9172399929212],
      [       1,       20, -20.2332992856809],
      [2.718282,       20, -12.0924397177234],
      [7.389056,       20, -5.42782740618019],
      [20.08554,       20, -2.77606790701183],
      [54.59815,       20, -11.0323742911444],
      [148.4132,       20, -48.9398968055992],
      [0.006737947,       80, -329.361016451248],
      [0.01831564,       80, -290.366802936792],
      [0.04978707,       80, -251.38253973953],
      [0.1353353,       80, -212.425310302283],
      [0.3678794,       80, -173.541591547749],
      [       1,       80, -134.857647483041],
      [2.718282,       80, -96.7167860218921],
      [7.389056,       80, -60.0521760052048],
      [20.08554,       80, -27.4004115088086],
      [54.59815,       80, -5.65672250671654],
      [148.4132,       80, -13.5642367360206],
      [0.006737947,      200, -923.452292385633],
      [0.01831564,      200, -824.458075238937],
      [0.04978707,      200, -725.473813715119],
      [0.1353353,      200, -626.516578812875],
      [0.3678794,      200, -527.632874205213],
      [       1,      200, -428.94892342557],
      [2.718282,      200, -330.808058178038],
      [7.389056,      200, -234.143452751062],
      [20.08554,      200, -141.49167826021],
      [54.59815,      200, -59.7479984856686],
      [148.4132,      200, -7.65549614467123]
      ].forEach(function(tuple) {
         var x, df, logp, rlogp;
         x = tuple[0];
         df = tuple[1];
         rlogp = tuple[2];
         logp = main.dchisq(df, true)(x);
         expect(utils.relativelyCloseTo(rlogp, logp, precision)).to.be.ok;
      });
   });

   it('pchisq, qchisq work', function() {
       [
       [0.006737947,        2, -5.69483119425833,    -0.0033689735],
       [0.01831564,        2, -4.69772253548653,      -0.00915782],
       [0.04978707,        2, -3.70556809507409,     -0.024893535],
       [0.1353353,        2, -2.72679010102147,      -0.06766765],
       [0.3678794,        2, -1.78370779751866,       -0.1839397],
       [       1,        2, -0.932752129567189,             -0.5],
       [2.718282,        2, -0.296899547527277,        -1.359141],
       [7.389056,        2, -0.0251733921859448,        -3.694528],
       [20.08554,        2, -4.35000595872362e-05,        -10.04277],
       [54.59815,        2, -1.39367774945875e-12,       -27.299075],
       [148.4132,        2, -5.92220039868364e-33,         -74.2066],
       [0.006737947,        3, -8.82642463602817, -0.000146812949994786],
       [0.01831564,        3, -7.32989536613004, -0.000655857213246413],
       [0.04978707,        3, -5.83931845109124, -0.0029150705125941],
       [0.1353353,        3, -4.3648467434476, -0.0127981520173338],
       [0.3678794,        3, -2.93360150711825, -0.0546727606927923],
       [       1,        3, -1.61571737153995, -0.221579828439849],
       [2.718282,        3, -0.574705508947227, -0.827526028765191],
       [7.389056,        3, -0.0623847438961791, -2.80546473887133],
       [20.08554,        3, -0.000162965933938784, -8.72205085510083],
       [54.59815,        3, -8.36446051905687e-12, -25.5070292762148],
       [148.4132,        3, -5.79504297583931e-32, -71.9257200830247],
       [0.006737947,        4, -12.0816872084172, -5.66227741079373e-06],
       [0.01831564,        4, -10.0855443031176, -4.16785702594091e-05],
       [0.04978707,        4, -8.09601993354231, -0.000304796096996099],
       [0.1353353,        4, -6.12442548776207, -0.00219114702996106],
       [0.3678794,        4, -4.20112081226959, -0.0150920938891684],
       [       1,        4, -2.40568139136037, -0.0945348918918356],
       [2.718282,        4, -0.931453378672986, -0.500843430271238],
       [7.389056,        4, -0.124092804976371, -2.14813042487487],
       [20.08554,        4, -0.000480466110648461, -7.64099408482555],
       [54.59815,        4, -3.94397911585147e-11, -23.9562458813974],
       [148.4132,        4, -4.45388556503641e-31, -69.886361006938],
       [0.006737947,        5, -15.4362477056362, -1.97752891825957e-07],
       [0.01831564,        5, -12.9403807994439, -2.39918946672461e-06],
       [0.04978707,        5, -10.4516084981723, -2.89021651454481e-05],
       [0.1353353,        5, -7.98207108508574, -0.000341589684158367],
       [0.3678794,        5, -5.56445282188275, -0.00383903625998366],
       [       1,        5, -3.28516983924399, -0.0381528793136381],
       [2.718282,        5, -1.35992860476223, -0.296627432436784],
       [7.389056,        5, -0.21477357700953, -1.64363647568753],
       [20.08554,        5, -0.00120509505549105, -6.72179931794294],
       [54.59815,        5, -1.57901254149125e-10, -22.5690512520977],
       [148.4132,        5, -2.90575839076799e-30, -68.0108583696411],
       [0.006737947,        6, -18.8737275277737, -6.35688395019024e-09],
       [0.01831564,        6, -15.8780676206074, -1.27128486518887e-07],
       [0.04978707,        6, -12.8898594285882, -2.52351404163745e-06],
       [0.1353353,        6, -9.92186519936918, -4.90907145485385e-05],
       [0.3678794,        6, -7.00851524977264, -0.000904559019440828],
       [       1,        6, -4.24138313545577, -0.0144921842182992],
       [2.718282,        6, -1.85331473939311, -0.170452468668726],
       [7.389056,        6, -0.33737709102109, -1.2505044081077],
       [20.08554,        6, -0.00267753212766922, -5.92419822393295],
       [54.59815,        6, -5.587516430094e-10, -21.3053160297755],
       [148.4132,        6, -1.67510406069728e-29, -66.2590924076906],
       [0.006737947,        7, -22.3823718366765, -1.90309040424811e-10],
       [0.01831564,        7, -18.8868729208931, -6.27386705172065e-09],
       [0.04978707,        7, -15.3991034756177, -2.05236395580733e-07],
       [0.1353353,        7, -11.9323094389905, -6.57453863645442e-06],
       [0.3678794,        7, -8.52227909521351, -0.000199005201724408],
       [       1,        7, -5.26459955827847, -0.00518488178225156],
       [2.718282,        7, -2.40554702261225, -0.0945482152029913],
       [7.389056,        7, -0.493518264557687, -0.942826716440797],
       [20.08554,        7, -0.00540241885284034, -5.22360848333197],
       [54.59815,        7, -1.79078771343108e-09, -20.1406102514193],
       [148.4132,        7, -8.74362168865518e-29, -64.6066432121389],
       [0.006737947,        8, -25.9533375794826, -5.35314124096855e-12],
       [0.01831564,        8, -21.9579674471013, -2.90921557498608e-10],
       [0.04978707,        8, -17.9705489751894, -1.56851886596397e-08],
       [0.1353353,        8, -14.0047148884785, -8.27617724406864e-07],
       [0.3678794,        8, -10.097338877843, -4.11898668138316e-05],
       [       1,        8, -6.34721274558465, -0.00175315844086944],
       [2.718282,        8, -3.01135064635544, -0.0504779932295867],
       [7.389056,        8, -0.683783062773108, -0.702599814578043],
       [20.08554,        8, -0.0100677342486407, -4.60344924152783],
       [54.59815,        8, -5.2843293860907e-09, -19.0585201183667],
       [148.4132,        8, -4.20080040423084e-28, -63.0371075242647],
       [0.006737947,        9, -29.5797325824876, -1.42457514646425e-13],
       [0.01831564,        9, -25.0844678087472, -1.27630376150053e-11],
       [0.04978707,        9, -20.5973364691979, -1.13420205485226e-09],
       [0.1353353,        9, -16.1322875657109, -9.85908364355899e-08],
       [0.3678794,        9, -11.7270816860349, -8.0722552810609e-06],
       [       1,        9, -7.48312422004386, -0.000562655563125331],
       [2.718282,        9, -3.66616997695093, -0.0259069379084965],
       [7.389056,        9, -0.90802970836921, -0.516371131654759],
       [20.08554,        9, -0.0175447959871104, -4.05175747019215],
       [54.59815,        9, -1.45268706387683e-08, -18.0472657622415],
       [148.4132,        9, -1.87964133729032e-27, -61.5387165302349],
       [0.006737947,       10, -33.2560350101997, -3.60651302579779e-15],
       [0.01831564,       10, -28.2608580264591, -5.3267866820469e-13],
       [0.04978707,       10, -23.2739659379612, -7.80270200816578e-11],
       [0.1353353,       10, -18.3095711377024, -1.11751739327637e-08],
       [0.3678794,       10, -13.4061721831165, -1.50582230004667e-06],
       [       1,       10, -8.66734403983606, -0.00017213044355067],
       [2.718282,       10, -4.36606979184513, -0.0127824084197093],
       [7.389056,       10, -1.16563577419325, -0.373565958289677],
       [20.08554,       10, -0.0288665856872806, -3.55946913047065],
       [54.59815,       10, -3.75353052890324e-08, -17.0979838912725],
       [148.4132,       10, -7.9024984798616e-27, -60.1026185380933],
       [0.006737947,       20,  -72.03894704141, -5.17467189081044e-32],
       [0.01831564,       20, -62.0442087739489, -1.13381577578287e-27],
       [0.04978707,       20, -52.0585124010722, -2.46192788914637e-23],
       [0.1353353,       20, -42.0973833729915, -5.21601340133494e-19],
       [0.3678794,       20, -32.2029859938102, -1.0337627573297e-14],
       [       1,       20, -22.4895505433391, -1.70967002949505e-10],
       [2.718282,       20, -13.2647231481226, -1.73461984026579e-06],
       [7.389056,       20, -5.3391236853068, -0.00481163260438419],
       [20.08554,       20, -0.602557819044077, -0.792767909230176],
       [54.59815,       20, -4.71412718721406e-05, -9.96238525062239],
       [148.4132,       20, -1.26463700069795e-21, -48.1195018278452],
       [0.006737947,       80, -338.049813731925, -1.5375572229724e-147],
       [0.01831564,       80, -298.055458945009, -3.59880967234968e-130],
       [0.04978707,       80, -258.070811825863, -8.34200572836297e-113],
       [0.1353353,       80, -218.112537903642, -1.88333836994162e-95],
       [0.3678794,       80, -178.225975169179, -3.95769618440266e-78],
       [       1,       80, -138.534260469319, -6.84439591828804e-61],
       [2.718282,       80, -99.3719816241683, -6.97102635785718e-44],
       [7.389056,       80, -61.6468546965104, -1.68698483522864e-27],
       [20.08554,       80, -27.8107453034184, -8.35500621158776e-13],
       [54.59815,       80, -4.31999142530116, -0.0133892346718982],
       [148.4132,       80, -5.2701326315127e-06, -12.1534576635229],
       [0.006737947,      200, -933.057429214767,               -0],
       [0.01831564,      200, -833.063154688737,               -0],
       [0.04978707,      200, -733.078737367898, -4.24570372097217e-319],
       [0.1353353,      200, -633.121078678137, -1.09397942339464e-275],
       [0.3678794,      200, -533.236221690228, -2.62090554099196e-232],
       [       1,      200, -433.549131064834, -5.15234273397188e-189],
       [2.718282,      200, -334.399681918619, -5.91647988789299e-146],
       [7.389056,      200, -236.711371771655, -1.57600443691237e-103],
       [20.08554,      200, -142.992236146372, -7.92977539866464e-63],
       [54.59815,      200, -60.0393854658925, -8.41833482922757e-27],
       [148.4132,      200, -5.99389645664611, -0.00249704263052815]
       ].forEach(function(tuple) {
         var x, df, rp, rq, p, q;
         x = tuple[0];
         df = tuple[1];
         rp = tuple[2];
         rq = tuple[3];
         p = main.pchisq(df, true, true)(x);
         q = main.pchisq(df, false, true)(x);
         if (rp > -30 && rp < -1e-10) {
            if (rp < -1e-3) {
               expect(utils.relativelyCloseTo(p, rp, precision)).to.be.ok;
            }
            expect(utils.relativelyCloseTo(Math.exp(p), Math.exp(rp), precision)).to.be.ok;
            expect(utils.relativelyCloseTo(x, main.qchisq(df, true, true)(p), 1e-8)).to.be.ok;
         }
         if (rq > -30 && rq < -1e-10) {
            if (rq < -1e-3) {
               expect(utils.relativelyCloseTo(q, rq, precision)).to.be.ok;
            }
            expect(utils.relativelyCloseTo(Math.exp(q), Math.exp(rq), precision)).to.be.ok;
            expect(utils.relativelyCloseTo(x, main.qchisq(df, false, true)(q), 1e-6)).to.be.ok;
         }
      });
   });
   it('chisq also exported as an object', function() {
      var o;
      o = main.chisq(.1);
      ['d', 'p', 'q', 'r'].forEach(function(s) {
         expect(o).to.respondTo(s);
      });
   });
   /* Rcode
      options(digits = 20)
      x = c(NaN, -Inf, -1.3, 0, 2, Inf)
      df = c(NaN, -Inf, -1.3, 0, 1.5, Inf)
      g = expand.grid(x=x, df=df)
      g$d = dchisq(g$x, g$df, log=TRUE)
      g$p = pchisq(g$x, g$df, log.p=TRUE)
      g$q = pchisq(g$x, g$df, lower.tail=FALSE, log.p=TRUE)
      s = paste(g$x, g$df, g$d, g$p, g$q, sep=", ", collapse="],\n[")
      s = paste("[", s, "]", sep="")
   */
   it('pchisq, dchisq handle inappropriate inputs', function() {
      [
      [NaN, NaN, NaN, NaN, NaN],
      [-Infinity, NaN, NaN, NaN, NaN],
      [-1.3, NaN, NaN, NaN, NaN],
      [0, NaN, NaN, NaN, NaN],
      [2, NaN, NaN, NaN, NaN],
      [Infinity, NaN, NaN, NaN, NaN],
      [NaN, -Infinity, NaN, NaN, NaN],
      [-Infinity, -Infinity, NaN, NaN, NaN],
      [-1.3, -Infinity, NaN, NaN, NaN],
      [0, -Infinity, NaN, NaN, NaN],
      [2, -Infinity, NaN, NaN, NaN],
      [Infinity, -Infinity, NaN, NaN, NaN],
      [NaN, -1.3, NaN, NaN, NaN],
      [-Infinity, -1.3, NaN, NaN, NaN],
      [-1.3, -1.3, NaN, NaN, NaN],
      [0, -1.3, NaN, NaN, NaN],
      [2, -1.3, NaN, NaN, NaN],
      [Infinity, -1.3, NaN, NaN, NaN],
      [NaN, 0, NaN, NaN, NaN],
      [-Infinity, 0, -Infinity, -Infinity, 0],
      [-1.3, 0, -Infinity, -Infinity, 0],
      [0, 0, Infinity, -Infinity, 0],
      [2, 0, -Infinity, 0, -Infinity],
      [Infinity, 0, -Infinity, 0, -Infinity],
      [NaN, 1.5, NaN, NaN, NaN],
      [-Infinity, 1.5, -Infinity, -Infinity, 0],
      [-1.3, 1.5, -Infinity, -Infinity, 0],
      [0, 1.5, Infinity, -Infinity, 0],
      [2, 1.5, -1.89642813199124, -0.301132078917963, -1.34699684526317],
      [Infinity, 1.5, -Infinity, 0, -Infinity],
      [NaN, Infinity, NaN, NaN, NaN],
      [-Infinity, Infinity, -Infinity, -Infinity, 0],
      [-1.3, Infinity, -Infinity, -Infinity, 0],
      [0, Infinity, -Infinity, -Infinity, 0],
      [2, Infinity, NaN, NaN, NaN],
      [Infinity, Infinity, -Infinity, 0, -Infinity]
      ].forEach(function(tuple) {
         var x, df, rlogp, rlogq, rlogd, logp, logq, logd;
         x = tuple[0];
         df = tuple[1];
         rlogd = tuple[2];
         rlogp = tuple[3];
         rlogq = tuple[4];
         logp = main.pchisq(df, true, true)(x);
         logq = main.pchisq(df, false, true)(x);
         logd = main.dchisq(df, true)(x);
         expect(utils.relativelyCloseTo(logp, rlogp)).to.equal(true);
         expect(utils.relativelyCloseTo(logq, rlogq)).to.equal(true);
         expect(utils.relativelyCloseTo(logd, rlogd)).to.equal(true);
      });
   });

   /* Rcode generating the tests:
      options(digits = 20)
      p = c(NaN, -Inf, -1, 0, 0.3, 1, Inf)
      df = c(NaN, -Inf, -1.3, 0, 1.5, Inf)
      g = expand.grid(p=p, df=df)
      g$x1 = qchisq(g$p, g$df, lower.tail=TRUE, log.p=FALSE)
      g$x2 = qchisq(g$p, g$df, lower.tail=FALSE, log.p=FALSE)
      s = paste(g$p, g$df, g$x1, g$x2, sep=", ", collapse="],\n[")
      s = paste("[", s, "]", sep="")
   */
   it('qchisq handles inappropriate inputs', function() {
      [
      [NaN, NaN, NaN, NaN],
      [-Infinity, NaN, NaN, NaN],
      [-1, NaN, NaN, NaN],
      [0, NaN, NaN, NaN],
      [0.3, NaN, NaN, NaN],
      [1, NaN, NaN, NaN],
      [Infinity, NaN, NaN, NaN],
      [NaN, -Infinity, NaN, NaN],
      [-Infinity, -Infinity, NaN, NaN],
      [-1, -Infinity, NaN, NaN],
      [0, -Infinity, 0, Infinity],
      [0.3, -Infinity, NaN, NaN],
      [1, -Infinity, Infinity, 0],
      [Infinity, -Infinity, NaN, NaN],
      [NaN, -1.3, NaN, NaN],
      [-Infinity, -1.3, NaN, NaN],
      [-1, -1.3, NaN, NaN],
      [0, -1.3, 0, Infinity],
      [0.3, -1.3, NaN, NaN],
      [1, -1.3, Infinity, 0],
      [Infinity, -1.3, NaN, NaN],
      [NaN, 0, NaN, NaN],
      [-Infinity, 0, NaN, NaN],
      [-1, 0, NaN, NaN],
      [0, 0, 0, Infinity],
      [0.3, 0, 0, 0],
      [1, 0, Infinity, 0],
      [Infinity, 0, NaN, NaN],
      [NaN, 1.5, NaN, NaN],
      [-Infinity, 1.5, NaN, NaN],
      [-1, 1.5, NaN, NaN],
      [0, 1.5, 0, Infinity],
      [0.3, 1.5, 0.401589210126614, 1.75379084277085],
      [1, 1.5, Infinity, 0],
      [Infinity, 1.5, NaN, NaN],
      [NaN, Infinity, NaN, NaN],
      [-Infinity, Infinity, NaN, NaN],
      [-1, Infinity, NaN, NaN],
      [0, Infinity, 0, Infinity],
      [0.3, Infinity, Infinity, Infinity],
      [1, Infinity, Infinity, 0],
      [Infinity, Infinity, NaN, NaN]
      ].forEach(function(tuple) {
         var p, df, x1r, x2r, x1, x2;
         p = tuple[0];
         df = tuple[1];
         x1r = tuple[2];
         x2r = tuple[3];
         x1 = main.qchisq(df, true, false)(p);
         x2 = main.qchisq(df, false, false)(p);
         expect(utils.relativelyCloseTo(x1r, x1)).to.equal(true);
         expect(utils.relativelyCloseTo(x2r, x2)).to.equal(true);
      });
   });
});


