(function(define) {'use strict';
define(function(require) {
   /*
    * bratioLog -- see Algorithm 708
    *
    * Return the logs of the tails of the incomplete beta function
    * of x, a, b.
    *
    * Based on "Significant Digit Computation of the Incomplete Beta Function
    * Ratios", by DiDonato and Morris
    */
   var C, stirlerr, log1p, gam1, algdiv, series, contFrac, logspaceAdd,
       phi, lbeta, erfc, gratR;

   C = require('../constants');
   stirlerr = require('./stirlerr').stirlerr;
   log1p = require('./log1p').log1p;
   gam1 = require('./gam1').gam1;
   phi = require('./phi').phi;
   erfc = require('./erf').erfc;
   lbeta = require('./lbeta').lbeta;
   series = require('../utils').series;
   contFrac = require('../utils').contFrac;
   logspaceAdd = require('../utils').logspaceAdd;
   gratR = require('./gratio').gratR;

   // Page 371 of DiDonato/Morris
   // algdiv(a, b) = lgamma(b) - lgamma(a + b)
   algdiv = (function(){
      var ds;
      ds = [
         0.0833333333333333,
        -0.00277777777760991,
         0.793650666825390e-3,
        -0.595202931351870e-3,
         0.837308034031215e-3,
        -0.125322962780713e-2
      ];
      return function algdiv(a, b) {
         var w, p, q, s;
         p = a / (a + b);
         q = b / (a + b);
         s = 1;
         w = p / b * series(function(i) {
            if (i === 0) { return ds[0]; }
            s += Math.pow(q, 2 * i - 1) + Math.pow(q, 2 * i);
            return ds[i] * s / Math.pow(b, 2 * i);
         }, 6);
         // (26)
         return w - (a + b - 0.5) * Math.log(1 + a / b) - a * (Math.log(b) - 1);
      };
   }());

   // bcorr
   // Computes stirlerr(a) + stirlerr(b) - stirlerr(a + b)
   // see DiDonato (23)
   function bcorr(a, b) {
      if (a > b) { return bcorr(b, a); }
      // a <= b
      return stirlerr(a) + algdiv(a, b) +
               (a + b - 0.5) * Math.log(1 + a / b) +
               a * (Math.log(b) - 1);
   }

   // brcomp -- see p. 370 - 372
   // Computes log(x^a * y^b / B(a, b))
   function brcomp(a, b, x, y) {
      var lambda;

      lambda = a - (a + b) * x;
      if (Math.min(a, b) < 8) {
         return a * Math.log(x) + b * Math.log(y) - lbeta(a, b);
      }
      return 0.5 * (Math.log(a) + Math.log(b) - C.log2pi - Math.log(a + b)) -
             a * phi(1 - lambda / a) - b * phi(1 + lambda / b) -
             bcorr(a, b);
   }

   // Assumes lp is the log of the lower prob
   // Returns a pair lower/upper in log space
   // fromUpper assumes given upper prob in log space
   function fromLower(lp) {
      return { lower: lp, upper: log1p(-Math.exp(lp)) };
   }
   function fromUpper(lq) {
      return { lower: log1p(-Math.exp(lq)), upper: lq };
   }

   // tailFlip exchanges roles of "lower" and "upper"
   function tailFlip(o) {
      return {
         lower: o.upper,
         upper: o.lower
      };
   }

   // bpser (7) Returns the log of the lower probability
   function bpser(a, b, x) {
      var c, ser;
      c = 1;
      ser = a * series(function(i) {
         if (i === 0) { return 0; }
         c = c * (1 - b / i) * x;
         return c / (a + i);
      });
      ser = ser > -1 ? log1p(ser) : -Infinity;
      return -lbeta(a, b) + a * Math.log(x) - Math.log(a) + ser;
   }

   // bup (8)
   // Assumes n >= 1
   /* eslint-disable max-params */
   function bup(a, b, x, y, n) {
      var logSeries;
      logSeries = Math.log(series(function(i, v) {
            if (i === 0) {
               return 1; }
            return v * x * (1 + (b - 1) / (a + i));
         }, n));
      return a * Math.log(x) + b * Math.log(y) - lbeta(a, b) - Math.log(a) +
         logSeries;
   }
   /* eslint-enable max-params */

   // bgrat (9)
   // log space
   // Assumes a > b
   function bgrat(a, b, x) {
      var T, fTs, u, logH, logM, lnx2s, lnx2sn, ps, cs, j, ser;

      T = a + 0.5 * (b - 1);
      fTs = 1 / (4 * T * T);
      u = -T * Math.log(x);
      logH = Math.log(b) + log1p(gam1(b)) + b * Math.log(u) - u;
      logM = logH - (algdiv(b, a) + b * Math.log(T));
      cs = [1];   // cs[n] = 1 / (2n + 1)!
      ps = [1];   // ps[n] = p_n in (9.3)
      lnx2s = Math.log(x) * 0.5;   // (ln(x) / 2) ^ 2
      lnx2s = lnx2s * lnx2s;
      lnx2sn = 1;   // lnx2s ^ n
      j = gratR(b, u, logH);
      ser =
      Math.log(
         series(function(n) {
            var mbmn = -n; // holds m*b - n
            if (n > 0) {
               j = (b + 2 * n - 2) * (b + 2 * n - 2 + 1) * j +
                     (u + b + 2 * n - 2 + 1) * lnx2sn;
               j = j * fTs;
               cs[n] = cs[n - 1] / (2 * n) / (2 * n + 1);
               ps[n] = (b - 1) * cs[n] + 1 / n * series(function(m) {
                  if (m === 0) { return 0; }
                  mbmn = mbmn + b;
                  return mbmn * cs[m] * ps[n - m];
               }, n);
               lnx2sn = lnx2sn * lnx2s;
            }
            return ps[n] * j;
         })
      );
      return logM + ser;
   }

   // bfrac (equation 10, DiDonato)
   // calculating in log space
   function bfrac(a, b, x, y) {
      var lambda;

      lambda = a - (a + b) * x;
      return brcomp(a, b, x, y) + Math.log(contFrac(
         function(i) { // betas from DiDonato
            var n = i - 1;
            return i === 0 ?
                   0 :
                   n + n * (b - n) * x / (a + 2 * n - 1) +
                                         (a + n) / (a + 2 * n + 1) *
                                                   (lambda + 1 + n * (1 + y));
         },
         function(i) { //alphas from DiDonato
            var n = i - 1;
            return i === 1 ?
                   1 :
                   (a + n - 1) * (a + b + n - 1) /
                   (a + 2 * n - 1) / (a + 2 * n - 1) *
                   n * (b - n) * x * x;
         }
      ));
   }

   // basym (11)
   // calculating in log space
   function basym(a, b, x, y) {
      var z, p, q, bg, bgn, Ls, as, bs, cs, es, r, ab, sgn, abnp1,
         root2z, root2znm1, ser;
      p = a / (a + b);
      q = b / (a + b);
      bg = a < b ? Math.sqrt(q / a) : Math.sqrt(p / b); // betaGamma (11.3)
      bgn = 1 / bg;                        // betaGamma^n
      z = Math.sqrt(a * phi(x / p) + b * phi(y / q));   // (11.12)
      Ls = [];
      as = [];
      cs = [];
      es = [];
      ab = a <= b ? a / b : b / a;
      abnp1 = 1;  // (-1)^n * (a/b)^(n+1) or (b/a)^(n+1) for (11.4 - 11.5)
      root2z = Math.sqrt(2) * z;
      root2znm1 = 1; // (sqrt(2) * z) ^ (n - 1) for (11.21)
      sgn = -1; // (-1)^n
      ser = series(function(n) {
               var i, seriesTerms;
               bgn = bgn * bg;
               abnp1 = -abnp1 * ab;
               sgn = -sgn;
               if (n === 0) {
                  es[0] = 1;
                  Ls[0] = Math.sqrt(Math.PI) / 4 * Math.exp(z * z) *
                           erfc(z);
               } else if (n === 1) {
                  as[1] = 2 / 3 * (1 - abnp1) * (a <= b ? q : sgn * p);
                  cs[1] = -0.5 * as[1];
                  es[1] = -cs[1];
                  Ls[1] = Math.pow(2, -1.5);
               } else {
                  root2znm1 = root2znm1 * root2z;
                  Ls[n] = Math.pow(2, -1.5) * root2znm1 +
                           (n - 1) * Ls[n - 2];     // (11.21)
                  as[n] = 2 / (n + 2) * (1 - abnp1) * (a <= b ? q : sgn * p);
                  r = (-n - 1) / 2;
                  seriesTerms = function(j) {
                        if (j === 0) { return 0; }
                        return ((i - j) * r - j) * bs[j] * as[i - j];
                     };
                  bs = [1, r * as[1]]; // first two bs for current r
                  for (i = 2; i <= n; i += 1) {
                     // add the series for b_n^r
                     bs[i] = r * as[i] + 1 / i * series(seriesTerms, i);
                  }
                  cs[n] = bs[n] / (n + 1); // (11.6)
                  es[n] = -series(function(k) {
                     return es[k] * cs[n - k];
                  }, n);
               }
               return es[n] * (Ls[n] * bgn);
             }, 100);
      return Math.log(2 / Math.sqrt(Math.PI)) -
             bcorr(a, b) -   // (11.20)
             z * z + Math.log(ser);
   }

   /* eslint-disable complexity */
   function bratioLog(a, b, x) {
      var p, q, y, bint, bbar, lambda;

      p = a / (a + b);
      q = b / (a + b);
      y = 1 - x;
      bint = Math.floor(b);
      bbar = b - bint;
      if (bbar === 0) {
         bbar += 1;
         bint -= 1;
      }
      lambda = a > b ? (a + b) * y - b :
                       a - (a + b) * x;

      // By convention, if 0 < x < 1 then a == 0 -> result = 1
      // and b == 0 -> result = 0. Can't have a == b == 0.
      if (a === 0) { return { lower: 0, upper: -Infinity }; }
      if (b === 0) { return { lower: -Infinity, upper: 0 }; }

      // R code handles specially the case where max(a, b) < eps * .001
      // but we will not (for now, anyway)

      // This logic follows pages 368-369 of algorithm 708
      /*  CASE where min(a, b) <= 1 */
      if (Math.min(a, b) <= 1) {
         if ( x > 0.5) { return tailFlip(bratioLog(b, a, y)); }
         // we are skipping the case, in R code, for b < min(eps, a*eps)
         // also skipping a < min(eps, b*eps)
         if ( Math.max(a, b) > 1) {
            // either a <= 1 < b or b <= 1 < a
            if (b <= 1) { return fromLower(bpser(a, b, x)); } // (12c)
            // a <= 1 < b
            if (x >= 0.3) { return fromUpper(bpser(b, a, y)); } // (13b)
            // x < 0.3
            if (x < 0.1 && Math.pow(x * b, a) <= 0.7) { return fromLower(bpser(a, b, x)); } // (12d)
            // x >= 0.1 or (bx)^a > 0.7
            // for now, assume bgrat will start with w0 = 0
            if (b > 15) { return fromUpper(bgrat(b, a, y)); } // (14a, 14b)
            // b <= 15
            return fromUpper(logspaceAdd(bup(b, a, y, x, 20),
                                         bgrat(b + 20, a, y))); // (15a, 15b)
         }
         // max(a, b) <= 1
         if (a >= Math.min(0.2, b)) { return fromLower(bpser(a, b, x)); } // (12a)
         // a < min(0.2, b)
         if (Math.pow(x, a) <= 0.9) { return fromLower(bpser(a, b, x)); } // (12b)
         // x^a > 0.9
         if (x >= 0.3) { return fromUpper(bpser(b, a, y)); } // (13a)
         // x < 0.3
         return fromUpper(logspaceAdd(bup(b, a, y, x, 20),
                                         bgrat(b + 20, a, y)));  // 15c
      }
      // min(a, b) > 1
      if (x > p) { return tailFlip(bratioLog(b, a, y)); }
      if (b < 40) {
         if (b * x <= 0.7 || lambda > 650) { return fromLower(bpser(a, b, x)); } // (16a)
         // b * x > 0.7
         if (x <= 0.7) {
            return fromLower(logspaceAdd(bup(bbar, a, y, x, bint),
                                         bpser(a, bbar, x)));
         } // (17a)
         // x > 0.7
         if (a > 15) {
            return fromLower(logspaceAdd(bup(bbar, a, y, x, bint),
                                         bgrat(a, bbar, x)));
         } // (18a)
         // a <= 15
         return fromLower(logspaceAdd(
                              logspaceAdd(bup(bbar, a, y, x, bint),
                                          bup(a, bbar, x, y, 20)),
                              bgrat(a + 20, bbar, x))); // (19a)
      }
      // b >= 40
      if (a <= b) {
         if (a <= 100) { return fromLower(bfrac(a, b, x, y)); } // (20a)
         // a > 100
         if (x < 0.97 * p) { return fromLower(bfrac(a, b, x, y)); } // (20b)
         // x >= 0.97 * p
         return fromLower(basym(a, b, x, y)); // (21a)
      }
      // a > b
      if (b <= 100) { return fromLower(bfrac(a, b, x, y)); } // (20c)
      // b > 100
      if (y > 1.03 * q) { return fromLower(bfrac(a, b, x, y)); } // (20d)
      // y <= 1.03 * q
      return fromLower(basym(a, b, x, y)); // (21c)
   }
   /* eslint-enable complexity */

   /**
    * TODO bratio
    *
    * Return the tails of the incomplete beta function
    * of x, a, b.
    * @memberof basicFunc
    */
   function bratio(a, b, x, lowerTail, logp) {
     var resLog;

     resLog = bratioLog(a, b, x);
     lowerTail = lowerTail !== false;
     logp = logp === true;
     resLog = lowerTail ? resLog.lower : resLog.upper;
     return logp ? resLog : Math.exp(resLog);
   }

   return {
      bratio: bratio
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
