(function(define) {'use strict';
define(function(require) {

   /*
    * bratioLog -- see Algorithm 708
    *
    * Return the logs of the tails of the incomplete beta function
    * of x, a, b.
    *
    * Based on "Significant Digit Computation of the Incomplete Beta Function Ratios",
    * DiDonato, Morris
    */
   var C, stirlerr, log1p, gam1, algdiv, series, logspaceAdd,
       lgamma, lbeta, gratio;

   C = require('../constants');
   stirlerr = require('./stirlerr').stirlerr;
   log1p = require('./log1p').log1p;
   gam1 = require('./gam1').gam1;
   lgamma = require('./lgamma').lgamma;
   lbeta = require('./lbeta').lbeta;
   series = require('../utils').series;
   logspaceAdd = require('../utils').logspaceAdd;
   gratio = require('./gratio').gratio;

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
      function algdiv(a, b) {
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
      }
   }());

   // bcorr
   // Computes stirlerr(a) + stirlerr(b) - stirlerr(a + b)
   function bcorr(a, b) {
      return a < b ?
             stirlerr(a) + algdiv(b, a + b) :
             stirlerr(b) + algdiv(a, a + b);
   }

   // Assumes lp is the log of the lower prob
   // Returns a pair lower/upper in log space
   // fromUpper assumes given upper prob in log space
   function fromLower(lp) {
      return { lower: lp, upper: log1p(-Math.exp(lp)) };
   }
   function fromLower(lq) {
      return { lower: log1p(-Math.exp(lq)), upper: lq };
   }

   // bpser (7) Returns the log of the lower probability
   function bpser(a, b, x) {
      var prob;

      return -lbeta(a, b) + a * Math.log(x) - Math.log(a) +
         log1p(a * series(function(i, v) {
            if (i === 0) { return 0; }
            if (i === 1) { return (1 - b) * x / (a + 1); }
            return v * (i - b) * x * (a + i - 1) / i / (a + i);
         }));
         // TODO: Check for need of stopping condition
   }

   // bup (8)
   // Assumes n >= 1
   function bup(a, b, x, y, n) {
      return a * Math.log(x) + b * Math.log(y) +
         Math.log(series(function(i) {
            return Math.pow(x, i) / Math.exp(lbeta(b, a + i)) / (a + i + 1);
         }, n));
      // TODO: Should it perhaps not return the log?
   }

   // bgrat (9)
   // log space
   // Assumes a > b
   // Not using the "grat_r" version as in r-code.
   function bgrat(a, b, x) {
      var T, fTs, u, logH, logM, lnx2s, lnx2sn, ps, cs, j;

      T = a + 0.5 * (b - 1);
      fTs = 1 / (4 * T * T);
      u = -T * Math.log(x);
      logH = Math.log(b) + log1p(gam1(b)) + b * Math.log(u) - u;
      logM = logH - (algdiv(b, a) + b * Math.log(T));
      cs = [1];   // cs[n] = 1 / (2n + 1)!
      ps = [1];   // ps[n] = p_n in (9.2)
      lnx2s = Math.log(x) * 0.5;   // (ln(x) / 2) ^ 2
      lnx2s = lnx2s * lnx2s;
      lnx2sn = 1;   // lnx2s ^ n
      j = gratio(b)(u) / Math.exp(logH);
      return logM + Math.log(
         series(function(n) {
            var mbmn = -n; // holds m*b - n
            if (n > 0) {
               cs[n] = cs[n - 1] / (2 * n) / (2 * n + 1);
               ps[n] = (b - 1) * cs[n] + 1 / n * series(function(m) {
                  if (m === 0) { return 0; }
                  mbmn = mbmn + b;
                  return mbmn * cs[m] * ps[n - m];
               }, n);
            }
            lnx2sn = lnx2sn * lnx2s;
            j = (b + 2 * n) * (b + 2 * n + 1) * j + (u + b + 2 * n + 1) * lnx2sn;
            j = j * fTs;
            return ps[i] * j;
         });
      );
   }

   function bratioLog(a, b, x) {
      var p, q, y, bint, bbar;

      p = a / (a + b);
      q = b / (a + b);
      y = 1 - x;
      bint = Math.floor(b);
      bbar = b - bint;

      // By convention, if 0 < x < 1 then a == 0 -> result = 1
      // and b == 0 -> result = 0. Can't have a == b == 0.
      if (a === 0) { return { lower: 0, upper: -Infinity }; }
      if (b === 0) { return { lower: -Infinity, upper: 0 }; }

      // R code handles specially the case where max(a, b) < eps * .001
      // but we will not (for now, anyway)

      // This logic follows pages 368-369 of algorithm 708
      /*  CASE where min(a, b) <= 1 */
      if (Math.min(a, b) <= 1) {
         if ( x > 0.5) { return bratioLog(b, a, y); }
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
            // TODO: bgrat here is supposed to take a w=bup(b, a, 1 - x, x, 20??)
            return fromUpper(bgrat(b + 20, a, y)); // (15a, 15b)
         }
         // max(a, b) <= 1
         if (a >= Math.min(0.2, b)) { return fromLower(bpser(a, b, x)); } // (12a)
         // a < min(0.2, b)
         if (Math.pow(x, a) <= 0.9) { return fromLower(bpser(a, b, x)); } // (12b)
         // x^a > 0.9
         if (x >= 0.3) { return fromUpper(bpser(b, a, y)); } // (13a)
         // x < 0.3
         // TODO: bgrat here is supposed to take a w=bup(b, a, 1 - x, x, 20??)
         return fromUpper(bgrat(b + 20, a, y));  // 15c
      }
      // min(a, b) > 1
      if (x > p) { return bratioLog(b, a, y); }
      if (b < 40) {
         if (b * x <= 0.7) { return fromLower(bpser(a, b, x)); } // (16a)
         // b * x > 0.7
         if (x <= 0.7) {
            return fromLower(logspaceAdd(bup(bbar, a, y, x, bint),
                                         bpser(a, bbar, x)));
         } // (17a)
         // x > 0.7
         // TODO need w0=bup(bbar, a, y, bint)
         if (a > 15) {
            return fromLower(logspaceAdd(bup(bbar, a, y, x, bint),
                                         bgrat(a, bbar, x)));
         } // (18a)
         // a <= 15
         // TODO: w0 =
         return fromLower(logspaceAdd(
                              logspaceAdd(bup(bbar, a, y, x, bint),
                                          bup(a, bbar, x, y, 20)),
                              bgrat(a + 20, bbar, x))); // (19a)
      }
      // b >= 40
      if (a <= b) {
         if (a <= 100) { return fromLower(bfrac(a, b, x)); } // (20a)
         // a > 100
         if (x < 0.97 * p) { return fromLower(bfrac(a, b, x)); } // (20b)
         // x >= 0.97 * p
         return fromLower(basym(a, b, x)); // (21a)
      }
      // a > b
      if (b <= 100) { return fromLower(bfrac(a, b, x)); } // (20c)
      // b > 100
      if (y > 1.03 * q) { return fromLower(bfrac(a, b, x)); } // (20d)
      // y <= 1.03 * q
      return fromLower(basym(a, b, x)); // (21c)
   }

   /*
    * bratio
    *
    * Return the tails of the incomplete beta function
    * of x, a, b.
    */
   function bratio(a, b, x) {
     var resLog;
     resLog = bratioLog(a, b, x);
     return {
        lower: Math.exp(resLog.lower),
        upper: Math.exp(resLog.upper)
     };
   }

    bratio.log = bratioLog;

   return {
      bratio: bratio
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
