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
   var C, stirlerr, lgamma, algdiv, series;

   C = require('../constants');
   stirlerr = require('../basicfunc/stirlerr').stirlerr;
   lgamma = require('../basicfunc/lgamma').lgamma;
   series = require('../utils').series;

   // Page 371 of DiDonato/Morris
   // Computes lgamma(b) - lgamma(a + b)
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

   // bpser (7)
   function bpser(a, b, x) {

   }

   function bratioLog(a, b, x) {
      var p, q, y, bint, bfrac;

      p = a / (a + b);
      q = b / (a + b);
      y = 1 - x;
      bint = Math.floor(b);
      bfrac = b - bint;

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
            if (b <= 1) { return bpser(a, b, x); } // (12c)
            // a <= 1 < b
            if (x >= 0.3) { return tailFlip(bpser(b, a, y)); } // (13b)
            // x < 0.3
            if (x < 0.1 && Math.pow(x * b, a) <= 0.7) { return bpser(a, b, x); } // (12d)
            // x >= 0.1 or (bx)^a > 0.7
            // for now, assume bgrat will start with w0 = 0
            if (b > 15) { return tailFlip(bgrat(b, a, y)); } // (14a, 14b)
            // b <= 15
            // TODO: bgrat here is supposed to take a w=bup(b, a, 1 - x, 20??)
            return tailFlip(bgrat(b + 20, a, y)); // (15a, 15b)
         }
         // max(a, b) <= 1
         if (a >= Math.min(0.2, b)) { return bpser(a, b, x); } // (12a)
         // a < min(0.2, b)
         if (Math.pow(x, a) <= 0.9) { return bpser(a, b, x); } // (12b)
         // x^a > 0.9
         if (x >= 0.3) { return tailFlip(bpser(b, a, y)); } // (13a)
         // x < 0.3
         // TODO: bgrat here is supposed to take a w=bup(b, a, 1 - x, 20??)
         return tailFlip(bgrat(b + 20, a, y));  // 15c
      }
      // min(a, b) > 1
      if (x > p) { return bratioLog(b, a, y); }
      if (b < 40) {
         if (b * x <= 0.7) { return bpser(a, b, x); } // (16a)
         // b * x > 0.7
         if (x <= 0.7) { return bup(bfrac, a, y, bint) + bpser(a, bfrac, x); } // (17a)
         // x > 0.7
         // TODO need w0=bup(bfrac, a, y, bint)
         if (a > 15) { return bgrat(a, bfrac, x); } // (18a)
         // a <= 15
         // TODO: w0 = bup(bfrac, a, y, bint) + bup(a, bfrac, x, 20)
         return bgrat(a + 20, bfrac, x); // (19a)
      }
      // b >= 40
      if (a <= b) {
         if (a <= 100) { return bfrac(a, b, x); } // (20a)
         // a > 100
         if (x < 0.97 * p) { return bfrac(a, b, x); } // (20b)
         // x >= 0.97 * p
         return basym(a, b, x); // (21a)
      }
      // a > b
      if (b <= 100) { return bfrac(a, b, x); } // (20c)
      // b > 100
      if (y > 1.03 * q) { return bfrac(a, b, x); } // (20d)
      // y <= 1.03 * q
      return basym(a, b, x); // (21c)
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
