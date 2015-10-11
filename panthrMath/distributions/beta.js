(function(define) {'use strict';
define(function(require) {

   var bratio, dbeta, solve;

   bratio = require('../basicFunc/bratio').bratio;
   solve = require('../utils').binSearchSolve;

   function pbeta(a, b, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      return function(x) {
         var lp;
         if (a <= 0 || b <= 0) { return NaN; }
         if (x > 0 && x < 1) {
            return bratio(a, b, x, lowerTail, logp);
         }
         lp = x <= 0 ? 0 : 1;
         lp = lowerTail ? lp : 1 - lp;
         return logp ? Math.log(lp) : lp;
      };
   }

   // inverse cdf
   // preliminary implementation uses binSearchSolve, similar to our qt
   function qbeta(a, b, lowerTail, logp) {
      var f;

      lowerTail = lowerTail !== false;
      logp = logp === true;
      f = pbeta(a, b, lowerTail, logp);

      /* eslint-disable complexity */
      return function(p) {
         // need reasonable endpoints
         var l, u, incr, realLogP;

         if (a <= 0 || b <= 0) { return NaN; }

         if (!logp && !(p >= 0 && p <= 1)) { return NaN; }
         realLogP = logp ? p : Math.log(p);
         realLogP = lowerTail ? realLogP : -realLogP;

         if (realLogP === -Infinity) { return 0; }
         if (realLogP === +Infinity) { return 1; }

         incr = 1e-1;
         l = 1e-1;
         while (f(l) > p === lowerTail) {
            l = l * incr;
            if (l < 1e-80) {
               return new Error('qbeta failed to find left endpoint');
            }
         }
         u = 1 - incr;
         while (f(u) < p === lowerTail) {
            incr *= .1;
            u = 1 - incr;
            if (incr < 1e-80) {
               return new Error('qbeta failed to find right endpoint');
            }
         }
         /* eslint-enable complexity */

         return solve(function(x) { return f(x); }, p, l, u);
      };
   }

   return {
      dbeta: dbeta,
      pbeta: pbeta,
      qbeta: qbeta
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
