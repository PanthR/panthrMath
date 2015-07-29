(function(define) {'use strict';
define(function(require) {

   var basicfunc, C;

   basicfunc = require('../basicfunc');
   C = require('../constants');

   /*
    * bratioLog -- see Algorithm 708
    *
    * Return the logs of the tails of the incomplete beta function
    * of x, a, b.
    */
   function bratioLog(a, b, x) {

      // By convention, if 0 < x < 1 then a == 0 -> result = 1
      // and b == 0 -> result = 0. Can't have a == b == 0.
      if (a === 0) { return { lower: 0, upper: -Infinity }; }
      if (b === 0) { return { lower: -Infinity, upper: 0 }; }

      // R code handles specially the case where max(a, b) < eps * .001
      // but we will not (for now, anyway)

      /*  CASE where min(a, b) <= 1 */
      if (Math.min(a, b) <= 1) {
         if ( x > 0.5) { return bratioLog(b, a, 1 - x); }
         // we are skipping the case, in R code, for b < min(eps, a*eps)
         // also skipping a < min(eps, b*eps)
         if ( Math.max(a, b) > 1) {
            if (b <= 1) { return bpser(a, b, x); }
            if (x >= 0.29) { return tailFlip(bpser(b, a, 1 - x)); }
            if (x < 0.1 && Math.pow(x * b, a) <= 0.7) { return bpser(a, b, x); }
            // for now, assume bgrat will start with w0 = 0
            if (b > 15) { return bgrat(b, a, 1 - x); }
         }
      }

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



   return {
      bratio: bratio
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
