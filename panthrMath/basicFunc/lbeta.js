(function(define) {'use strict';
define(function(require) {

   /*
    * log beta computation
    * Based on R's code
    */
    /**
     * TODO
     * @memberof basicFunc
     */
   var stirlerr, sqrt2pi, log1p, lgamma;

   stirlerr = require('./stirlerr').stirlerr;
   log1p = require('./log1p').log1p;
   lgamma = require('./lgamma').lgamma;
   sqrt2pi = require('../constants').sqrt2pi;

   /**
    * TODO
    * @memberof basicFunc
    */
   function lbeta(a, b) {
      if (a > b) { return lbeta(b, a); }
      // a <= b
      if (a < 0) { return NaN; }
      if (a === 0) { return Infinity; }
      if (b === Infinity) { return -Infinity; }
      if (a > 10) {
         return -0.5 * Math.log(b) + Math.log(sqrt2pi) +
                (stirlerr(a) + stirlerr(b) - stirlerr(a + b)) +
                (a - 0.5) * Math.log(a / (a + b)) +
                b * log1p(-a / (a + b));
      }
      if (b >= 10) {
         return lgamma(a) + (stirlerr(b) - stirlerr(a + b)) +
                a - a * Math.log(a + b) +
                (b - 0.5) * log1p(-a / (a + b));
      }
      return lgamma(a) + (lgamma(b) - lgamma(a + b));
   }

   function beta(a, b) {
      return Math.exp(lbeta(a, b));
   }

   return {
      lbeta: lbeta,
      beta: beta
   };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
