(function(define) {'use strict';
define(function(require) {

   var C, bd0, stirlerr;

   C = require('../constants');
   bd0 = require('./bd0').bd0;
   stirlerr = require('./stirlerr').stirlerr;

   // returns the log of the poisson distribution
   // Based on dpois from Loader (2000).

   /**
    * Computes the logarithm of the Poisson pdf
    * $$\log\left(\frac{\lambda^x}{\Gamma(x+1)}e^{-\lambda}\right)
    * = -\lambda + x \ln(\lambda) - \ln(\Gamma(x+1))$$
    * where $x>0$ and $\lambda>0$.
    *
    * Based on:  *Fast and Accurate Computation of Binomial Probabilities*,
    * by Catherine Loader, 2000
    * @fullName lpoisson(lambda)(x)
    * @memberof basicFunc
    */
   function lpoisson(lambda) {
      if (lambda === 0) {
         return function(x) {
            return x === 0 ? 0 : -Infinity;
         };
      }
      return function(x) {
         if (x === 0) { return -lambda; }
         return -stirlerr(x) - bd0(x, lambda) - 0.5 * Math.log(C.twopi * x);
      };
   }

   return {
      lpoisson: lpoisson
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
