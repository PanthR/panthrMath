(function(define) {'use strict';
define(function(require) {

   // Poisson distribution
   // No input validation provided.

   var C, bd0, stirlerr;

   C = require('../constants');
   bd0 = require('../basicFunc/bd0').bd0;
   stirlerr = require('../basicFunc/stirlerr').stirlerr;

   // returns the log of the poisson distribution
   // Based on dpois from Loader (2000).
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

   // density / pdf
   function dpois(lambda) {
      return function(x) {
         return Math.exp(lpoisson(lambda)(x));
      };
   }

   dpois.log = lpoisson;

   return {
      dpois: dpois
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
