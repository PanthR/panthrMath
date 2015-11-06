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
   function dpois(lambda, logp) {
      logp = logp === true;
      return function(x) {
         return logp ? lpoisson(lambda)(x) : Math.exp(lpoisson(lambda)(x));
      };
   }

   return {
      pois: function(lambda) {
         return {
            d: function(x, logp) { return dpois(lambda, logp)(x); },
            p: function(q, lowerTail, logp) {
               return ppois(lambda, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qpois(lambda, lowerTail, logp)(p);
            },
            r: function(n) { return rpois(lambda)(n); }
         };
      },
      dpois: dpois
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
