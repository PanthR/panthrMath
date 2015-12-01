(function(define) {'use strict';
define(function(require) {

   // Poisson distribution
   // No input validation provided.

   var lpoisson, pgamma, pWrap, qnorm, discInvCdf;

   lpoisson = require('../basicFunc/lpoisson').lpoisson;
   pgamma = require('./gamma').pgamma;
   pWrap = require('../utils').pWrap;
   discInvCdf = require('../utils').discInvCdf;
   qnorm = require('./normal').qnorm;

   // density / pdf
   function dpois(lambda, logp) {
      logp = logp === true;
      return function(x) {
         return logp ? lpoisson(lambda)(x) : Math.exp(lpoisson(lambda)(x));
      };
   }

   function ppois(lambda, lowerTail, logp) {
      lowerTail = lowerTail !== false;
      logp = logp === true;

      if (!(lambda >= 0)) { return function(x) { return NaN; }; }

      return function(x) {
         var ret;

         if (x >= 0 && lambda > 0) {
            return pgamma(Math.floor(x + 1e-10) + 1, 1, !lowerTail, logp)(lambda);
         }

         ret = lowerTail ? 1 : 0;
         if (x < 0) { ret = lowerTail ? 0 : 1; }
         return logp ? Math.log(ret) : ret;
      };
   }

   function qpois(lambda, lowerTail, logp) {

      // TODO:  make a separate calculation for !lowerTail

      var goodParams, mu, sigma, gamma;

      logp = logp === true;
      lowerTail = lowerTail !== false;
      goodParams = lambda >= 0 && lambda < Infinity;

      if (!goodParams) { return function(prob) { return NaN; }; }
      if (lambda === 0) {
         return pWrap(lowerTail, logp, function(prob) { return 0; });
      }

      mu = lambda;
      sigma = Math.sqrt(lambda);
      gamma = 1 / sigma;

      return pWrap(true, logp, function(prob) {
         var z, ret;
         z = qnorm(0, 1, lowerTail)(prob); // initial value
         if (z < -10) { z = -10; }
         if (z > 10) { z = 10; }
         ret = Math.floor(mu + sigma * (z + gamma * (z * z - 1) / 6) + 0.5);
         if (prob === 1) { return lowerTail ? Infinity : 0; }
         if (lowerTail) {
            return discInvCdf(0, 1e10, ret, prob, ppois(lambda));
         }
         return -discInvCdf(-1e10, 0, -ret, prob, function(x) {
            return ppois(lambda, lowerTail)(-x);
         });
      });
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
      dpois: dpois,
      ppois: ppois,
      qpois: qpois
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
