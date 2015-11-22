(function(define) {'use strict';
define(function(require) {

   // Gamma distribution

   var lpoisson, gratio;

   lpoisson = require('../basicFunc/lpoisson').lpoisson;
   gratio = require('../basicFunc/gratio');

   // helper function
   function dgammaLog(a, s) {
      if (a < 1) {
         return function(x) {
            if (x === 0) { return Infinity; }
            return Math.log(a / x) + lpoisson(x / s)(a);
         };
      }
      if (a === 1) {
         return function(x) {
            return -x / s - Math.log(s);
         };
      }
      // a > 1
      return function(x) {
         if (x === 0) { return -Infinity; }
         return -Math.log(s) + lpoisson(x / s)(a - 1);
      };
   }

   // density function
   // 1 / (s^a * gamma(a)) * x ^ (a - 1) * exp(-x / s)
   function dgamma(a, s, logp) {
      logp = logp === true;
      return function(x) {
         return logp ? dgammaLog(a, s)(x) : Math.exp(dgammaLog(a, s)(x));
      };
   }

   // cumulative distribution function
   // a = shape, s = scale (rate = 1 / s)
   function pgamma(a, s, lowerTail, logp) {
      var f;

      if (s <= 0 || a <= 0) { return function(x) { return NaN; }; }
      logp = logp === true;
      lowerTail = lowerTail !== false;
      f = lowerTail ? gratio.gratio(a) : gratio.gratioc(a);

      return function(x) {
         var p;

         p = x > 0 ? f(x / s) : lowerTail ? 0 : 1;
         return logp ? Math.log(p) : p;
      };
   }

   // inverse cdf
   // a = shape, s = scale (rate = 1 / s)
   function qgamma(a, s, lowerTail, logp) {
      var gaminv;

      if (s <= 0 || a <= 0) { return function(x) { return NaN; }; }
      logp = logp === true;
      lowerTail = lowerTail !== false;
      gaminv = gratio.gaminv(a); // function of p, lowerTail

      return function(p) {
         p = logp ? Math.exp(p) : p;
         if (!(p >= 0 && p <= 1)) { return NaN; }
         return s * gaminv(p, lowerTail);
      };
   }

   return {
      gammadistr: function(a, s) {
         return {
            d: function(x, logp) { return dgamma(a, s, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pgamma(a, s, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qgamma(a, s, lowerTail, logp)(p);
            },
            r: function(n) { return rgamma(a, s)(n); }
         };
      },
      dgamma: dgamma,
      pgamma: pgamma,
      qgamma: qgamma
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
