(function(define) {'use strict';
define(function(require) {

   /**
    * Provides density function, cumulative distribution function,
    * quantile function, and random number generator
    * for the Gamma distribution.
    *
    * @module distributions.gamma
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var lpoisson, gratio, exponential, rgen;

   lpoisson = require('../basicFunc/lpoisson').lpoisson;
   gratio = require('../basicFunc/gratio');
   exponential = require('../rgen/exponential');
   rgen = require('../rgen/rgen');


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
   /**
    * TODO
    * @memberof gamma
    */
   function dgamma(a, s, logp) {
      logp = logp === true;
      return function(x) {
         return logp ? dgammaLog(a, s)(x) : Math.exp(dgammaLog(a, s)(x));
      };
   }

   // cumulative distribution function
   // a = shape, s = scale (rate = 1 / s)
   /**
    * TODO
    * @memberof gamma
    */
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
   /**
    * TODO
    * @memberof gamma
    */
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

   // rgamma stuff...
   // Client calls one of three methods depending on the value of a.
   // a for 'shape', s for 'scale'
   /**
    * TODO
    * @memberof gamma
    */
   function rgamma(a, s) {
      if (a < 1) {
         return gammaBest(a, s);
      }
      if (a > 1) {
         return gammaCheng(a, s);
      }
      return exponential(1 / s);
   }

   // a < 1:  Best/Ahrens/Dieter Algorithm
   function gammaBest(a, s) {
      var t, b;
      t = .07 + .75 * Math.sqrt(1 - a);
      b = 1 + Math.exp(-t) * a / t;

      return function() {
         var x, y, u1, u2, v;
         /* eslint-disable no-constant-condition */
         while (true) {
         /* eslint-enable */
            u1 = rgen.random();
            u2 = rgen.random();
            v = b * u1;
            if (v <= 1) {
               x = t * Math.pow(v, 1 / a);
               if (u2 <= (2 - x) / (2 + x) || u2 <= Math.exp(-x)) {
                  return x * s;
               }
            } else {
               x = -Math.log(t * (b - v) / a);
               y = x / t;
               if (u2 * (a + y * (1 - a)) <= 1 ||
                   u2 <= Math.pow(y, a - 1))
               {
                  return x * s;
               }
            }
         }
      };
   }

   // a > 1:  Cheng/Feast Algorithm
   function gammaCheng(a, s) {
      var c1, c2, c3, c4;
      c1 = a - 1;
      c2 = (a - 1 / (6 * a)) / c1;
      c3 = 2 / c1;
      c4 = c3 + 2;

      return function() {
         var u1, u2, w;
         /* eslint-disable no-constant-condition */
         while (true) {
         /* eslint-enable */
            u1 = rgen.random();
            u2 = rgen.random();
            w = c2 * u2 / u1;
            if (c3 * u1 + w + 1 / w <= c4) {
               return c1 * w * s;
            }
            if (c3 * Math.log(u1) - Math.log(w) + w < 1) {
               return c1 * w * s;
            }
         }
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
      qgamma: qgamma,
      rgamma: rgamma
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
