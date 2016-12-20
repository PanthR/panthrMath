(function(define) {
'use strict';
define(function(require) {

   var utils, lbeta, lgamma, kSmallMax;

   utils = require('../utils');
   lbeta = require('./lbeta').lbeta;
   lgamma = require('./lgamma').lgamma;

   kSmallMax = 30;

   /**
    * Computes the logarithm of the absolute value of the binomial coefficient.
    * Valid for any real n and for integer k.
    * k will be rounded if it is not an integer.
    */
   function lchoose(n, k) {
      k = Math.round(k);
      if (utils.hasNaN(n, k)) { return NaN; }
      if (k < 0) { return -Infinity; }
      if (k === 0) { return 0; }
      if (k === 1) { return Math.log(Math.abs(n)); }
      if (n < 0) { return lchoose(-n + k - 1, k); }
      if (n === Math.round(n)) {
         if (n < k) { return -Infinity; }
         if (n - k < 2) { return lchoose(n, n - k); }

         return lfastchoose(n, k);
      }
      if (n < k - 1) { return lfastchoose2(n, k); }

      return lfastchoose(n, k);
   }

   /**
    * Computes the binomial coefficient "n choose k".
    * Valid for any real n and for integer k.
    * k will be rounded if it is not an integer.
    */
   function choose(n, k) {
      var ret, j;

      k = Math.round(k);
      if (utils.hasNaN(n, k)) { return NaN; }
      if (k < kSmallMax) {
         if (n - k < k && n >= 0 && n === Math.round(n)) { k = n - k; }
         if (k < 0) { return 0; }
         if (k === 0) { return 1; }
         ret = n;
         for (j = 2; j <= k; j += 1) {
            ret *= (n - j + 1) / j;
         }

         return n === Math.round(n) ? Math.round(ret) : ret;
      }
      if (n < 0) {
         ret = choose(-n + k - 1, k);

         return k % 2 === 1 ? -ret : ret;
      }
      if (n === Math.round(n)) {
         if (n < k) { return 0; }
         if (n - k < kSmallMax) { return choose(n, n - k); }

         return Math.round(Math.exp(lfastchoose(n, k)));
      }
      if (n < k - 1) {
         return signGamma(n - k + 1) * Math.exp(lfastchoose2(n, k));
      }

      return Math.exp(lfastchoose(n, k));
   }

   // Helper functions
   function lfastchoose(n, k) {
      return -Math.log(n + 1) - lbeta(n - k + 1, k + 1);
   }

   function lfastchoose2(n, k) {
      return lgamma(n + 1) - lgamma(k + 1) - lgamma(-(n - k + 1));
   }

   // Computes the sign of gamma(x). Based on R's lgammafn_sgn
   function signGamma(x) {
      return x < 0 && Math.floor(-x) % 2 === 0 ? -1 : 1;
   }

   return {
      choose: choose,
      lchoose: lchoose
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
