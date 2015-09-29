(function(define) {'use strict';
define(function(require) {

   var C, bd0, stirlerr, sqrt2pi, qt, bratio;

   C = require('../constants');
   bd0 = require('../basicFunc/bd0').bd0;
   stirlerr = require('../basicFunc/stirlerr').stirlerr;
   sqrt2pi = C.sqrt2pi;
   bratio = require('../basicFunc/bratio').bratio;

   /*
    * Return the log-of-density function for student's t distribution.
    * `n` is the degrees of freedom.
    * See: R source code.
    */
   function dtlog(n) {
      var t;
      t = -bd0(n / 2, (n + 1) / 2) + stirlerr((n + 1) / 2) - stirlerr(n / 2);
      return function(x) {
         var x2n, logx2n, u;
         x2n = x * x / n;
         logx2n = Math.log(1 + x2n) / 2;
         u = x2n > 0.2 ? n * logx2n
                       : x * x / 2 - bd0(n / 2, (n + x * x) / 2);
         return t - u - (logx2n + Math.log(sqrt2pi));
      };
   }

   // density function
   // dt(n) returns a function for calculating cumulative distribution
   function dt(n, logp) {
      logp = logp === true;
      var dtl;
      dtl = dtlog(n);
      return function(x) {
         return logp ? dtl(x) : Math.exp(dtl(x));
      };
   }

   // cumulative distribution function, log version
   // based on pt.c from R implementation
   function ptlog(df, lowerTail) {
      return function(x) {
         return Math.log(0.5) +
            df > x * x ?
               bratio(0.5, df / 2, x * x / (df + x * x), !lowerTail, true) :
               bratio(df / 2, 0.5, 1 / (1 + x / df * x), lowerTail, true);
      };
   }

   // cumulative distribution function
   function pt(df, lowerTail, logp) {
      var ptl;

      lowerTail = lowerTail !== false;
      logp = logp === true;
      ptl = ptlog(df, lowerTail);
      return function(x) {
         return logp ? ptl(x) : Math.exp(ptl(x));
      };
   }

   return {
      dt: dt,
      pt: pt,
      qt: qt
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
