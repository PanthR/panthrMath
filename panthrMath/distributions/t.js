(function(define) {'use strict';
define(function(require) {

   var C, bd0, stirlerr, sqrt2pi, bratio, log1p, solve;
   // var expm1,  qnorm;

   C = require('../constants');
   bd0 = require('../basicFunc/bd0').bd0;
   stirlerr = require('../basicFunc/stirlerr').stirlerr;
   sqrt2pi = C.sqrt2pi;
   bratio = require('../basicFunc/bratio').bratio;
   log1p = require('../basicFunc/log1p').log1p;
   // expm1 = require('../basicFunc/expm1').expm1;
   solve = require('../utils').binSearchSolve;
   // qnorm = require('./normal').qnorm;

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
         var val;

         val = df > x * x ?
               bratio(0.5, df / 2, x * x / (df + x * x), false, true) :
               bratio(df / 2, 0.5, 1 / (1 + x / df * x), true, true);

         if (x <= 0) { lowerTail = !lowerTail; }

         return lowerTail ? log1p(-0.5 * Math.exp(val))
                          : val - Math.log(2);
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

   // inverse cdf
   // From qt.c in R
   function qt(df, lowerTail, logp) {
      lowerTail = lowerTail !== false;
      logp = logp === true;

      return function(p) {
         // var pp;

         // two-tailed prob, pp = 2 * min(p, 1-p)
         // pp = 2 * (logp ? Math.min(Math.exp(p), -expm1(p))
                        // : Math.min(p, 1 - p));

         // Using the solver on whole range (possibly inefficient but works)
         // if (df < 1) {
            return solve(function(x) {
               return pt(df, lowerTail, logp)(x);
            }, p);
         // }

         // if (df > 1e20) {
         //    return qnorm(0, 1, lowerTail, logp)(p);
         // }

         // TODO: Could consider optimizing for df close to 1 or 2

      };
   }

   return {
      tdistr: function(df) {
         return {
            d: function(x, logp) { return dt(df, logp)(x); },
            p: function(q, lowerTail, logp) {
               return pt(df, lowerTail, logp)(q);
            },
            q: function(p, lowerTail, logp) {
               return qt(df, lowerTail, logp)(p);
            },
            r: function(n) { return rt(df)(n); }
         };
      },
      dt: dt,
      pt: pt,
      qt: qt
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
