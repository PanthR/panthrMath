(function(define) {'use strict';
define(function(require) {

   /**
    * Provides density function, cdf, inverse cdf, and random
    * number generator for the student's t distribution.
    *
    * @module distributions.t
    * @memberof distributions
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   var C, bd0, stirlerr, sqrt2pi, bratio, log1p, solve, rgen;
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
   rgen = require('../rgen/rgen');

   /*
    * Return the log-of-density function for student's t distribution.
    * `df` is the degrees of freedom.
    * See: R source code.
    */
   function dtlog(df) {
      var t;
      t = -bd0(df / 2, (df + 1) / 2) + stirlerr((df + 1) / 2) - stirlerr(df / 2);
      return function(x) {
         var x2n, logx2n, u;
         x2n = x * x / df;
         logx2n = Math.log(1 + x2n) / 2;
         u = x2n > 0.2 ? df * logx2n
                       : x * x / 2 - bd0(df / 2, (df + x * x) / 2);
         return t - u - (logx2n + Math.log(sqrt2pi));
      };
   }

   // density function
   // dt(df) returns a function for calculating cumulative distribution
   /**
    * TODO
    * @fullName dt(df, logp)(x)
    * @memberof t
    */
   function dt(df, logp) {
      logp = logp === true;
      var dtl;
      dtl = dtlog(df);
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
   /**
    * TODO
    * @fullName pt(df, lowerTail, logp)(x)
    * @memberof t
    */
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
   /**
    * TODO
    * @fullName qt(df, lowerTail, logp)(p)
    * @memberof t
    */
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

   /**
    * TODO
    * @memberof t
    */
   function rt(df) {
         return function() {
            var u1, u2, rsq;
            do {
               u1 = rgen.random() * 2 - 1;
               u2 = rgen.random() * 2 - 1;
               rsq = u1 * u1 + u2 * u2;
            } while (rsq >= 1);
            return u1 * Math.sqrt(df * ( Math.pow(rsq, -2 / df) - 1) / rsq);
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
      qt: qt,
      rt: rt
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
