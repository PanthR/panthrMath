(function(define) {'use strict';
define(function(require) {

   var C, bd0, stirlerr, twopi, sqrt2pi, Rational, pt, qt;

   C = require('../constants');
   bd0 = require('../basicFunc/bd0').bd0;
   stirlerr = require('../basicFunc/stirlerr').stirlerr;
   twopi = C.twopi;
   sqrt2pi = C.sqrt2pi;
   Rational = require('../rational');

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
   function dt(n) {
      var dtl;
      dtl = dtlog(n);
      return function(x) {
         return Math.exp(dtl(x));
      };
   }



   return {
      dt: dt,
      dtlog: dtlog,
      pt: pt,
      qt: qt
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
