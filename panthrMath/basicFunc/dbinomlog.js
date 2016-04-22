(function(define) {
'use strict';
define(function(require) {

   var C, stirlerr, bd0, log1p;

   C = require('../constants');
   stirlerr = require('./stirlerr').stirlerr;
   bd0 = require('./bd0').bd0;
   log1p = require('./log1p').log1p;

   /**
    * Performs the raw binomial log computation, even when $x$ is not an
    * integer.
    *
    * Based on:  *Fast and Accurate Computation of Binomial Probabilities*,
    * by Catherine Loader, 2000
    * @memberof basicFunc
    */
   function dbinomLog(size, p) {
      if (p === 0) { return function(x) { return x === 0 ? 0 : -Infinity; }; }
      if (p === 1) { return function(x) { return x === size ? 0 : -Infinity; }; }
      return function(x) {
         if (x === 0) {
            if (size === 0) { return 0; }
            return size * log1p(-p);
         }
         if (x === size ) { return size * Math.log(p); }
         return stirlerr(size) - stirlerr(x) - stirlerr(size - x) -
            bd0(x, size * p) - bd0(size - x, size * (1 - p)) +
            0.5 * Math.log(size / (C.twopi * x * (size - x)));
      };
   }

   return {
      dbinomLog: dbinomLog
   };

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
