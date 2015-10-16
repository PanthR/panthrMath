(function(define) {'use strict';
define(function(require) {

   function dunif(min, max, logp) {
      logp = logp === true;

      if (!(min < max)) { return function(x) { return NaN; }; }

      return function(x) {
         var p;
         p = x < min || x > max ? 0 : 1 / (max - min);

         return logp ? Math.log(p) : p;
      };
   }

   function punif(min, max, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (!(min < max)) { return function(x) { return NaN; }; }

      return function(x) {
         var p;

         p = x <= min ? 0 :
             x >= max ? 1 :
                        (x - min) / (max - min);
         p = lowerTail ? p : 1 - p;

         return logp ? Math.log(p) : p;
      };
   }

   function qunif(min, max, lowerTail, logp) {
      logp = logp === true;
      lowerTail = lowerTail !== false;

      if (!(min < max)) { return function(x) { return NaN; }; }

      return function(p) {
         p = logp ? Math.exp(p) : p;

         if (p < 0 || p > 1) { return NaN; }
         return lowerTail ? min + p * (max - min) : max - p * (max - min);
      };
   }

   return {
      dunif: dunif,
      punif: punif,
      qunif: qunif
   };
});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
