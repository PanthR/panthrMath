(function(define) {
'use strict';
define(function(require) {

   /**
    * Computes the expression $$bd0(x, y) = x \ln(x/y) + y - x$$
    *
    * Based on:  *Fast and Accurate Computation of Binomial Probabilities*,
    * by Catherine Loader, 2000
    * @memberof basicFunc
    */
   function bd0(x, y) {
      var ej, j, s, sTemp, v;

      if (Math.abs(x - y) < 0.1 * (x + y)) {
         v = (x - y) / (x + y);
         s = (x - y) * v;
         ej = 2 * x * v;
         for (j = 1; j < 10000; j += 1) {
            ej *= v * v;
            sTemp = s + ej / (2 * j + 1);
            if (sTemp === s) {
               return sTemp;
            }
            s = sTemp;
         }
         throw new Error('bd0(x, y) made it to 10000 ...');
      }
      return x * Math.log(x / y) + y - x;
   }

   return { bd0: bd0 };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
