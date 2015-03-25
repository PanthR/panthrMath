(function(define) {'use strict';
define(function(require) {

   // creates a polynomial function from the array of coefficients
   // for example, if the function is ax^3 + bx^2 + cx + d
   // the array is [a, b, c, d], that is, the array indexing is
   // the reverse of the usual coefficient indexing
   function Polynomial(coefs) {
      this.coefs = coefs;
   }

   Polynomial.prototype.evalAt = function evalAt(x) {
      return this.coefs.reduce(function(acc, coef) {
         return acc * x + coef;
      }, 0);
   };


   return Polynomial;

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
