(function(define) {'use strict';
define(function(require) {

   // creates a polynomial function from the array of coefficients
   // for example, if the function is ax^3 + bx^2 + cx + d
   // the array is [a, b, c, d], that is, the array indexing is
   // the reverse of the usual coefficient indexing

   /**
    * Polynomial Class
    * @memberof panthrMath
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   function Polynomial(coefs) {
      if (coefs instanceof Polynomial) { return coefs; }
      this.coefs = coefs;
   }
   /**
    * TODO
    * @memberof panthrMath
    */
   Polynomial.new = function(coefs) {
      return new Polynomial(coefs);
   };
   /**
    * TODO
    * @memberof panthrMath
    */
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
