(function(define) {
'use strict';
define(function(require) {

   /**
    * A small wrapper class for representing polynomial functions.
    * @memberof panthrMath
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
   function Polynomial(coefs) {
      if (coefs instanceof Polynomial) { return coefs; }
      this.coefs = coefs;
   }
   /**
    * Creates a polynomial function from `coefs`, an array of coefficients.
    * For example, if the function is $ax^3 + bx^2 + cx + d$
    * the array is `[a, b, c, d]`, that is, the array indexing is
    * the reverse of the usual coefficient indexing.
    * @memberof panthrMath
    */
   Polynomial.new = function(coefs) {
      return new Polynomial(coefs);
   };
   /**
    * Evaluates the polynomial at a given value `x`.
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
