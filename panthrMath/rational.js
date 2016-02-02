(function(define) {
'use strict';
define(function(require) {

   var Polynomial;

   Polynomial = require('./polynomial');

   /**
    * A small wrapper class for representing rational functions.
    * @memberof panthrMath
    * @author Haris Skiadas <skiadas@hanover.edu>, Barb Wahl <wahl@hanover.edu>
    */
    function Rational(num, denom) {
      this.num = new Polynomial(num);
      this.denom = new Polynomial(denom);
   }
   /**
    * Creates a rational function.
    *
    * `num` and `denom` can be instances of `Polynomial`
    * or arrays of coefficients for `Polynomial.new`.
    *
    * Common factors between the numerator and denominator are *not*
    * eliminated.
    * @memberof panthrMath
    */
   Rational.new = function(num, denom) {
      return new Rational(num, denom);
   };
   /**
    * Evaluates the rational function at the value `x`.
    * @memberof panthrMath
    */
   Rational.prototype.evalAt = function evalAt(x) {
      return this.num.evalAt(x) / this.denom.evalAt(x);
   };

   return Rational;

});

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
