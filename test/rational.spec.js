var Rational = require('../panthrMath/rational');
var Polynomial = require('../panthrMath/polynomial');
var chai = require('chai');
var expect = chai.expect;

describe('Rational class', function() {
   it('has an evalAt', function() {
      var P = Polynomial.new([1, 2, 3]);
      var Q = Polynomial.new([2, -1, 2, 4]);
      var R = Rational.new(P, Q);
      function f(x) {
         return P.evalAt(x) / Q.evalAt(x);
      }
      for (var i = 0; i < 100; i += 1) {
         v = Math.random() * 200 - 100;
         expect(R.evalAt(v)).to.be.closeTo(f(v),.0001);
      }
   });
});
