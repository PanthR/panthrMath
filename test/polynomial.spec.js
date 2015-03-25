var Polynomial = require('../panthrMath/polynomial');
var chai = require('chai');
var expect = chai.expect;

describe('Polynomial class', function() {
   it('evalAt', function() {
      var arr = [3, 4, 0, 2, 1];
      var P = new Polynomial(arr);
      var v;
      function f(x) {
         return 3*x*x*x*x + 4*x*x*x + 2*x + 1;
      }
      for (var i = 0; i < 100; i += 1) {
         v = Math.random() * 200 - 100;
         expect(P.evalAt(v)).to.be.closeTo(f(v),.0001);
      }
   });
});
