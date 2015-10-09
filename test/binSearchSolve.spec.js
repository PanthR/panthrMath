var chai = require('chai');
var expect = chai.expect;
var precision = 1e-10;
var utils = require('../panthrMath/utils');
var solve = utils.binSearchSolve;

describe('binSearchSolve function', function() {
   var f = function(x) { return Math.exp(x); };
   it('works when interval is given', function() {
      var i, a, b, L, c;
      for (i = 0; i < 100; i += 1) {
         a = Math.random() * 100 - 100;
         b = Math.random() * 100 + a;
         L = Math.random() * (f(b) - f(a)) + f(a);
         c = solve(f, L, a, b);
         expect(utils.relativelyCloseTo(f(c), L)).to.be.ok;
      }
   });
   it('Can determine interval on its own', function() {
      var i, a, b, L, c;
      for (i = 0; i < 100; i += 1) {
         a = Math.random() * 100 - 100;
         b = Math.random() * 100 + a;
         L = Math.random() * (f(b) - f(a)) + f(a);
         c = solve(f, L);
         expect(utils.relativelyCloseTo(f(c), L)).to.be.ok;
      }
   });
   it('Works for a decreasing function', function() {
      var i, a, b, L, c, g;
      g = function(x) { return Math.exp(-x); }
      for (i = 0; i < 100; i += 1) {
         a = Math.random() * 100 - 100;
         b = Math.random() * 100 + a;
         L = Math.random() * (g(b) - g(a)) + g(a);
         c = solve(g, L);
         expect(utils.relativelyCloseTo(g(c), L)).to.be.ok;
      }
   });

});
