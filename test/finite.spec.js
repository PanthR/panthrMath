var chai = require('chai');
var expect = chai.expect;
var precision = 1e-10;
var utils = require('../panthrMath/utils');
var main = require('..');

describe('finite distribution', function() {
   it('specified via list of values and probabilities', function() {
      var d;
      d = main.finite({ xs: [2, 5, 7], ws: [1, 2, 1] });
      ['d', 'p', 'q', 'r'].forEach(function(s) {
         expect(d).to.respondTo(s);
      });
      [[1, 0], [2, 0.25], [5, 0.5], [6, 0], [7, 0.25], [8, 0]].forEach(
         function(pair) {
         expect(d.d(pair[0])).to.equal(pair[1]);
         expect(d.d(pair[0], true)).to.equal(Math.log(pair[1]));
      });
      [[1, 0], [2, 0.25], [3, 0.25], [5, 0.75],
       [6, 0.75], [7, 1], [8, 1]].forEach(function(pair) {
         expect(d.p(pair[0])).to.equal(pair[1]);
         expect(d.p(pair[0], true, true)).to.equal(Math.log(pair[1]));
         expect(d.p(pair[0], false)).to.equal(1 - pair[1]);
         expect(d.p(pair[0], false, true)).to.equal(Math.log(1 - pair[1]));
      });
      [[0, 2], [0.2, 2], [0.25, 2], [0.3, 5], [0.45, 5], [0.7, 5],
       [0.75, 5], [0.8, 7], [1, 7]].forEach(function(pair) {
         expect(d.q(pair[0])).to.equal(pair[1]);
         expect(d.q(Math.log(pair[0]), true, true)).to.equal(pair[1]);
         expect(d.q(1-pair[0], false)).to.equal(pair[1]);
         expect(d.q(Math.log(1-pair[0]), false, true)).to.equal(pair[1]);
      });
   });
   it('understands non-integer xs', function() {
      var d;
      d = main.finite({ xs: [-0.2, 5.1, 7], ws: [1, 2, 1] });
      ['d', 'p', 'q', 'r'].forEach(function(s) {
         expect(d).to.respondTo(s);
      });
      [[-1, 0], [-0.2, 0.25], [5.1, 0.5], [6, 0], [7, 0.25], [8, 0]].forEach(function(p) {
         expect(d.d(p[0])).to.equal(p[1]);
         expect(d.d(p[0], true)).to.equal(Math.log(p[1]));
      });
      [[-1, 0], [-0.2, 0.25], [3, 0.25], [5.1, 0.75],
       [6, 0.75], [7, 1], [8, 1]].forEach(function(p) {
         expect(d.p(p[0])).to.equal(p[1]);
         expect(d.p(p[0], true, true)).to.equal(Math.log(p[1]));
         expect(d.p(p[0], false)).to.equal(1 - p[1]);
         expect(d.p(p[0], false, true)).to.equal(Math.log(1 - p[1]));
      });
      [[0, -0.2], [0.2, -0.2], [0.25, -0.2], [0.3, 5.1], [0.45, 5.1], [0.7, 5.1],
       [0.75, 5.1], [0.8, 7], [1, 7]].forEach(function(p) {
         expect(d.q(p[0])).to.equal(p[1]);
         expect(d.q(Math.log(p[0]), true, true)).to.equal(p[1]);
         expect(d.q(1-p[0], false)).to.equal(p[1]);
         expect(d.q(Math.log(1-p[0]), false, true)).to.equal(p[1]);
      });
   });
   it('specified via a function', function() {
      var n, p, f, cdf, d, i, q;
      n = 20;
      for (var times = 0; times < 20; times += 1) {
         p = Math.random();
         f = main.dbinom(n, p);
         cdf = main.pbinom(n, p);
         d = main.finite({ f: f, min: 0, max: n });
         for (i = 0; i <= n; i += 1) {
            expect(utils.relativelyCloseTo(d.d(i), f(i))).to.be.ok;
            expect(utils.relativelyCloseTo(d.p(i), cdf(i))).to.be.ok;
            expect(utils.relativelyCloseTo(d.p(i, false),
               main.pbinom(n, p, false)(i))).to.be.ok;
            if (d.p(i, false) > 1e-8) {
               expect(d.q(d.p(i))).to.equal(i);
            }
            q = main.pbinom(n, p, false)(i);
            if (d.p(i, true) > 1e-8) {
               expect(d.q(d.p(i, false), false)).to.equal(i);
            }
         }
      }
   });
});
