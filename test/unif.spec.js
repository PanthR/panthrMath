var main = require('..');
var chai = require('chai');
var expect = chai.expect;
var utils = require('../panthrMath/utils');

var precision = 1e-12;

describe('Uniform Distribution', function() {
   it('dunif', function() {
      [
         [4, 5, 7, 0],
         [6, 5, 7, 0.5],
         [8, 5, 7, 0]
      ].forEach(function(tuple) {
         var x, a, b, logp, rlogp;
         x = tuple[0];
         a = tuple[1];
         b = tuple[2];
         rlogp = tuple[3];
         logp = main.dunif(a, b)(x);
         expect(utils.relativelyCloseTo(rlogp, logp, precision)).to.be.ok;
      });
   });

   it('punif and qunif', function() {
      [
         [4, 5, 7, -Infinity, 0],
         [6, 5, 7, Math.log(0.5), Math.log(0.5)],
         [6, 5, 8, Math.log(1/3), Math.log(2/3)],
         [8, 5, 7, 0, -Infinity]
      ].forEach(function(tuple) {
         var x, a, b, logp, rlogp, rlogq, xfromqbeta;
         x = tuple[0];
         a = tuple[1];
         b = tuple[2];
         rlogp = tuple[3];
         rlogq = tuple[4];
         logp = main.punif(a, b, true, true)(x);
         expect(utils.relativelyCloseTo(rlogp, logp, precision)).to.be.ok;
         expect(utils.relativelyCloseTo(Math.exp(rlogp),
            main.punif(a, b)(x), precision)).to.be.ok;
         expect(utils.relativelyCloseTo(rlogq,
            main.punif(a, b, false, true)(x), precision)).to.be.ok;
         expect(utils.relativelyCloseTo(Math.exp(rlogq),
            main.punif(a, b, false)(x), precision)).to.be.ok;
         if (logp > -Infinity && logp < 0) {
            expect(utils.relativelyCloseTo(
               main.qunif(a,b, true, true)(rlogp),
               x, precision)).to.be.ok;
            expect(utils.relativelyCloseTo(
               main.qunif(a,b, false, true)(rlogq),
               x, precision)).to.be.ok;
         }
      });
   });
});
