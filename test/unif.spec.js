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
   it('is also exported as an object', function() {
      var o;
      o = main.unif(0, 37);
      ['d', 'p', 'q', 'r'].forEach(function(s) {
         expect(o).to.respondTo(s);
         if (s !== 'r') {
            expect(o[s](.2)).to.equal(main[s + 'unif'](0, 37)(.2));
         }
      });
   });
   describe('runif', function() {
      it('returns values in [a, b]', function() {
         var i, f, a, b;
         a = Math.random();
         b = a + 5 * Math.random();
         f = main.runif(a, b);
         for(i = 0; i < 100; i += 1) {
            expect(f()).to.be.within(a, b);
         }
      });
   });
   /* Rcode generating the tests:
      options(digits = 20)
      x = c(NaN, -Inf, -1.3, 1.5, Inf)
      min = c(NaN, -Inf, -1.3, 1.5, Inf)
      max = c(NaN, -Inf, -1, 1.5, 3.7, Inf)
      g = expand.grid(x=x, min=min, max=max)
      g$d = dunif(g$x, g$min, g$max, log=TRUE)
      g$p = punif(g$x, g$min, g$max, log.p=TRUE)
      g$q = punif(g$x, g$min, g$max, lower.tail=FALSE, log.p=TRUE)
      s = paste(g$x, g$min, g$max, g$d, g$p, g$q, sep=", ", collapse="],\n[")
      s = paste("[", s, "]", sep="")
   */
   it('punif, dunif handle inappropriate inputs', function() {
      [
      [NaN, NaN, NaN, NaN, NaN, NaN],
      [-Infinity, NaN, NaN, NaN, NaN, NaN],
      [-1.3, NaN, NaN, NaN, NaN, NaN],
      [1.5, NaN, NaN, NaN, NaN, NaN],
      [Infinity, NaN, NaN, NaN, NaN, NaN],
      [NaN, -Infinity, NaN, NaN, NaN, NaN],
      [-Infinity, -Infinity, NaN, NaN, NaN, NaN],
      [-1.3, -Infinity, NaN, NaN, NaN, NaN],
      [1.5, -Infinity, NaN, NaN, NaN, NaN],
      [Infinity, -Infinity, NaN, NaN, NaN, NaN],
      [NaN, -1.3, NaN, NaN, NaN, NaN],
      [-Infinity, -1.3, NaN, NaN, NaN, NaN],
      [-1.3, -1.3, NaN, NaN, NaN, NaN],
      [1.5, -1.3, NaN, NaN, NaN, NaN],
      [Infinity, -1.3, NaN, NaN, NaN, NaN],
      [NaN, 1.5, NaN, NaN, NaN, NaN],
      [-Infinity, 1.5, NaN, NaN, NaN, NaN],
      [-1.3, 1.5, NaN, NaN, NaN, NaN],
      [1.5, 1.5, NaN, NaN, NaN, NaN],
      [Infinity, 1.5, NaN, NaN, NaN, NaN],
      [NaN, Infinity, NaN, NaN, NaN, NaN],
      [-Infinity, Infinity, NaN, NaN, NaN, NaN],
      [-1.3, Infinity, NaN, NaN, NaN, NaN],
      [1.5, Infinity, NaN, NaN, NaN, NaN],
      [Infinity, Infinity, NaN, NaN, NaN, NaN],
      [NaN, NaN, -Infinity, NaN, NaN, NaN],
      [-Infinity, NaN, -Infinity, NaN, NaN, NaN],
      [-1.3, NaN, -Infinity, NaN, NaN, NaN],
      [1.5, NaN, -Infinity, NaN, NaN, NaN],
      [Infinity, NaN, -Infinity, NaN, NaN, NaN],
      [NaN, -Infinity, -Infinity, NaN, NaN, NaN],
      [-Infinity, -Infinity, -Infinity, NaN, NaN, NaN],
      [-1.3, -Infinity, -Infinity, NaN, NaN, NaN],
      [1.5, -Infinity, -Infinity, NaN, NaN, NaN],
      [Infinity, -Infinity, -Infinity, NaN, NaN, NaN],
      [NaN, -1.3, -Infinity, NaN, NaN, NaN],
      [-Infinity, -1.3, -Infinity, NaN, NaN, NaN],
      [-1.3, -1.3, -Infinity, NaN, NaN, NaN],
      [1.5, -1.3, -Infinity, NaN, NaN, NaN],
      [Infinity, -1.3, -Infinity, NaN, NaN, NaN],
      [NaN, 1.5, -Infinity, NaN, NaN, NaN],
      [-Infinity, 1.5, -Infinity, NaN, NaN, NaN],
      [-1.3, 1.5, -Infinity, NaN, NaN, NaN],
      [1.5, 1.5, -Infinity, NaN, NaN, NaN],
      [Infinity, 1.5, -Infinity, NaN, NaN, NaN],
      [NaN, Infinity, -Infinity, NaN, NaN, NaN],
      [-Infinity, Infinity, -Infinity, NaN, NaN, NaN],
      [-1.3, Infinity, -Infinity, NaN, NaN, NaN],
      [1.5, Infinity, -Infinity, NaN, NaN, NaN],
      [Infinity, Infinity, -Infinity, NaN, NaN, NaN],
      [NaN, NaN, -1, NaN, NaN, NaN],
      [-Infinity, NaN, -1, NaN, NaN, NaN],
      [-1.3, NaN, -1, NaN, NaN, NaN],
      [1.5, NaN, -1, NaN, NaN, NaN],
      [Infinity, NaN, -1, NaN, NaN, NaN],
      [NaN, -Infinity, -1, NaN, NaN, NaN],
      [-Infinity, -Infinity, -1, -Infinity, NaN, NaN],
      [-1.3, -Infinity, -1, -Infinity, NaN, NaN],
      [1.5, -Infinity, -1, -Infinity, NaN, NaN],
      [Infinity, -Infinity, -1, -Infinity, NaN, NaN],
      [NaN, -1.3, -1, NaN, NaN, NaN],
      [-Infinity, -1.3, -1, -Infinity, -Infinity, 0],
      [-1.3, -1.3, -1, 1.20397280432594, -Infinity, 0],
      [1.5, -1.3, -1, -Infinity, 0, -Infinity],
      [Infinity, -1.3, -1, -Infinity, 0, -Infinity],
      [NaN, 1.5, -1, NaN, NaN, NaN],
      [-Infinity, 1.5, -1, NaN, NaN, NaN],
      [-1.3, 1.5, -1, NaN, NaN, NaN],
      [1.5, 1.5, -1, NaN, NaN, NaN],
      [Infinity, 1.5, -1, NaN, NaN, NaN],
      [NaN, Infinity, -1, NaN, NaN, NaN],
      [-Infinity, Infinity, -1, NaN, NaN, NaN],
      [-1.3, Infinity, -1, NaN, NaN, NaN],
      [1.5, Infinity, -1, NaN, NaN, NaN],
      [Infinity, Infinity, -1, NaN, NaN, NaN],
      [NaN, NaN, 1.5, NaN, NaN, NaN],
      [-Infinity, NaN, 1.5, NaN, NaN, NaN],
      [-1.3, NaN, 1.5, NaN, NaN, NaN],
      [1.5, NaN, 1.5, NaN, NaN, NaN],
      [Infinity, NaN, 1.5, NaN, NaN, NaN],
      [NaN, -Infinity, 1.5, NaN, NaN, NaN],
      [-Infinity, -Infinity, 1.5, -Infinity, NaN, NaN],
      [-1.3, -Infinity, 1.5, -Infinity, NaN, NaN],
      [1.5, -Infinity, 1.5, -Infinity, NaN, NaN],
      [Infinity, -Infinity, 1.5, -Infinity, NaN, NaN],
      [NaN, -1.3, 1.5, NaN, NaN, NaN],
      [-Infinity, -1.3, 1.5, -Infinity, -Infinity, 0],
      [-1.3, -1.3, 1.5, -1.02961941718116, -Infinity, 0],
      [1.5, -1.3, 1.5, -1.02961941718116, 0, -Infinity],
      [Infinity, -1.3, 1.5, -Infinity, 0, -Infinity],
      [NaN, 1.5, 1.5, NaN, NaN, NaN],
      [-Infinity, 1.5, 1.5, NaN, -Infinity, 0],
      [-1.3, 1.5, 1.5, NaN, -Infinity, 0],
      [1.5, 1.5, 1.5, NaN, 0, -Infinity],
      [Infinity, 1.5, 1.5, NaN, 0, -Infinity],
      [NaN, Infinity, 1.5, NaN, NaN, NaN],
      [-Infinity, Infinity, 1.5, NaN, NaN, NaN],
      [-1.3, Infinity, 1.5, NaN, NaN, NaN],
      [1.5, Infinity, 1.5, NaN, NaN, NaN],
      [Infinity, Infinity, 1.5, NaN, NaN, NaN],
      [NaN, NaN, 3.7, NaN, NaN, NaN],
      [-Infinity, NaN, 3.7, NaN, NaN, NaN],
      [-1.3, NaN, 3.7, NaN, NaN, NaN],
      [1.5, NaN, 3.7, NaN, NaN, NaN],
      [Infinity, NaN, 3.7, NaN, NaN, NaN],
      [NaN, -Infinity, 3.7, NaN, NaN, NaN],
      [-Infinity, -Infinity, 3.7, -Infinity, NaN, NaN],
      [-1.3, -Infinity, 3.7, -Infinity, NaN, NaN],
      [1.5, -Infinity, 3.7, -Infinity, NaN, NaN],
      [Infinity, -Infinity, 3.7, -Infinity, NaN, NaN],
      [NaN, -1.3, 3.7, NaN, NaN, NaN],
      [-Infinity, -1.3, 3.7, -Infinity, -Infinity, 0],
      [-1.3, -1.3, 3.7, -1.6094379124341, -Infinity, 0],
      [1.5, -1.3, 3.7, -1.6094379124341, -0.579818495252942, -0.82098055206983],
      [Infinity, -1.3, 3.7, -Infinity, 0, -Infinity],
      [NaN, 1.5, 3.7, NaN, NaN, NaN],
      [-Infinity, 1.5, 3.7, -Infinity, -Infinity, 0],
      [-1.3, 1.5, 3.7, -Infinity, -Infinity, 0],
      [1.5, 1.5, 3.7, -0.78845736036427, -Infinity, 0],
      [Infinity, 1.5, 3.7, -Infinity, 0, -Infinity],
      [NaN, Infinity, 3.7, NaN, NaN, NaN],
      [-Infinity, Infinity, 3.7, NaN, NaN, NaN],
      [-1.3, Infinity, 3.7, NaN, NaN, NaN],
      [1.5, Infinity, 3.7, NaN, NaN, NaN],
      [Infinity, Infinity, 3.7, NaN, NaN, NaN],
      [NaN, NaN, Infinity, NaN, NaN, NaN],
      [-Infinity, NaN, Infinity, NaN, NaN, NaN],
      [-1.3, NaN, Infinity, NaN, NaN, NaN],
      [1.5, NaN, Infinity, NaN, NaN, NaN],
      [Infinity, NaN, Infinity, NaN, NaN, NaN],
      [NaN, -Infinity, Infinity, NaN, NaN, NaN],
      [-Infinity, -Infinity, Infinity, -Infinity, NaN, NaN],
      [-1.3, -Infinity, Infinity, -Infinity, NaN, NaN],
      [1.5, -Infinity, Infinity, -Infinity, NaN, NaN],
      [Infinity, -Infinity, Infinity, -Infinity, NaN, NaN],
      [NaN, -1.3, Infinity, NaN, NaN, NaN],
      [-Infinity, -1.3, Infinity, -Infinity, NaN, NaN],
      [-1.3, -1.3, Infinity, -Infinity, NaN, NaN],
      [1.5, -1.3, Infinity, -Infinity, NaN, NaN],
      [Infinity, -1.3, Infinity, -Infinity, NaN, NaN],
      [NaN, 1.5, Infinity, NaN, NaN, NaN],
      [-Infinity, 1.5, Infinity, -Infinity, NaN, NaN],
      [-1.3, 1.5, Infinity, -Infinity, NaN, NaN],
      [1.5, 1.5, Infinity, -Infinity, NaN, NaN],
      [Infinity, 1.5, Infinity, -Infinity, NaN, NaN],
      [NaN, Infinity, Infinity, NaN, NaN, NaN],
      [-Infinity, Infinity, Infinity, NaN, NaN, NaN],
      [-1.3, Infinity, Infinity, NaN, NaN, NaN],
      [1.5, Infinity, Infinity, NaN, NaN, NaN],
      [Infinity, Infinity, Infinity, NaN, NaN, NaN]
      ].forEach(function(tuple) {
         var x, min, max, rlogp, rlogq, rlogd, logp, logq, logd;
         x = tuple[0];
         min = tuple[1];
         max = tuple[2];
         rlogd = tuple[3];
         rlogp = tuple[4];
         rlogq = tuple[5];
         logp = main.punif(min, max, true, true)(x);
         logq = main.punif(min, max, false, true)(x);
         logd = main.dunif(min, max, true)(x);
         expect(utils.relativelyCloseTo(logp, rlogp)).to.equal(true);
         expect(utils.relativelyCloseTo(logq, rlogq)).to.equal(true);
         expect(utils.relativelyCloseTo(logd, rlogd)).to.equal(true);
      });
   });

   /* Rcode generating the tests:
      options(digits = 20)
      p = c(NaN, -Inf, -1, 0, 0.3, 1, Inf)
      min = c(NaN, -Inf, -1.3, 1.5, Inf)
      max = c(NaN, -Inf, -1, 1.5, 3.7, Inf)
      g = expand.grid(p=p, min=min, max=max)
      g$x1 = qunif(g$p, g$min, g$max, lower.tail=TRUE, log.p=FALSE)
      g$x2 = qunif(g$p, g$min, g$max, lower.tail=FALSE, log.p=FALSE)
      s = paste(g$p, g$min, g$max, g$x1, g$x2, sep=", ", collapse="],\n[")
      s = paste("[", s, "]", sep="")
   */
   it('qunif handles inappropriate inputs', function() {
      [
      [NaN, NaN, NaN, NaN, NaN],
      [-Infinity, NaN, NaN, NaN, NaN],
      [-1, NaN, NaN, NaN, NaN],
      [0, NaN, NaN, NaN, NaN],
      [0.3, NaN, NaN, NaN, NaN],
      [1, NaN, NaN, NaN, NaN],
      [Infinity, NaN, NaN, NaN, NaN],
      [NaN, -Infinity, NaN, NaN, NaN],
      [-Infinity, -Infinity, NaN, NaN, NaN],
      [-1, -Infinity, NaN, NaN, NaN],
      [0, -Infinity, NaN, NaN, NaN],
      [0.3, -Infinity, NaN, NaN, NaN],
      [1, -Infinity, NaN, NaN, NaN],
      [Infinity, -Infinity, NaN, NaN, NaN],
      [NaN, -1.3, NaN, NaN, NaN],
      [-Infinity, -1.3, NaN, NaN, NaN],
      [-1, -1.3, NaN, NaN, NaN],
      [0, -1.3, NaN, NaN, NaN],
      [0.3, -1.3, NaN, NaN, NaN],
      [1, -1.3, NaN, NaN, NaN],
      [Infinity, -1.3, NaN, NaN, NaN],
      [NaN, 1.5, NaN, NaN, NaN],
      [-Infinity, 1.5, NaN, NaN, NaN],
      [-1, 1.5, NaN, NaN, NaN],
      [0, 1.5, NaN, NaN, NaN],
      [0.3, 1.5, NaN, NaN, NaN],
      [1, 1.5, NaN, NaN, NaN],
      [Infinity, 1.5, NaN, NaN, NaN],
      [NaN, Infinity, NaN, NaN, NaN],
      [-Infinity, Infinity, NaN, NaN, NaN],
      [-1, Infinity, NaN, NaN, NaN],
      [0, Infinity, NaN, NaN, NaN],
      [0.3, Infinity, NaN, NaN, NaN],
      [1, Infinity, NaN, NaN, NaN],
      [Infinity, Infinity, NaN, NaN, NaN],
      [NaN, NaN, -Infinity, NaN, NaN],
      [-Infinity, NaN, -Infinity, NaN, NaN],
      [-1, NaN, -Infinity, NaN, NaN],
      [0, NaN, -Infinity, NaN, NaN],
      [0.3, NaN, -Infinity, NaN, NaN],
      [1, NaN, -Infinity, NaN, NaN],
      [Infinity, NaN, -Infinity, NaN, NaN],
      [NaN, -Infinity, -Infinity, NaN, NaN],
      [-Infinity, -Infinity, -Infinity, NaN, NaN],
      [-1, -Infinity, -Infinity, NaN, NaN],
      [0, -Infinity, -Infinity, NaN, NaN],
      [0.3, -Infinity, -Infinity, NaN, NaN],
      [1, -Infinity, -Infinity, NaN, NaN],
      [Infinity, -Infinity, -Infinity, NaN, NaN],
      [NaN, -1.3, -Infinity, NaN, NaN],
      [-Infinity, -1.3, -Infinity, NaN, NaN],
      [-1, -1.3, -Infinity, NaN, NaN],
      [0, -1.3, -Infinity, NaN, NaN],
      [0.3, -1.3, -Infinity, NaN, NaN],
      [1, -1.3, -Infinity, NaN, NaN],
      [Infinity, -1.3, -Infinity, NaN, NaN],
      [NaN, 1.5, -Infinity, NaN, NaN],
      [-Infinity, 1.5, -Infinity, NaN, NaN],
      [-1, 1.5, -Infinity, NaN, NaN],
      [0, 1.5, -Infinity, NaN, NaN],
      [0.3, 1.5, -Infinity, NaN, NaN],
      [1, 1.5, -Infinity, NaN, NaN],
      [Infinity, 1.5, -Infinity, NaN, NaN],
      [NaN, Infinity, -Infinity, NaN, NaN],
      [-Infinity, Infinity, -Infinity, NaN, NaN],
      [-1, Infinity, -Infinity, NaN, NaN],
      [0, Infinity, -Infinity, NaN, NaN],
      [0.3, Infinity, -Infinity, NaN, NaN],
      [1, Infinity, -Infinity, NaN, NaN],
      [Infinity, Infinity, -Infinity, NaN, NaN],
      [NaN, NaN, -1, NaN, NaN],
      [-Infinity, NaN, -1, NaN, NaN],
      [-1, NaN, -1, NaN, NaN],
      [0, NaN, -1, NaN, NaN],
      [0.3, NaN, -1, NaN, NaN],
      [1, NaN, -1, NaN, NaN],
      [Infinity, NaN, -1, NaN, NaN],
      [NaN, -Infinity, -1, NaN, NaN],
      [-Infinity, -Infinity, -1, NaN, NaN],
      [-1, -Infinity, -1, NaN, NaN],
      [0, -Infinity, -1, NaN, NaN],
      [0.3, -Infinity, -1, NaN, NaN],
      [1, -Infinity, -1, NaN, NaN],
      [Infinity, -Infinity, -1, NaN, NaN],
      [NaN, -1.3, -1, NaN, NaN],
      [-Infinity, -1.3, -1, NaN, NaN],
      [-1, -1.3, -1, NaN, NaN],
      [0, -1.3, -1, -1.3, -1],
      [0.3, -1.3, -1, -1.21, -1.09],
      [1, -1.3, -1, -1, -1.3],
      [Infinity, -1.3, -1, NaN, NaN],
      [NaN, 1.5, -1, NaN, NaN],
      [-Infinity, 1.5, -1, NaN, NaN],
      [-1, 1.5, -1, NaN, NaN],
      [0, 1.5, -1, NaN, NaN],
      [0.3, 1.5, -1, NaN, NaN],
      [1, 1.5, -1, NaN, NaN],
      [Infinity, 1.5, -1, NaN, NaN],
      [NaN, Infinity, -1, NaN, NaN],
      [-Infinity, Infinity, -1, NaN, NaN],
      [-1, Infinity, -1, NaN, NaN],
      [0, Infinity, -1, NaN, NaN],
      [0.3, Infinity, -1, NaN, NaN],
      [1, Infinity, -1, NaN, NaN],
      [Infinity, Infinity, -1, NaN, NaN],
      [NaN, NaN, 1.5, NaN, NaN],
      [-Infinity, NaN, 1.5, NaN, NaN],
      [-1, NaN, 1.5, NaN, NaN],
      [0, NaN, 1.5, NaN, NaN],
      [0.3, NaN, 1.5, NaN, NaN],
      [1, NaN, 1.5, NaN, NaN],
      [Infinity, NaN, 1.5, NaN, NaN],
      [NaN, -Infinity, 1.5, NaN, NaN],
      [-Infinity, -Infinity, 1.5, NaN, NaN],
      [-1, -Infinity, 1.5, NaN, NaN],
      [0, -Infinity, 1.5, NaN, NaN],
      [0.3, -Infinity, 1.5, NaN, NaN],
      [1, -Infinity, 1.5, NaN, NaN],
      [Infinity, -Infinity, 1.5, NaN, NaN],
      [NaN, -1.3, 1.5, NaN, NaN],
      [-Infinity, -1.3, 1.5, NaN, NaN],
      [-1, -1.3, 1.5, NaN, NaN],
      [0, -1.3, 1.5, -1.3, 1.5],
      [0.3, -1.3, 1.5, -0.46, 0.66],
      [1, -1.3, 1.5, 1.5, -1.3],
      [Infinity, -1.3, 1.5, NaN, NaN],
      [NaN, 1.5, 1.5, NaN, NaN],
      [-Infinity, 1.5, 1.5, NaN, NaN],
      [-1, 1.5, 1.5, NaN, NaN],
      [0, 1.5, 1.5, 1.5, 1.5],
      [0.3, 1.5, 1.5, 1.5, 1.5],
      [1, 1.5, 1.5, 1.5, 1.5],
      [Infinity, 1.5, 1.5, NaN, NaN],
      [NaN, Infinity, 1.5, NaN, NaN],
      [-Infinity, Infinity, 1.5, NaN, NaN],
      [-1, Infinity, 1.5, NaN, NaN],
      [0, Infinity, 1.5, NaN, NaN],
      [0.3, Infinity, 1.5, NaN, NaN],
      [1, Infinity, 1.5, NaN, NaN],
      [Infinity, Infinity, 1.5, NaN, NaN],
      [NaN, NaN, 3.7, NaN, NaN],
      [-Infinity, NaN, 3.7, NaN, NaN],
      [-1, NaN, 3.7, NaN, NaN],
      [0, NaN, 3.7, NaN, NaN],
      [0.3, NaN, 3.7, NaN, NaN],
      [1, NaN, 3.7, NaN, NaN],
      [Infinity, NaN, 3.7, NaN, NaN],
      [NaN, -Infinity, 3.7, NaN, NaN],
      [-Infinity, -Infinity, 3.7, NaN, NaN],
      [-1, -Infinity, 3.7, NaN, NaN],
      [0, -Infinity, 3.7, NaN, NaN],
      [0.3, -Infinity, 3.7, NaN, NaN],
      [1, -Infinity, 3.7, NaN, NaN],
      [Infinity, -Infinity, 3.7, NaN, NaN],
      [NaN, -1.3, 3.7, NaN, NaN],
      [-Infinity, -1.3, 3.7, NaN, NaN],
      [-1, -1.3, 3.7, NaN, NaN],
      [0, -1.3, 3.7, -1.3, 3.7],
      [0.3, -1.3, 3.7, 0.2, 2.2],
      [1, -1.3, 3.7, 3.7, -1.3],
      [Infinity, -1.3, 3.7, NaN, NaN],
      [NaN, 1.5, 3.7, NaN, NaN],
      [-Infinity, 1.5, 3.7, NaN, NaN],
      [-1, 1.5, 3.7, NaN, NaN],
      [0, 1.5, 3.7, 1.5, 3.7],
      [0.3, 1.5, 3.7, 2.16, 3.04],
      [1, 1.5, 3.7, 3.7, 1.5],
      [Infinity, 1.5, 3.7, NaN, NaN],
      [NaN, Infinity, 3.7, NaN, NaN],
      [-Infinity, Infinity, 3.7, NaN, NaN],
      [-1, Infinity, 3.7, NaN, NaN],
      [0, Infinity, 3.7, NaN, NaN],
      [0.3, Infinity, 3.7, NaN, NaN],
      [1, Infinity, 3.7, NaN, NaN],
      [Infinity, Infinity, 3.7, NaN, NaN],
      [NaN, NaN, Infinity, NaN, NaN],
      [-Infinity, NaN, Infinity, NaN, NaN],
      [-1, NaN, Infinity, NaN, NaN],
      [0, NaN, Infinity, NaN, NaN],
      [0.3, NaN, Infinity, NaN, NaN],
      [1, NaN, Infinity, NaN, NaN],
      [Infinity, NaN, Infinity, NaN, NaN],
      [NaN, -Infinity, Infinity, NaN, NaN],
      [-Infinity, -Infinity, Infinity, NaN, NaN],
      [-1, -Infinity, Infinity, NaN, NaN],
      [0, -Infinity, Infinity, NaN, NaN],
      [0.3, -Infinity, Infinity, NaN, NaN],
      [1, -Infinity, Infinity, NaN, NaN],
      [Infinity, -Infinity, Infinity, NaN, NaN],
      [NaN, -1.3, Infinity, NaN, NaN],
      [-Infinity, -1.3, Infinity, NaN, NaN],
      [-1, -1.3, Infinity, NaN, NaN],
      [0, -1.3, Infinity, NaN, NaN],
      [0.3, -1.3, Infinity, NaN, NaN],
      [1, -1.3, Infinity, NaN, NaN],
      [Infinity, -1.3, Infinity, NaN, NaN],
      [NaN, 1.5, Infinity, NaN, NaN],
      [-Infinity, 1.5, Infinity, NaN, NaN],
      [-1, 1.5, Infinity, NaN, NaN],
      [0, 1.5, Infinity, NaN, NaN],
      [0.3, 1.5, Infinity, NaN, NaN],
      [1, 1.5, Infinity, NaN, NaN],
      [Infinity, 1.5, Infinity, NaN, NaN],
      [NaN, Infinity, Infinity, NaN, NaN],
      [-Infinity, Infinity, Infinity, NaN, NaN],
      [-1, Infinity, Infinity, NaN, NaN],
      [0, Infinity, Infinity, NaN, NaN],
      [0.3, Infinity, Infinity, NaN, NaN],
      [1, Infinity, Infinity, NaN, NaN],
      [Infinity, Infinity, Infinity, NaN, NaN]
      ].forEach(function(tuple) {
         var p, min, max, x1r, x2r, x1, x2;
         p = tuple[0];
         min = tuple[1];
         max = tuple[2];
         x1r = tuple[3];
         x2r = tuple[4];
         x1 = main.qunif(min, max, true, false)(p);
         x2 = main.qunif(min, max, false, false)(p);
         expect(utils.relativelyCloseTo(x1r, x1)).to.equal(true);
         expect(utils.relativelyCloseTo(x2r, x2)).to.equal(true);
      });
   });
});

