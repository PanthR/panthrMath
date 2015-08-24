var chai = require('chai');
var expect = chai.expect;
var precision = 1e-10;
var contFrac = require('../panthrMath/utils').contFrac;

chai.use(function(_chai, utils) {
   var Assertion = _chai.Assertion;
   Assertion.addMethod('relativelyCloseTo', function(x0, delta) {
      var x = utils.flag(this, 'object');
      var denom = Math.max(Math.abs(x0), Math.abs(x));
      if (denom !== 0) {
         var res = Math.abs(x - x0) / denom;
         new Assertion(res).to.be.below(delta);
      }
   });
});

describe('contFrac function', function() {
   it('passes correct arguments to callbacks', function() { // n is an integer

      function buildCB(index, val) {
         return function(i, v) {
            expect(i).to.equal(index);
            if (i > 0) {
               expect(v).to.equal(val);
            }
            index += 1;
            val = Math.random();
            return val;
         }
      }
      contFrac(buildCB(0), buildCB(1, 1), 50);
   });
   it('approximates the cont frac when stop is an integer', function() {
      function f() { return 1; }
      expect(contFrac(f, f, 5)).to.equal(8 / 5);
   });
   it('approximates the cont frac when stop is not given', function() {
      // golden ratio has an = bn = 1
      function f() { return 1; }
      expect(contFrac(f, f)).to.be
         .relativelyCloseTo((1 + Math.sqrt(5))/2, precision);
      // Square root of 2
      expect(contFrac(function(i) { return i === 0 ? 1 : 2; }, f)).to.be
         .relativelyCloseTo(Math.sqrt(2), precision);
   });
   it('works when given a function for stopping condition', function() {
      var count;
      count = 0;
      contFrac(function() { count += 1; return count; },
         function() { return 0; },
         function() { return count === 5; });
      expect(count).to.equal(5);
   });
});
