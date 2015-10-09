var chai = require('chai');
var expect = chai.expect;
var precision = 1e-10;
var contFrac = require('../panthrMath/utils').contFrac;
var utils = require('../panthrMath/utils');

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
      expect(utils.relativelyCloseTo(
         contFrac(f, f), (1 + Math.sqrt(5))/2, precision)).to.be.ok;
      // Square root of 2
      expect(utils.relativelyCloseTo(
         contFrac(function(i) { return i === 0 ? 1 : 2; }, f), Math.sqrt(2), precision)).to.be.ok;
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
