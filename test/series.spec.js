var chai = require('chai');
var expect = chai.expect;
var precision = 1e-10;
var series = require('../panthrMath/utils').series;

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

describe('series function', function() {
   it('passes correct arguments to callback', function() { // n is an integer
      var index, val;
      index = 0;
      function callback(i, v) {
         expect(i).to.equal(index);
         if (i > 0) {
            expect(v).to.equal(val);
         }
         index += 1;
         val = Math.random();
         return val;
      }
      series(callback, 50);

   });
   it('adds terms', function() { // n is an integer
      expect(series(function(i) { return i + 1; }, 50)).to.equal(50 * 51 / 2);
   });
   it('doesn\'t stop too early', function() { // n is missing
      // calculate exp(1) as sum of 1 / n!
      expect(series(function(i, v) {
         return i === 0 ? 1 : v / i;
      })).to.be.relativelyCloseTo(Math.exp(1), precision);
   });
   it('works when given a function for stopping condition', function() {
      var count;
      count = 0;
      function stop() { return count === 5; }
      series(function() { count += 1; return count; }, function() { return count === 5; });
      expect(count).to.equal(5);
   });
});
