var t = require('../panthrMath/distributions/t');
var chai = require('chai');
var expect = chai.expect;

var precision = 1e-5;

chai.use(function(_chai, utils) {
   var Assertion = _chai.Assertion;

   Assertion.addMethod('relativelyCloseTo', function(x0, delta) {
      var x = utils.flag(this, 'object');
      var denom = Math.max(Math.abs(x0), Math.abs(x));
      var res = Math.abs(x - x0) / denom;
      new Assertion(res).to.be.below(delta);
   });
});

describe('T Distribution', function() {
   it('dtlog', function() {
      [[-5, -5.93333292517818],
       [-4, -5.00442403409698],
       [-3, -3.92746674386584],
       [-2, -2.71369720441159],
       [-1, -1.53868813129725],
       [ 0, -0.980829253011726],
       [ 1, -1.53868813129725],
       [ 2, -2.71369720441159],
       [ 3, -3.92746674386584],
       [ 4, -5.00442403409698],
       [ 5, -5.93333292517818]
      ].forEach(function(pair) {
      expect(t.dt(4, true)(pair[0])) // dt.log
         .to.be.relativelyCloseTo(pair[1], precision);
      });
   });
   it('dt', function() {
      [[-5, 0.00264963621655722],
       [-4, 0.00670820393249937],
       [-3, 0.0196934980908365],
       [-2, 0.0662912607362388],
       [-1, 0.21466252583998],
       [ 0, 0.375],
       [ 1, 0.21466252583998],
       [ 2, 0.0662912607362388],
       [ 3, 0.0196934980908365],
       [ 4, 0.00670820393249937],
       [ 5, 0.00264963621655722]
      ].forEach(function(pair) {
      expect(t.dt(4)(pair[0]))
         .to.be.relativelyCloseTo(pair[1], precision);
      });
   });
});
