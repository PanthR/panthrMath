var chai = require('chai');
var expect = chai.expect;
var precision = 1e-8;
var log1p = require('../panthrMath/basicFunc/log1p').log1p;

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

describe('log1p function', function() {
   it('works', function() {
      [[-0.1,-0.105360515657826],
       [-0.01,-0.0100503358535015],
       [-0.001,-0.00100050033358353],
       [-1e-04,-0.000100005000333347],
       [-1e-05,-1.00000500002878e-05],
       [-1e-06,-1.00000050002909e-06],
       [-1e-07,-1.00000004947365e-07],
       [0.1,0.0953101798043249],
       [0.01,0.00995033085316809],
       [0.001,0.000999500333083423],
       [1e-04,9.99950003332973e-05],
       [1e-05,9.99995000039884e-06],
       [1e-06,9.99999499918067e-07],
       [1e-07,9.9999995058387e-08],
       [1,0.693147180559945],
       [2,1.09861228866811],
       [3,1.38629436111989],
       [4,1.6094379124341],
       [5,1.79175946922805]
      ].forEach(function(pair) {
      expect(log1p(pair[0]))
         .to.be.relativelyCloseTo(pair[1], precision);
      });
   });
});
