var chai = require('chai');
var expect = chai.expect;
var precision = 1e-10;
var log1p = require('../panthrMath/basicFunc/log1p').log1p;
var utils = require('../panthrMath/utils');

describe('log1p function', function() {
   it('works', function() {
      [[-0.1,-0.105360515657826],
       [-0.01,-0.0100503358535014],
       [-0.001,-0.00100050033358353],
       [-1e-04,-0.000100005000333358],
       [-1e-05,-1.00000500003333e-05],
       [-1e-06,-1.00000050000033e-06],
       [-1e-07,-1.00000005e-07],
       [0.1,0.0953101798043249],
       [0.01,0.00995033085316808],
       [0.001,0.000999500333083533],
       [1e-04,9.99950003333083e-05],
       [1e-05,9.99995000033333e-06],
       [1e-06,9.99999500000333e-07],
       [1e-07,9.99999950000003e-08],
       [1,0.693147180559945],
       [2,1.09861228866811],
       [3,1.38629436111989],
       [4,1.6094379124341],
       [5,1.79175946922805]
      ].forEach(function(pair) {
      expect(utils.relativelyCloseTo(log1p(pair[0]),pair[1], precision)).to.be.ok;
      });
   });
});
