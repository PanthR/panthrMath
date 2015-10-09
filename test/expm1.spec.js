var chai = require('chai');
var expect = chai.expect;
var precision = 1e-10;
var expm1 = require('../panthrMath/basicFunc/expm1').expm1;
var utils = require('../panthrMath/utils');

describe('expm1 function', function() {
   it('works', function() {
      [[-0.1,-0.0951625819640404],
       [-0.01,-0.00995016625083195],
       [-0.001,-0.000999500166625008],
       [-1e-04,-9.99950001666625e-05],
       [-1e-05,-9.99995000016667e-06],
       [-1e-06,-9.99999500000167e-07],
       [-1e-07,-9.99999950000002e-08],
       [0.1,0.105170918075648],
       [0.01,0.0100501670841681],
       [0.001,0.00100050016670834],
       [1e-04,0.000100005000166671],
       [1e-05,1.00000500001667e-05],
       [1e-06,1.00000050000017e-06],
       [1e-07,1.00000005e-07],
       [1,1.71828182845905],
       [2,6.38905609893065],
       [3,19.0855369231877],
       [4,53.5981500331442],
       [5,147.413159102577]
      ].forEach(function(pair) {
      expect(utils.relativelyCloseTo(expm1(pair[0]), pair[1], precision)).to.be.ok;
      });
   });
});
