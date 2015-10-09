var chai = require('chai');
var expect = chai.expect;
var precision = 1e-10;
var erf = require('../panthrMath/basicFunc/erf').erf;
var erfc = require('../panthrMath/basicFunc/erf').erfc;
var utils = require('../panthrMath/utils');

describe('erf function', function() {
   it('erf', function() {
      [[-5,-0.999999999998463],
       [-4,-0.999999984582742],
       [-3,-0.999977909503001],
       [-2,-0.995322265018953],
       [-1,-0.842700792949715],
       [1,0.842700792949715],
       [2,0.995322265018953],
       [3,0.999977909503001],
       [4,0.999999984582742],
       [5,0.999999999998463],
       [1e-05,1.1283791670591e-05]
      ].forEach(function(pair) {
      expect(utils.relativelyCloseTo(erf(pair[0]), pair[1], precision)).to.be.ok;
      });
   });
   it('erfc', function() {
      [[-5,1.99999999999846],
       [-4,1.99999998458274],
       [-3,1.999977909503],
       [-2,1.99532226501895],
       [-1,1.84270079294972],
       [1,0.157299207050285],
       [2,0.0046777349810473],
       [3,2.2090496998585e-05],
       [4,1.541725790028e-08],
       [5,1.537459794428e-12],
       [1e-05,0.999988716208329]
      ].forEach(function(pair) {
      expect(utils.relativelyCloseTo(erfc(pair[0]), pair[1], precision)).to.be.ok;
      });
   });
});
