var chai = require('chai');
var expect = chai.expect;
var precision = 1e-8;
var lgamma = require('../panthrMath/basicFunc/lgamma').lgamma;

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

describe('lgamma function', function() {
   it('works', function() {
      [
       [0.9,0.0663762397347429],
       [0.99,0.00585480676470978],
       [0.999,0.000578038532891351],
       [0.9999,5.7729791561249e-05],
       [0.99999,5.77223889613319e-06],
       [0.999999,5.77216487440187e-07],
       [0.9999999,5.77215747415977e-08],
       [1.1,-0.0498724412598397],
       [1.01,-0.00569030794606961],
       [1.001,-0.000576393598283262],
       [1.0001,-5.7713342220443e-05],
       [1.00001,-5.77207440282475e-06],
       [1.000001,-5.77214842422725e-07],
       [1.0000001,-5.77215583114071e-08],
       // [1,0],
       // [2,0],
       [3,0.693147180559945],
       [4,1.79175946922805],
       [5,3.17805383034795],
       [6,4.78749174278205],
       [7,6.5792512120101],
       [8,8.52516136106541],
       [9,10.6046029027453],
       [10,12.8018274800815],
       [0.1,2.25271265173421],
       [0.01,4.59947987804202],
       [0.001,6.90717888538385],
       [1e-04,9.21028265863396],
       [1e-05,11.5129196928958],
       [1e-06,13.8155099807494],
       [1e-07,16.1180955932368],
       [1.9,-0.0389842759230833],
       [1.99,-0.00419552908879171],
       [1.999,-0.00042246180069209],
       [1.9999,-4.22752087721528e-05],
       [1.99999,-4.22781110431458e-06],
       [1.999999,-4.22784012688221e-07],
       [1.9999999,-4.22784302067691e-08],
       [2.1,0.0454377385444851],
       [2.01,0.00426002290709843],
       [2.001,0.000423106734800085],
       [2.0001,4.22816581126612e-05],
       [2.00001,4.22787559765649e-06],
       [2.000001,4.22784657616089e-07],
       [2.0000001,4.22784367459758e-08],
       [-0.8,1.74720737374499],
       [-1.8,1.15942070884288],
       [-2.8,0.129801291661716],
       [-3.8,-1.20519977507062],
       [-4.8,-2.77381569298447],
       [-5.8,-4.53167361053684],
       [-6.8,-6.4485962227189],
       [-7.8,-8.50271995641445],
       [-8.8,-10.6774716778986],
       [-9.8,-12.9598540635751]
      ].forEach(function(pair) {
      console.log(pair, lgamma(pair[0]));
      expect(lgamma(pair[0]))
         .to.be.relativelyCloseTo(pair[1], precision);
      });
   });
});
