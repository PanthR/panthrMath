var normal = require('../panthrMath/distributions/normal');
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

describe('Normal Distribution', function() {
   it('dnormlog', function() {
      [[0, -0.918938],
       [1, -1.418939],
       [2, -2.918939],
       [3, -5.418939],
       [4, -8.918939],
       [-1, -1.418939],
       [-2, -2.918939],
       [-3, -5.418939],
       [-4, -8.918939],
       [-5, -13.418939]
      ].forEach(function(pair) {
      expect(normal.dnormLog(0, 1)(pair[0]))
         .to.be.relativelyCloseTo(pair[1], precision);
      });
      var mu = 2.3, sigma = 0.0012;
      expect(normal.dnormLog(mu,sigma)(2.2912))
         .to.be.relativelyCloseTo(-21.08239, precision);
   });
   it('dnorm', function() {
      [[0, 3.989423e-01],
       [1, 2.419707e-01],
       [2, 5.399097e-02],
       [3, 4.431848e-03],
       [4, 1.338302e-04],
       [-1, 2.419707e-01],
       [-2, 5.399097e-02],
       [-3, 4.431848e-03],
       [-4, 1.338302e-04],
       [-5, 1.486720e-06]
      ].forEach(function(pair) {
      expect(normal.dnorm(0, 1)(pair[0]))
         .to.be.relativelyCloseTo(pair[1], precision);
      });
      var mu = 2.3, sigma = 0.0012;
      expect(normal.dnorm(mu,sigma)(2.2912))
         .to.be.relativelyCloseTo(6.982851e-10, precision);
   });
});
