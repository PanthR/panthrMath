(function(define) {'use strict';
define(function(require) {

   // error term in Stirling's approximation
   // log(n!) - log( sqrt(2*pi*n) * (n/e)^n )

   var cs, precomputed, lgamma, C;

   lgamma = require('./lgamma').lgamma;
   C = require('../constants');

   // coefficients in Stirling's expansion for log(Gamma)
   cs = [1 / 12, 1 / 360, 1 / 1260, 1 / 1680, 1 / 1188];
   precomputed = [
      NaN, /* place holder only */
      0.1534264097200273452913848,  /* 0.5 */
      0.0810614667953272582196702,  /* 1.0 */
      0.0548141210519176538961390,  /* 1.5 */
      0.0413406959554092940938221,  /* 2.0 */
      0.03316287351993628748511048, /* 2.5 */
      0.02767792568499833914878929, /* 3.0 */
      0.02374616365629749597132920, /* 3.5 */
      0.02079067210376509311152277, /* 4.0 */
      0.01848845053267318523077934, /* 4.5 */
      0.01664469118982119216319487, /* 5.0 */
      0.01513497322191737887351255, /* 5.5 */
      0.01387612882307074799874573, /* 6.0 */
      0.01281046524292022692424986, /* 6.5 */
      0.01189670994589177009505572, /* 7.0 */
      0.01110455975820691732662991, /* 7.5 */
      0.010411265261972096497478567, /* 8.0 */
      0.009799416126158803298389475, /* 8.5 */
      0.009255462182712732917728637, /* 9.0 */
      0.008768700134139385462952823, /* 9.5 */
      0.008330563433362871256469318, /* 10.0 */
      0.007934114564314020547248100, /* 10.5 */
      0.007573675487951840794972024, /* 11.0 */
      0.007244554301320383179543912, /* 11.5 */
      0.006942840107209529865664152, /* 12.0 */
      0.006665247032707682442354394, /* 12.5 */
      0.006408994188004207068439631, /* 13.0 */
      0.006171712263039457647532867, /* 13.5 */
      0.005951370112758847735624416, /* 14.0 */
      0.005746216513010115682023589, /* 14.5 */
      0.005554733551962801371038690  /* 15.0 */
   ];

   // error term in Stirling's approximation
   // log(n!) - log( sqrt(2*pi*n) * (n/e)^n )
   // Code adapted from R source, based off Loader (2000)
   // For n > 15, uses series
   // For n <= 15 which are integers or half-integers, uses stored values
   // For other n < 15, uses lgamma.

   /**
    * TODO
    * @memberof basicFunc
    */
   function stirlerr(n) {
      var nsq = n * n;
      if (n <= 15) {
         if (2 * n === Math.floor(2 * n)) {
            return precomputed[2 * n];
         }
         return lgamma(n + 1) -
              (n + 0.5) * Math.log(n) + n - Math.log(C.sqrt2pi);
      }
      if (n > 500) { return (cs[0] - cs[1] / nsq) / n; }
      if (n > 80) { return (cs[0] - (cs[1] - cs[2] / nsq) / nsq) / n; }
      if (n > 35) {return (cs[0] - (cs[1] - (cs[2] - cs[3] / nsq) / nsq) / nsq) / n; }
      return (cs[0] - (cs[1] - (cs[2] - (cs[3] - cs[4] / nsq) / nsq) / nsq) / nsq) / n;
   }

   return { stirlerr: stirlerr };

});  // end define

}(typeof define === 'function' && define.amd ? define : function(factory) {
   'use strict';
   module.exports = factory(require);
}));
