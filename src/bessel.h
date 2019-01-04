#ifndef BESSEL_H
#define BESSEL_H

double maxi(double x, double y);

double BESSI0(double X);
double BESSI1(double X);



inline double maxi(double x, double y)
{
  if(x > y) return(x);
  return(y);
}


// Returns the modified Bessel function I_0 (x) for any real x
inline double BESSI0(double x)
{
return R::bessel_i(fabs(x), 0, 1);
}

// Returns the modified Bessel function I_1 (x) for any real x .
inline double BESSI1(double x)
{
  double ax = fabs(x), ans = R::bessel_i(ax, 1, 1);
  return x < 0.0 ? -ans : ans;
}

inline double BESSI(int n, double x)
{
  double ax = fabs(x), ans = R::bessel_i(ax, n, 1);

  return x < 0.0 && (n & 1) ? -ans : ans;
}

// // Returns the modified Bessel function I_0 (x) for any real x
// inline double BESSI0(double x)
// {
//   double ax,ans=0;
//   double y;
//   if ((ax=fabs(x)) < 3.75) { // Polynomial fit
//     y=x/3.75;
//     y*=y;
//     ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
//                                            +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
//   } else {
//     y=3.75/ax;
//     ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
//                                             +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
//                                                                                  +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
//                                                                                  +y*0.392377e-2))))))));
//   }
//   return ans;
// }
//
//
// // Returns the modified Bessel function I_1 (x) for any real x .
// inline double BESSI1(double x)
// {
//   double ax,ans;
//   double y;
//   if ((ax=fabs(x)) < 3.75) { //Polynomial fit.
//     y=x/3.75;
//     y*=y;
//     ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
//                                                  +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
//   } else {
//     y=3.75/ax;
//     ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
//                                            -y*0.420059e-2));
//     ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
//                                          +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
//     ans *= (exp(ax)/sqrt(ax));
//   }
//   return x < 0.0 ? -ans : ans;
// }
//
// //   /*----------------------------------------------------------------------
// //   !     This subroutine calculates the first kind modified Bessel function
// //   !     of integer order n, for any REAL x. The codes are adapted from
// //   !     Numerical recipes in C (Vol. 2) by Press, W. H., Teukolsky, S. A.,
// //   !     Vetterling, W. T., & Flannery, B. P. (1996), Cambridge university
// //   !     press, p. 236-240
// //   ------------------------------------------------------------------------*/
// inline double BESSI(int n, double x)
// {
//   int j;
//   double ACC = 40.0;  //Make larger to increase accuracy.
//   double BIGNO = 1.0e10, BIGNI = 1.0e-10;
//   double bi,bim,bip,tox,ans;
//   if (n == 0) return BESSI0(x);
//   else if (n == 1) return BESSI1(x);
//   else if (x == 0.0) return 0.0;
//   else {
//     tox=2.0/fabs(x);
//     bip=ans=0.0;
//     bi=1.0;
//     for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) { // Downward recurrence from even m
//       bim=bip+j*tox*bi;
//       bip=bi;
//       bi=bim;
//       if (fabs(bi) > BIGNO) {  // Renormalize to prevent overflows.
//         ans *= BIGNI;
//         bi *= BIGNI;
//         bip *= BIGNI;
//       }
//       if (j == n) ans=bip;
//     }
//     ans *= BESSI0(x)/bi; // Normalize with BESSI0.
//     return x < 0.0 && (n & 1) ? -ans : ans;
//   }
// }



inline int sgn(double val) {
  if (val > 0) {
    return 1;
  } else if (val < 0) {
    return -1;
  } else {
    return 0;
  }
}

inline double plus(double x) {
  return fmax(x, 1e-10);
}


#endif


