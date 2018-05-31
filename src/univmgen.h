#ifndef UNIVMGEN_H
#define UNIVMGEN_H

double prncp_reg_C(double x) ;
double runivm_single_onepar(double k, double mu);

// converts an angle into [0, 2pi]
inline double prncp_reg_C(double x)
{
  return x - 2 * M_PI * floor(x / (2 * M_PI));
}

// using the algorithm provided on p. 43, Mardia,
// Jupp

inline double runivm_single_onepar(double k, double mu)
{


  if (k > 1) {
    double a = 1 + sqrt(1 + 4 * k * k),
      b = (a - sqrt(2 * a))/(2 * k),
      r = (1 + b * b)/(2 * b);
    double U1, U2, U3;
    int check = 0;
    double z, f, c;
    double res;


    while (check < 1) {
      U1 = R::unif_rand();
      z = cos(M_PI * U1);
      f = (1. + r * z)/(r + z);
      c = k * (r - f);

      U2 = R::unif_rand();
      if (c * (2 - c) - U2 > 0) {
        U3 = unif_rand();
        if (U3 > 0.50) res = acos(f) + mu;
        else res = -acos(f) + mu;
        check++;
      } else if (log(c / U2) + 1 - c >= 0) {
        U3 = unif_rand();
        if (U3 > 0.50) res = acos(f) + mu;
        else res = -acos(f) + mu;
        check++;
      }
    }
    // now convert the final result into [0, 2 * M_PI]
    return prncp_reg_C(res);


  }

  else {
    int accpt = 0;
    double U, x;

    while (accpt < 1) {
      x = R::unif_rand()*2*M_PI;
      U = R::unif_rand();
      if (log(U) <= (k*cos(x-mu)-k)) {
        accpt++;
      }
    }

    return x;
  }

}

#endif
