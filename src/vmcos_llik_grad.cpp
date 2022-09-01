#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
//
// #include "thread_num.h"

#include "bessel.h"
#include "univmgen.h"

using namespace Rcpp;

// typedef std::vector<arma::mat> stdvec_mat;
typedef std::vector<arma::vec> stdvec_vec;



// Returns the modified Bessel function I_0 (x) for any real x
// same as the original function in bessel.h, but this one exports to R
// [[Rcpp::export]]
double BESSI0_C(double x)
{
  // double ax,ans=0;
  // double y;
  // if ((ax=fabs(x)) < 3.75) { // Polynomial fit
  //   y=x/3.75;
  //   y*=y;
  //   ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
  //                                          +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  // } else {
  //   y=3.75/ax;
  //   ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
  //                                           +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
  //                                                                                +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
  //                                                                                +y*0.392377e-2))))))));
  // }
  // return ans;
  return R::bessel_i(fabs(x), 0, 1);
}



// [[Rcpp::export]]
double const_vmcos_anltc(double k1, double k2, double k3)
{
  int p = 1;
  double start, temp;
  start = BESSI0(k1) * BESSI0(k2) * BESSI0(k3);
  temp = 2 * BESSI(p, k1) * BESSI(p, k2) * BESSI(p, k3);
  start += temp;
  while((fabs(temp) > 1e-7) && (p <= 10000)) {
    p++;
    temp = 2 * BESSI(p, k1) * BESSI(p, k2) * BESSI(p, k3);
    start += temp;
  }
  return  4 * M_PI * M_PI * start;
}


// [[Rcpp::export]]
double const_vmcos_mc(double k1, double k2, double k3, arma::mat uni_rand,
                      bool return_log = false)
{
  double x = uni_rand(0, 0) * 2 * M_PI,
    y = uni_rand(0, 1) * 2 * M_PI,
    temp = k1 * cos(x) + k2 * cos(y) + k3 * cos(x - y),
    expon, sum_exp = 1.0;
  int nsim = uni_rand.n_rows, i;

  for(i = 1; i < nsim; i++) {
    x = uni_rand(i, 0) * 2 * M_PI;
    y = uni_rand(i, 1) * 2 * M_PI;
    expon = k1 * cos(x) + k2 * cos(y) + k3 * cos(x - y);
    sum_exp += exp(expon - temp);
  }

  double out;
  if (return_log) {
    out = temp + log(4 * M_PI * M_PI * sum_exp / nsim);
  } else {
    out = 4 * M_PI * M_PI * exp(temp) * sum_exp / nsim;
  }

  return out;
}


// [[Rcpp::export]]
double const_vmcos(double k1, double k2, double k3, arma::mat  uni_rand,
                   bool return_log = false)
{
  if((k3 >= -5 && fmin(k1, k2) >= 0.01 &&
     fmax(fmax(k1, k2), fabs(k3)) <= 25) ||
     fabs(k3) < 1e-4) {
    // if(k3 >= 0) {
    double out = const_vmcos_anltc(k1, k2, k3);
    if (return_log) {
      return log(out);
    } else {
      return out;
    }
  } else {
    return const_vmcos_mc(k1, k2, k3, uni_rand, return_log);
  }
}


// // not used
// // [[Rcpp::export]]
// double d_const_vmcos_k1_anltc(double k1, double k2, double k3)
// {
//   int p = 1;
//   double start, temp;
//   start = BESSI1(k1) * BESSI0(k2) * BESSI0(k3);
//   temp = (BESSI(p-1, k1) + BESSI(p+1, k1))  * BESSI(p, k2) * BESSI(p, k3);
//   start += temp;
//   while((fabs(temp) > 1e-7) && (p <= 10000)) {
//     p++;
//     temp = (BESSI(p-1, k1) + BESSI(p+1, k1))  * BESSI(p, k2) * BESSI(p, k3);
//     start += temp;
//   }
//   return  4 * M_PI * M_PI * start;
// }


// // [[Rcpp::export]]
// arma::vec d_const_vmcos_anltc(double k1, double k2, double k3)
// {
//   arma::vec result;
//   result << d_const_vmcos_k1_anltc(k1, k2, k3)
//          << d_const_vmcos_k1_anltc(k2, k3, k1)
//          << d_const_vmcos_k1_anltc(k3, k1, k2)
//          << arma::endr;
//   return result;
// }


// [[Rcpp::export]]
arma::vec d_const_vmcos_anltc(double k1, double k2, double k3)
{
  int p = 1;
  double
    start1 = BESSI1(k1) * BESSI0(k2) * BESSI0(k3),
      start2 = BESSI1(k2) * BESSI0(k1) * BESSI0(k3),
      start3 = BESSI1(k3) * BESSI0(k2) * BESSI0(k1),
      I_p_k1 = BESSI(p, k1), I_p_k2 = BESSI(p, k2), I_p_k3 = BESSI(p, k3),
      I_p_pls1_k1 = BESSI(p+1, k1), I_p_pls1_k2 = BESSI(p+1, k2), I_p_pls1_k3 = BESSI(p+1, k3),
      I_p_mins1_k1 = BESSI(p-1, k1), I_p_mins1_k2 = BESSI(p-1, k2), I_p_mins1_k3 = BESSI(p-1, k3),

      temp1 = (I_p_pls1_k1 + I_p_mins1_k1)  * I_p_k2 * I_p_k3,
      temp2 = (I_p_pls1_k2 + I_p_mins1_k2)  * I_p_k1 * I_p_k3,
      temp3 = (I_p_pls1_k3 + I_p_mins1_k3)  * I_p_k1 * I_p_k2;

  start1 += temp1;
  start2 += temp2;
  start3 += temp3;

  while((std::max(fabs(temp1), std::max(fabs(temp2), fabs(temp3))) > 1e-7) && (p <= 10000)) {
    p++;
    // shift p on right by 1
    I_p_mins1_k1 = I_p_k1;
    I_p_mins1_k2 = I_p_k2;
    I_p_mins1_k3 = I_p_k3;
    I_p_k1 = I_p_pls1_k1;
    I_p_k2 = I_p_pls1_k2;
    I_p_k3 = I_p_pls1_k3;
    I_p_pls1_k1 = BESSI(p+1, k1);
    I_p_pls1_k2 = BESSI(p+1, k2);
    I_p_pls1_k3 = BESSI(p+1, k3);

    temp1 = (I_p_pls1_k1 + I_p_mins1_k1)  * I_p_k2 * I_p_k3;
    temp2 = (I_p_pls1_k2 + I_p_mins1_k2)  * I_p_k1 * I_p_k3;
    temp3 = (I_p_pls1_k3 + I_p_mins1_k3)  * I_p_k1 * I_p_k2;

    start1 += temp1;
    start2 += temp2;
    start3 += temp3;
  }
  arma::vec out = {
    4 * M_PI * M_PI * start1,
    4 * M_PI * M_PI * start2,
    4 * M_PI * M_PI * start3
  };

  return  out;
}


// [[Rcpp::export]]
arma::vec d_const_vmcos_mc(double k1, double k2, double k3, arma::mat uni_rand, int ncores = 1)
{
  double x = uni_rand(0, 0) * 2 * M_PI,
    y = uni_rand(0, 1) * 2 * M_PI,
    expon, sum_exp1 = cos(x), sum_exp2 = cos(y), sum_exp3 = cos(x-y),
    temp = k1 * sum_exp1 + k2 * sum_exp2 + k3 * sum_exp3;

  int nsim = uni_rand.n_rows, i;

  for(i = 1; i < nsim; i++) {
    x = uni_rand(i, 0) * 2 * M_PI;
    y = uni_rand(i, 1) * 2 * M_PI;
    expon = k1 * cos(x) + k2 * cos(y) + k3 * cos(x - y);
    sum_exp1 += exp(expon - temp) * cos(x);
    sum_exp2 += exp(expon - temp) * cos(y);
    sum_exp3 += exp(expon - temp) * cos(x-y);
  }
  arma::vec result = {
    4 * M_PI * M_PI * exp(temp) * sum_exp1 / nsim,
    4 * M_PI * M_PI * exp(temp) * sum_exp2 / nsim,
    4 * M_PI * M_PI * exp(temp) * sum_exp3 / nsim
  };


  return result;
}


// [[Rcpp::export]]
arma::vec d_const_vmcos(arma::vec par, arma::mat uni_rand, int ncores = 1)
{
  double k1 = par[0], k2 = par[1], k3 = par[2];
  if ((k3 >= -5 && fmin(k1, k2) >= 0.01 &&
      fmax(fmax(k1, k2), fabs(k3)) <= 25) ||
      fabs(k3) < 1e-4) {
    // if(k3 >= 0) {
    return d_const_vmcos_anltc(k1, k2, k3);
  } else {
    return d_const_vmcos_mc(k1, k2, k3, uni_rand, ncores);
  }
}


// [[Rcpp::export]]
arma::vec log_const_vmcos_all(arma::mat par_mat, arma::mat uni_rand) {
  int K = par_mat.n_cols;
  arma::vec all_lconsts(K);
  for(int j = 0; j < K; j++)
    all_lconsts[j] = const_vmcos(par_mat(0,j), par_mat(1,j), par_mat(2,j), uni_rand, true);
  return all_lconsts;
}





// [[Rcpp::export]]
double ldcosnum(double x, double y, arma::vec par)
{
  double k1 = par[0], k2 = par[1], k3 = par[2], mu = par[3], nu = par[4];
  if(k1 <= 0 || k2 <= 0 )
    return NA_REAL; else
      return (k1*cos(x-mu) + k2*cos(y-nu) + k3 * cos(x-mu-y+nu));
}

// [[Rcpp::export]]
arma::vec dcos_onex_manypar(arma::vec x, arma::vec k1, arma::vec k2, arma::vec k3,
                            arma::vec mu1, arma::vec mu2, arma::vec l_const_all)
{
  int n = k1.size();

  arma::mat all_par(5, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k1[i];
    all_par(1,i) = k2[i];
    all_par(2,i) = k3[i];
    all_par(3,i) = mu1[i];
    all_par(4,i) = mu2[i];
  }

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldcosnum(x[0], x[1], all_par.col(i));
  }

  return arma::exp(ld_num - l_const_all);
}

// [[Rcpp::export]]
arma::vec dcos_manyx_onepar(arma::mat x, double k1, double k2, double k3,
                            double mu1, double mu2, double l_const)
{
  int n = x.n_rows;

  // arma::vec par(5);
  // par[0] = k1; par[1] = k2; par[2] = k3; par[3] = mu1; par[4] = mu2;
  arma::vec par = { k1, k2, k3, mu1, mu2 };



  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldcosnum(x(i,0), x(i,1), par);
  }

  return arma::exp(ld_num - l_const);
}


// [[Rcpp::export]]
arma::vec dcos_manyx_manypar(arma::mat x, arma::vec k1, arma::vec k2, arma::vec k3,
                             arma::vec mu1, arma::vec mu2, arma::vec l_const_all)
{
  int n = k1.size();

  arma::mat all_par(5, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k1[i];
    all_par(1,i) = k2[i];
    all_par(2,i) = k3[i];
    all_par(3,i) = mu1[i];
    all_par(4,i) = mu2[i];
  }

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldcosnum(x(i,0), x(i,1), all_par.col(i));
  }

  return arma::exp(ld_num - l_const_all);
}



// generates a single observation
arma::rowvec2 rcos_unimodal_single(double k1, double k2, double k3,
                                   double mu1, double mu2,
                                   double kappa_opt, double log_I0_kappa_opt,
                                   double logK,
                                   double log_const_vmcos) {
  double x, y, U1;


  double log_prop_den, log_target_den,
  log_2pi=log(2*M_PI), cos_y_mu2;

  // do marginal y then conditional x|y generation


  // generating from marginal of y:
  // for marginal of y, proposal is vm(kappa_optim, mu2)
  double k13;
  int accpt_y = 0, ntry = 1; // will be 1 if accepted in rejection sampling from y marginal
  // First generate y from the marginal distribution given in Eqn. 5 in Mardia 2007
  while (accpt_y < 1) {
    y = runivm_single_onepar(kappa_opt, mu2);
    cos_y_mu2 = cos(y-mu2);
    k13 = sqrt(k1*k1 + k3*k3 + 2*k1*k3*cos_y_mu2);
    U1 = R::unif_rand();
    log_target_den = -log_const_vmcos+log_2pi+log(BESSI0(k13))+k2*cos_y_mu2;
    log_prop_den = -log_2pi - log_I0_kappa_opt + kappa_opt*cos_y_mu2;
    // Rcout << ntry << "\t";
    ntry++;
    if (log(U1) <= log_target_den-log_prop_den-logK)
      accpt_y++;
  }


  double mean_x = mu1+atan(k3*sin(y-mu2)/(k1 + k3*cos_y_mu2));

  // generating x|y
  // generate x from vM(k13, mean_x)
  x = runivm_single_onepar(k13, mean_x);

  arma::rowvec2 final_sample = { x, y };

  return final_sample;
}



// [[Rcpp::export]]
arma::mat rcos_unimodal(int n, double k1, double k2, double k3,
                        double mu1, double mu2,
                        double kappa_opt, double log_I0_kappa_opt,
                        double logK,
                        double log_const_vmcos)
{
  arma::mat out(n, 2);
  for (int i = 0; i < n; i++)
    out.row(i) =
      rcos_unimodal_single(k1, k2, k3, mu1, mu2, kappa_opt,
                           log_I0_kappa_opt, logK, log_const_vmcos);

  return out;
}



arma::rowvec2 rcos_bimodal_single(double k1, double k2, double k3,
                                  double mu1, double mu2,
                                  double kappa_opt, double log_I0_kappa_opt,
                                  double logK,
                                  double log_const_vmcos,
                                  double mode_1, double mode_2,
                                  double vmpropn, double unifpropn)

{
  double x, y, U1, U0;


  double log_prop_den, log_target_den,
  log_2pi=log(2*M_PI), cos_y_mu2;

  // do marginal y then conditional x|y generation


  // generating from marginal of y:
  // for marginal of y, proposal is vmpropn-vmpropn-unifpropn
  // mixture of vm(kappa_optim, mode_1), vm(kappa_optim, mode_2)
  // unif(0, 2*pi)
  double k13;

  int accpt_y = 0, ntry = 1; // will be 1 if accepted in rejection sampling from y marginal
  // First generate y from the marginal distribution given in Eqn. 5 in Mardia 2007
  while (accpt_y < 1) {
    U0 = R::unif_rand();

    if (U0 < vmpropn) {
      y = runivm_single_onepar(kappa_opt, mode_1);
    } else if (U0 < 2*vmpropn) {
      y = runivm_single_onepar(kappa_opt, mode_2);
    } else {
      y = R::runif(0, 2*M_PI);
    }
    cos_y_mu2 = cos(y-mu2);
    k13 = sqrt(k1*k1 + k3*k3 + 2*k1*k3*cos_y_mu2);
    log_target_den = -log_const_vmcos
      +log_2pi+log(BESSI0(k13))+k2*cos_y_mu2;
    log_prop_den =
    -log_2pi +  log(exp(log(vmpropn) - log(BESSI0(kappa_opt)) +
    kappa_opt*cos(y-mode_1) +
    log(1+exp(kappa_opt*(cos(y-mode_2) - cos(y-mode_1)))))
                      + unifpropn);
    // Rcout << ntry << "\t";
    ntry++;
    if (ntry % 100 == 0) checkUserInterrupt();
    U1 = R::unif_rand();
    if (log(U1) <= log_target_den-log_prop_den-logK)
      accpt_y++;
  }



  double mean_x = mu1+atan(k3*sin(y-mu2)/(k1 + k3*cos_y_mu2));

  // generating x|y
  // generate x from vM(k13, mean_x)
  x = runivm_single_onepar(k13, mean_x);

  arma::rowvec2 final_sample = { x, y };

  return final_sample;
}


//  [[Rcpp::export]]
arma::mat rcos_bimodal(int n, double k1, double k2, double k3,
                       double mu1, double mu2,
                       double kappa_opt, double log_I0_kappa_opt,
                       double logK,
                       double log_const_vmcos,
                       double mode_1, double mode_2,
                       double vmpropn, double unifpropn)

{
  arma::mat out(n, 2);
  for (int i = 0; i < n; i++)
    out.row(i) =
      rcos_bimodal_single(k1, k2, k3, mu1, mu2, kappa_opt,
                          log_I0_kappa_opt, logK, log_const_vmcos,
                          mode_1, mode_2, vmpropn, unifpropn);

  return out;
}


// generates a single observation using naive 2d rejection
arma::rowvec rcos_single_onepar(double k1, double k2, double k3,
                                double mu1, double mu2, double I_upper_bd) {
  arma::vec par(5);
  double x, y, U1;

  arma::rowvec final_sample(2);

  // brute force bivariate rejection sampling
  // with numerically computed upper bound I_upper_bd (passed from R)
  int accpt = 0;
  while (accpt < 1) {
    x = R::runif(0, 2*M_PI);
    y = R::runif(0, 2*M_PI);
    U1 = R::unif_rand();
    if (log(U1) <= (k1*cos(x-mu1)+k2*cos(y-mu2)+k3*cos(x-y-mu1+mu2))-I_upper_bd) {
      accpt++;
    }
  }

  final_sample[0] = x;
  final_sample[1] = y;

  return final_sample;
}



//  [[Rcpp::export]]
arma::mat rcos_onepar(int n, double k1, double k2, double k3,
                      double mu1, double mu2, double I_upper_bd) {
  if(n == 1){
    return(rcos_single_onepar(k1, k2, k3, mu1, mu2, I_upper_bd));
  } else {
    arma::mat all_sample(n, 2);
    for(int i = 0; i < n; i++)
      all_sample.row(i) = rcos_single_onepar(k1, k2, k3, mu1, mu2, I_upper_bd);
    return(all_sample);
  }
}


// // not needed
// //  [[Rcpp::export]]
// arma::mat rcos_manypar(arma::vec k1, arma::vec k2, arma::vec k3,
//                        arma::vec mu1, arma::vec mu2, arma::vec I_upper_bd) {
//   int n = k1.size();
//   arma::mat all_sample(n, 2);
//   for(int i = 0; i < n; i++)
//     all_sample.row(i) = rcos_single_onepar(k1[i], k2[i], k3[i], mu1[i], mu2[i], I_upper_bd[i]);
//   return(all_sample);
// }


// [[Rcpp::export]]
arma::mat mem_p_cos(arma::mat data, arma::mat par, arma::vec pi,
                    arma::vec log_c_von)
{
  int n = data.n_rows, K = par.n_cols, j;
  double row_total;
  arma::mat den(n, K);
  for(int i = 0; i < n; i++){
    row_total = 0;
    for(j = 0; j < K; j++){
      den(i, j) = pi[j]*exp(ldcosnum(data(i, 0), data(i, 1), par.col(j)) - log_c_von[j]);;
      row_total += den(i, j);
    }
    if(row_total < 1e-50)
      row_total = 1e-50;
    for(j = 0; j < K; j++)
      den(i, j) /= row_total;
  }
  return(den);
}

// // not used
// // [[Rcpp::export]]
// double llik_vmcos_full(arma::mat data, arma::mat par,
//                        arma::vec pi, arma::vec log_c, int ncores = 1)
// {
//
//   int n = data.n_rows, K = pi.size(), j;
//   double temp, log_sum = 0.0;
//   arma::vec log_pi = log(pi);
//   if(K > 1) {
//     for(int i = 0; i < n; i++) {
//       temp = 0;
//       for(j = 0; j < K; j++)
//         temp += exp(ldcosnum(data(i,0), data(i,1), par.col(j)) - log_c[j] + log_pi[j] );
//       log_sum += log(temp);
//     }
//   } else {
//     for(int i = 0; i < n; i++)
//       log_sum += ldcosnum(data(i,0), data(i,1), par);
//     log_sum -= n*log_c[0];
//   }
//
//   return(log_sum);
//
// }



// [[Rcpp::export]]
arma::vec llik_vmcos_contri_C(arma::mat data, arma::mat par, arma::vec pi, arma::vec log_c)
{

  int n = data.n_rows, K = pi.size() ,j;
  double temp;
  arma::vec llik_contri = arma::vec(n);
  arma::vec log_pi = log(pi);

  if(K > 1) {
    for(int i = 0; i < n; i++) {
      temp = 0;
      for(j = 0; j < K; j++)
        temp += exp(ldcosnum(data(i,0), data(i,1), par.col(j)) - log_c[j] + log_pi[j] );
      llik_contri[i] = log(temp);
    }
  }
  else {
    for(int i = 0; i < n; i++) {
      llik_contri[i] = ldcosnum(data(i,0), data(i,1), par) - log_c[0];
    }
  }
  return(llik_contri);

}



// [[Rcpp::export]]
double llik_vmcos_one_comp(arma::mat data, arma::vec par_vec,
                           double log_c)
{

  int n = data.n_rows;
  double log_sum = 0.0;

  for(int i = 0; i < n; i++)
    log_sum += ldcosnum(data(i,0), data(i,1), par_vec);
  log_sum -= n*log_c;
  return (log_sum);

}


// // not used
// // [[Rcpp::export]]
// arma::vec grad_log_vmcos_one_comp_i(double x, double y, arma::vec par,
//                                     double c_vmcos, arma::vec del_const_vmcos)
// {
//   double k1 = par[0], k2 = par[1], k3 = par[2], mu1 = par[3], mu2 = par[4];
//   arma::vec all_entries(5);
//   all_entries[0] = (cos(x - mu1) - del_const_vmcos[0]/c_vmcos);
//   all_entries[1] = (cos(y - mu2) - del_const_vmcos[1]/c_vmcos);
//   all_entries[2] = (cos(x - mu1 - y + mu2) - del_const_vmcos[2]/c_vmcos);
//   all_entries[3] = (k1*sin(x - mu1) + k3*sin(x - mu1 - y + mu2));
//   all_entries[4] = (k2*sin(y - mu2) - k3*sin(x - mu1 - y + mu2));
//
//   return all_entries;
// }
//
//
// // not used
// // [[Rcpp::export]]
// arma::mat grad_vmcos_all_comp(arma::mat data, arma::mat par_mat,
//                               arma::vec pi, arma::mat uni_rand, int ncores = 1)
// {
//   int n = data.n_rows, K = pi.size(), j;
//   arma::vec l_c_vmcos = log_const_vmcos_all(par_mat, uni_rand, ncores);
//   arma::vec c_vmcos = arma::exp(l_c_vmcos), log_pi = log(pi);
//
//   stdvec_vec del_const_vmcos_all(K);
//   stdvec_vec par_mat_cols(K);
//
//   for(j = 0; j < K; j++) {
//     del_const_vmcos_all[j] = d_const_vmcos(par_mat.col(j), uni_rand, ncores);
//     par_mat_cols[j] = par_mat.col(j);
//   }
//
//   stdvec_mat grad_sum_part(ncores);
//   for(int g = 0; g < ncores; g++)
//     grad_sum_part[g] = arma::zeros(5, K);
//
//   double denom, temp;
//
//   for(int i = 0; i < n; i++) {
//     arma::mat grad_temp(5, K);
//     int g = get_thread_num_final();
//     // int g = 0;
//     denom = 0;
//     for(j = 0; j < K; j++) {
//       temp = exp(ldcosnum(data(i, 0), data(i, 1), par_mat_cols[j]) - l_c_vmcos[j] + log_pi[j]);
//       denom += temp;
//       grad_temp.col(j) = temp*grad_log_vmcos_one_comp_i(data(i, 0), data(i, 1), par_mat_cols[j], c_vmcos[j], del_const_vmcos_all[j]);
//     }
//     grad_sum_part[g] += grad_temp/denom;
//   }
//
//   arma::mat grad_sum = arma::zeros(5, K);
//   for(int g = 0; g < ncores; g++)
//     grad_sum += grad_sum_part[g];
//
//   return (grad_sum);
// }


// [[Rcpp::export]]
arma::vec grad_llik_vmcos_C(arma::mat data,
                            arma::vec par,
                            arma::mat uni_rand) {
  int n = data.n_rows;
  double k1 = par[0], k2 = par[1], k3 = par[2], mu1 = par[3], mu2 = par[4];


  arma::vec grad_llik_wo_const = arma::zeros(6);
  // first 5 entries = grad, last entry = llik

  for(int i = 0; i < n; i++) {
    double x_mu1 = (data(i, 0) - mu1), y_mu2 = (data(i, 1) - mu2),
      x_mu1_y_mu2 = x_mu1 - y_mu2;
    double  cos_x_mu1 = cos(x_mu1), sin_x_mu1 = sin(x_mu1),
      cos_y_mu2 = cos(y_mu2), sin_y_mu2 = sin(y_mu2),
      cos_x_mu1_y_mu2 = cos(x_mu1_y_mu2), sin_x_mu1_y_mu2 = sin(x_mu1_y_mu2);

    grad_llik_wo_const[0] += cos_x_mu1; //(cos(x - mu1) - del_const_vmcos[0]/c_vmcos);
    grad_llik_wo_const[1] += cos_y_mu2; // - del_const_vmcos[1]/c_vmcos);
    grad_llik_wo_const[2] += cos_x_mu1_y_mu2; // - del_const_vmcos[2]/c_vmcos);
    grad_llik_wo_const[3] += k1*sin_x_mu1 + k3*sin_x_mu1_y_mu2;
    grad_llik_wo_const[4] += k2*sin_y_mu2 - k3*sin_x_mu1_y_mu2;
    grad_llik_wo_const[5] += k1*cos_x_mu1 + k2*cos_y_mu2 + k3*cos_x_mu1_y_mu2;
  }

  double c_vmcos = const_vmcos(k1, k2, k3, uni_rand);
  arma::vec del_const_vmcos = d_const_vmcos(par, uni_rand);

  // adjust const terms in the first three pars
  for(int j = 0; j < 3; j++) {
    grad_llik_wo_const[j] -= n * del_const_vmcos[j]/c_vmcos;
  }

  // subtract constant term from llik
  grad_llik_wo_const[5] -= n * log(c_vmcos);

  return(grad_llik_wo_const);
}



// // not used
// // [[Rcpp::export]]
// double vmcosmix(double x, double y, arma::mat par, arma::vec pi, arma::vec log_c_von)
// {
//   double res = 0;
//   int K = par.n_cols;
//
//   for(int j = 0; j < K; j++)
//     res += exp(ldcosnum(x, y, par.col(j)) - log_c_von[j]) * pi[j];
//   return(res);
// } //likelihood contribution (density) of a single point in a mixture model
//
// // not used
// // [[Rcpp::export]]
// arma::vec vmcosmix_manyx(arma::mat x, arma::mat par, arma::vec pi, arma::vec log_c_von)
// {
//   int n = x.n_rows;
//   arma::vec result(n);
//   for(int i = 0; i < n; i++)
//     result[i] = vmcosmix(x(i,0), x(i,1), par, pi, log_c_von);
//   return(result);
// }



// [[Rcpp::export]]
List vmcos_var_corr_anltc(double k1, double k2, double k3)
{

  // double c = 0, c_k1 = 0, c_k2 = 0, c_k3 = 0, c_k1k2 = 0,
  //   c_k1k1 = 0, c_k2k2 = 0;

  int p = 0;

  double
    I_p_k1 = BESSI(p, k1),
      I_p_k2 = BESSI(p, k2),
      I_p_k3 = BESSI(p, k3),
      I_pplus1_k1 = BESSI(p+1, k1),
      I_pplus1_k2 = BESSI(p+1, k2),
      I_pplus1_k3 = BESSI(p+1, k3),
      I_pplus2_k1 = BESSI(p+2, k1),
      I_pplus2_k2 = BESSI(p+2, k2),
      I_pmins1_k1 = I_pplus1_k1,
      I_pmins1_k2 = I_pplus1_k2,
      I_pmins1_k3 = I_pplus1_k3,
      I_pmins2_k1 = I_pplus2_k1,
      I_pmins2_k2 = I_pplus2_k2;

  double
    c = 0.5 * I_p_k1 * I_p_k2 * I_p_k3,
      c_k1 = I_pplus1_k1 * I_p_k2 * I_p_k3,
      c_k2 = I_pplus1_k2 * I_p_k1 * I_p_k3,
      c_k3 = I_pplus1_k3 * I_p_k2 * I_p_k1,
      c_k1k1 = (I_p_k1 + I_pplus2_k1) * I_p_k2 * I_p_k3,
      c_k2k2 = (I_p_k2 + I_pplus2_k2) * I_p_k1 * I_p_k3,
      c_k1k2 = 2 * I_pplus1_k1 * I_pplus1_k2 * I_p_k3;


  double den_min = std::min(c, std::min(k1, k2)),
    temp_c = den_min;

  while(temp_c/den_min > 1e-6) {
    p++;

    I_pmins2_k1 = I_pmins1_k1;
    I_pmins2_k2 = I_pmins1_k2;
    I_pmins1_k1 = I_p_k1;
    I_pmins1_k2 = I_p_k2;
    I_pmins1_k3 = I_p_k3;
    I_p_k1 = I_pplus1_k1;
    I_p_k2 = I_pplus1_k2;
    I_p_k3 = I_pplus1_k3;
    I_pplus1_k1 = I_pplus2_k1;
    I_pplus1_k2 = I_pplus2_k2;
    I_pplus1_k3 = BESSI(p+1, k3);
    I_pplus2_k1 = BESSI(p+2, k1);
    I_pplus2_k2 = BESSI(p+2, k2);

    temp_c = I_p_k1 * I_p_k2 * I_p_k3;
    c += temp_c;
    c_k1 += (I_pplus1_k1 + I_pmins1_k1)  * I_p_k2 * I_p_k3;
    c_k2 += (I_pplus1_k2 + I_pmins1_k2)  * I_p_k1 * I_p_k3;
    c_k3 += (I_pplus1_k3 + I_pmins1_k3)  * I_p_k2 * I_p_k1;
    c_k1k2 += (I_pplus1_k1 + I_pmins1_k1) * (I_pplus1_k2 + I_pmins1_k2) * I_p_k3;
    c_k1k1 += (I_pmins2_k1 + 2*I_p_k1 + I_pplus2_k1)  * I_p_k2 * I_p_k3;
    c_k2k2 += (I_pmins2_k2 + 2*I_p_k2 + I_pplus2_k2)  * I_p_k1 * I_p_k3;

  }


  double pi_sq_times_4 = 4 * M_PI * M_PI;

  c *= pi_sq_times_4 * 2;
  c_k1 *= pi_sq_times_4;
  c_k2 *= pi_sq_times_4;
  c_k3 *= pi_sq_times_4;
  c_k1k2 *= pi_sq_times_4 * 0.5;
  c_k1k1 *= pi_sq_times_4 * 0.5;
  c_k2k2 *= pi_sq_times_4 * 0.5;

  double
    rho_js = (fabs(c_k3 - c_k1k2) < 1e-10) ? 0 :
    sgn(c_k3 - c_k1k2) *
      fmin(1, exp(log(plus(fabs(c_k3 - c_k1k2)))
                    - 0.5*log(plus(c-c_k1k1))
                    - 0.5*log(plus(c-c_k2k2)))
      ),


      // rho_fl = rho_js * c_k1k2 / sqrt(c_k1k1 * c_k2k2),
      rho_fl = (fabs(c_k1k2) < 1e-10) ?  0 :
    rho_js * sgn(c_k1k2) *
      fmin(1, exp(log(plus(fabs(c_k1k2)))
                    - 0.5*log(plus(c_k1k1))
                    - 0.5*log(plus(c_k2k2)))
      );

  // double
  //   rho_js = (c_k3 - c_k1k2)/sqrt((c-c_k1k1) * (c-c_k2k2)),
  //     rho_fl = rho_js * c_k1k2 / sqrt(c_k1k1 * c_k2k2);

  double
    // var1 = 1 - (c_k1/c);
    var1 = fmin(1,
                1 - sgn(c_k1) *
                  exp(log(plus(fabs(c_k1)))
                        - log(plus(c)))
    ),
    // var2 = 1 - (c_k2/c);
    var2 = fmin(1, 1 - sgn(c_k2)*exp(log(plus(fabs(c_k2)))
                                       - log(plus(c)))
    );


  return List::create(Rcpp::Named("var1") = var1,
                      Rcpp::Named("var2") = var2,
                      Rcpp::Named("rho_fl") = rho_fl,
                      Rcpp::Named("rho_js") = rho_js);

}



// [[Rcpp::export]]
List vmcos_var_corr_mc(double k1, double k2, double k3,
                       arma::mat uni_rand, int ncores = 1)
{
  double x, y, expon, cos_x, cos_y, cos_x_y,
  sum_exp, sum_exp1, sum_exp2, sum_exp3, sum_exp11,
  sum_exp12, sum_exp22,temp, // exp_temp,
  exp_expon_temp;

  int nsim = uni_rand.n_rows;

  x = uni_rand(0, 0) * 2 * M_PI;
  y = uni_rand(0, 1) * 2 * M_PI;
  cos_x = cos(x);
  cos_y = cos(y);
  cos_x_y = cos(x-y);
  sum_exp = 1.0;
  sum_exp1 = cos_x;
  sum_exp2 = cos_y;
  sum_exp3 = cos_x_y;
  sum_exp11 = cos_x*cos_x;
  sum_exp22 = cos_y*cos_y;
  sum_exp12 = cos_x*cos_y;
  temp = k1 * cos_x + k2 * cos_y + k3 * cos_x_y;
  // exp_temp = exp(temp);

  for(int i = 1; i < nsim; i++) {
    x = uni_rand(i, 0) * 2 * M_PI;
    y = uni_rand(i, 1) * 2 * M_PI;
    cos_x = cos(x);
    cos_y = cos(y);
    cos_x_y = cos(x-y);
    expon = k1 * cos_x + k2 * cos_y + k3 * cos_x_y;
    exp_expon_temp = exp(expon - temp);

    sum_exp += exp_expon_temp;
    sum_exp1 += exp_expon_temp *  cos_x;
    sum_exp2 += exp_expon_temp *  cos_y;
    sum_exp3 += exp_expon_temp *  cos_x_y;
    sum_exp11 += exp_expon_temp *  cos_x*cos_x;
    sum_exp22 += exp_expon_temp *  cos_y*cos_y;
    sum_exp12 += exp_expon_temp *  cos_x*cos_y;

  }



  double pisq_4_exp_temp_over_nsim = 4 * M_PI * M_PI / nsim;

  // all of the following need to be multiplied by exp(expon)
  double c = sum_exp * pisq_4_exp_temp_over_nsim,
    c_k1 = sum_exp1 * pisq_4_exp_temp_over_nsim,
    c_k2 = sum_exp2 * pisq_4_exp_temp_over_nsim,
    c_k3 = sum_exp3 * pisq_4_exp_temp_over_nsim,
    c_k1k2 = sum_exp12 * pisq_4_exp_temp_over_nsim,
    c_k1k1 = sum_exp11 * pisq_4_exp_temp_over_nsim,
    c_k2k2 = sum_exp22 * pisq_4_exp_temp_over_nsim;


  double
    rho_js = (fabs(c_k3 - c_k1k2) < 1e-10) ? 0 :
    sgn(c_k3 - c_k1k2) *
      fmin(1, exp(log(plus(fabs(c_k3 - c_k1k2)))
                    - 0.5*log(plus(c-c_k1k1))
                    - 0.5*log(plus(c-c_k2k2)))
      ),


      // rho_fl = rho_js * c_k1k2 / sqrt(c_k1k1 * c_k2k2),
      rho_fl = (fabs(c_k1k2) < 1e-10) ?  0 :
    rho_js * sgn(c_k1k2) *
      fmin(1, exp(log(plus(fabs(c_k1k2)))
                    - 0.5*log(plus(c_k1k1))
                    - 0.5*log(plus(c_k2k2)))
      );

  // double
  //   rho_js = (c_k3 - c_k1k2)/sqrt((c-c_k1k1) * (c-c_k2k2)),
  //     rho_fl = rho_js * c_k1k2 / sqrt(c_k1k1 * c_k2k2);

  double
    // var1 = 1 - (c_k1/c);
    var1 = fmin(1,
                1 - sgn(c_k1) *
                  exp(log(plus(fabs(c_k1)))
                        - log(plus(c)))
    ),
    // var2 = 1 - (c_k2/c);
    var2 = fmin(1, 1 - sgn(c_k2)*exp(log(plus(fabs(c_k2)))
                                       - log(plus(c)))
    );

  return List::create(Rcpp::Named("var1") = var1,
                      Rcpp::Named("var2") = var2,
                      Rcpp::Named("rho_fl") = rho_fl,
                      Rcpp::Named("rho_js") = rho_js);
}




// [[Rcpp::export]]
List vmcos_var_cor_singlepar_cpp(double k1, double k2, double k3,
                                 arma::mat uni_rand, int ncores = 1)
{
  if(k3 >= 0 && k1 <= 50 && k2 <= 50 && k3 <= 50) {
    return vmcos_var_corr_anltc(k1, k2, k3);
  } else {
    return vmcos_var_corr_mc(k1, k2, k3, uni_rand, ncores);
  }

}



// [[Rcpp::export]]
arma::vec ldcos_onex_manypar(arma::vec x, arma::vec k1, arma::vec k2, arma::vec k3,
                             arma::vec mu1, arma::vec mu2, arma::vec l_const_all)
{
  int n = k1.size();

  arma::mat all_par(5, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k1[i];
    all_par(1,i) = k2[i];
    all_par(2,i) = k3[i];
    all_par(3,i) = mu1[i];
    all_par(4,i) = mu2[i];
  }

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldcosnum(x[0], x[1], all_par.col(i));
  }

  return (ld_num - l_const_all);
}

// [[Rcpp::export]]
arma::vec ldcos_manyx_onepar(arma::mat x, double k1, double k2, double k3,
                             double mu1, double mu2, double l_const)
{
  int n = x.n_rows;

  arma::vec par = { k1, k2, k3, mu1, mu2 };
  // par[0] = k1; par[1] = k2; par[2] = k3; par[3] = mu1; par[4] = mu2;

  arma::vec out(n);
  for(int i = 0; i < n; i++) {
    out[i] = ldcosnum(x(i,0), x(i,1), par) - l_const;
  }

  return (out);
}


// [[Rcpp::export]]
arma::vec ldcos_manyx_manypar(arma::mat x, arma::vec k1, arma::vec k2, arma::vec k3,
                              arma::vec mu1, arma::vec mu2, arma::vec l_const_all)
{
  int n = k1.size();

  arma::mat all_par(5, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k1[i];
    all_par(1,i) = k2[i];
    all_par(2,i) = k3[i];
    all_par(3,i) = mu1[i];
    all_par(4,i) = mu2[i];
  }

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldcosnum(x(i,0), x(i,1), all_par.col(i));
  }

  return (ld_num - l_const_all);
}
