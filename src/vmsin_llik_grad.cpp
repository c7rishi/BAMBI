#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
//
// #include "thread_num.h"
#include "bessel.h"
#include "univmgen.h"

// typedef std::vector<arma::mat> stdvec_mat;
typedef std::vector<arma::vec> stdvec_vec;

// ---------------------------------------------------------------------
// VON MISES SINE
// ---------------------------------------------------------------------

// [[Rcpp::export]]
double const_vmsin(double k1, double k2, double lambda)
{

  double start = BESSI(0,k1) * BESSI(0,k2), rat = pow(lambda, 2) / (4*k1*k2);
  int p = 1;
  double temp = Rf_choose(2*p, p) * pow(rat, p) * BESSI(p, k1) * BESSI(p, k2);
  start += temp;
  while((temp/start > 1e-6) && (p <= 10000)) {
    p++;
    temp = Rf_choose(2*p, p) * pow(rat, p) * BESSI(p, k1) * BESSI(p, k2);
    start += temp;
  }
  return (4 * M_PI * M_PI * start);
}



// // [[Rcpp::export]]
// double d_const_vmsin_lambda(double k1, double k2, double lambda)
// {
//   double result = 0;
//   if(lambda != 0){
//     double start = 0, rat = pow(lambda, 2) / (4*k1*k2);
//     int p = 1;
//     double temp = Rf_choose(2*p, p) * p * pow(rat, p) * BESSI(p, k1) * BESSI(p, k2);
//     start += temp;
//     while((temp/start > 1e-6) && (p <= 10000)) {
//       p++;
//       temp = Rf_choose(2*p, p) * p * pow(rat, p) * BESSI(p, k1) * BESSI(p, k2);
//       start += temp;
//     }
//     result = (8/lambda) * M_PI * M_PI * start;
//   }
//   return (result);
// }
//
// // [[Rcpp::export]]
// double d_const_vmsin_k1(double k1, double k2, double lambda)
// {
//   double start = BESSI(1,k1) * BESSI(0,k2) * 2;
//   if(lambda != 0){
//     double rat = pow(lambda, 2) / (4*k1*k2);
//     int p = 1;
//     double temp = Rf_choose(2*p, p) * pow(rat, p) * BESSI(p, k2) * ( BESSI((p+1), k1) + BESSI((p-1), k1) - (2*p/k1)*BESSI(p, k1));
//     start += temp;
//     while((temp/start > 1e-6) && (p <= 10000)) {
//       p++;
//       temp = Rf_choose(2*p, p) * pow(rat, p) * BESSI(p, k2) * ( BESSI((p+1), k1) + BESSI((p-1), k1) - (2*p/k1)*BESSI(p, k1));
//       start += temp;
//     }
//   }
//   return (2 * M_PI * M_PI * start);
// }
//
// // [[Rcpp::export]]
// double d_const_vmsin_k2(double k1, double k2, double lambda)
// {
//   return (d_const_vmsin_k1(k2, k1, lambda));
// }



// [[Rcpp::export]]
arma::vec d_const_vmsin(arma::vec par)
{
  double k1 = par[0], k2 = par[1], lambda = par[2], c = 0, c_k1 = 0, c_k2 = 0, c_lambda = 0, rat = lambda*lambda / (4*k1*k2);

  int p = 0;

  double choose = Rf_choose(2*p, p),
    rat_pow_p = pow(rat, p),
    I_p_k1 = BESSI(p, k1),
    I_p_k2 = BESSI(p, k2),
    I_pplus1_k1 = BESSI(p+1, k1),
    I_pplus1_k2 = BESSI(p+1, k2);
  // I_pplus2_k1 = BESSI(p+2, k1),
  // I_pplus2_k2 = BESSI(p+2, k2);

  double temp_c = choose *  rat_pow_p * I_p_k1 * I_p_k2,
    temp_c_k1 = choose * rat_pow_p * I_pplus1_k1 * I_p_k2,
    temp_c_k2 = choose * rat_pow_p * I_p_k1 * I_pplus1_k2,
    temp_c_lambda = p * choose *  rat_pow_p * I_p_k1 * I_p_k2;
  // temp_c_k1k2 = choose * rat_pow_p * I_pplus1_k1 * I_pplus1_k2,
  // temp_c_k1k1 = choose * rat_pow_p * (I_pplus1_k1/k1 + I_pplus2_k1) * I_p_k2,
  // temp_c_k2k2 = choose * rat_pow_p * (I_pplus1_k2/k2 + I_pplus2_k2) * I_p_k1;

  c += temp_c;
  c_k1 += temp_c_k1;
  c_k2 += temp_c_k2;
  if(lambda != 0) c_lambda += temp_c_lambda;
  // c_k1k2 += temp_c_k1k2;
  // c_k1k1 += temp_c_k1k1;
  // c_k2k2 += temp_c_k2k2;

  while(temp_c/c > 1e-6) {
    p++;

    choose = Rf_choose(2*p, p);
    rat_pow_p = pow(rat, p);
    I_p_k1 = I_pplus1_k1;
    I_p_k2 = I_pplus1_k2;
    I_pplus1_k1 = BESSI(p+1, k1);
    I_pplus1_k2 = BESSI(p+1, k2);
    // I_pplus2_k1 = BESSI(p+2, k1);
    // I_pplus2_k2 = BESSI(p+2, k2);

    temp_c = choose *  rat_pow_p * I_p_k1 * I_p_k2;
    temp_c_k1 = choose * rat_pow_p * I_pplus1_k1 * I_p_k2;
    temp_c_k2 = choose * rat_pow_p * I_p_k1 * I_pplus1_k2;
    if(lambda != 0) temp_c_lambda = p * choose *  rat_pow_p * I_p_k1 * I_p_k2;
    // temp_c_k1k2 = choose * rat_pow_p * I_pplus1_k1 * I_pplus1_k2;
    // temp_c_k1k1 = choose * rat_pow_p * (I_pplus1_k1/k1 + I_pplus2_k1) * I_p_k2;
    // temp_c_k2k2 = choose * rat_pow_p * (I_pplus1_k2/k2 + I_pplus2_k2) * I_p_k1;

    c += temp_c;
    c_k1 += temp_c_k1;
    c_k2 += temp_c_k2;
    if(lambda != 0) c_lambda += temp_c_lambda;
    // c_k1k2 += temp_c_k1k2;
    // c_k1k1 += temp_c_k1k1;
    // c_k2k2 += temp_c_k2k2;

  }


  double pi_sq_times_4 = 4 * M_PI * M_PI;

  arma::vec out;

  if(lambda != 0) {
    out = {
      c_k1 * pi_sq_times_4,
      c_k2 * pi_sq_times_4,
      c_lambda * pi_sq_times_4 * 2/lambda
    };

  } else {
    out = {
      c_k1 * pi_sq_times_4,
      c_k2 * pi_sq_times_4,
      0.0
    };
  }

  return out;

}




// [[Rcpp::export]]
arma::vec log_const_vmsin_all(arma::mat par_mat) {
  int K = par_mat.n_cols;
  arma::vec all_lconsts(K);
  for(int j = 0; j < K; j++)
    all_lconsts[j] = log(const_vmsin(par_mat(0,j), par_mat(1,j), par_mat(2,j)));
  return all_lconsts;
}


// [[Rcpp::export]]
double ldsinnum(double x, double y, arma::vec par)
{
  double k1 = par[0], k2 = par[1], lambda = par[2], mu = par[3], nu = par[4];
  return (k1*cos(x-mu) + k2*cos(y-nu) + lambda*sin(x-mu)*sin(y-nu));
}



// [[Rcpp::export]]
arma::vec dsin_onex_manypar(arma::vec x, arma::vec k1, arma::vec k2, arma::vec k3,
                            arma::vec mu1, arma::vec mu2)
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

  arma::vec l_const_all = log_const_vmsin_all(all_par);
  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldsinnum(x[0], x[1], all_par.col(i));
  }

  return arma::exp(ld_num - l_const_all);
}

// [[Rcpp::export]]
arma::vec dsin_manyx_onepar(arma::mat x, double k1, double k2, double k3,
                            double mu1, double mu2)
{
  int n = x.n_rows;

  double l_const = log(const_vmsin(k1, k2, k3));
  arma::vec par(5);
  par[0] = k1; par[1] = k2; par[2] = k3; par[3] = mu1; par[4] = mu2;

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldsinnum(x(i,0), x(i,1), par);
  }

  return arma::exp(ld_num - l_const);
}

// [[Rcpp::export]]
arma::vec dsin_manyx_manypar(arma::mat x, arma::vec k1, arma::vec k2, arma::vec k3,
                             arma::vec mu1, arma::vec mu2)
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

  arma::vec l_const_all = log_const_vmsin_all(all_par);

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldsinnum(x(i,0), x(i,1), all_par.col(i));
  }

  return arma::exp(ld_num - l_const_all);
}


// // generates a single observation
// arma::rowvec rsin_single_onepar(double k1, double k2, double k3,
//                                 double mu1, double mu2, double I_upper_bd) {
//   arma::vec par(5);
//   par << k1 << k2 << k3 << mu1 << mu2 << arma::endr;
//
//   double x, y, U1, sin_x_minus_mu1, a;
//
//   // First generate x from the marginal distribution given in Eqn. 2.2 in Singh et al 2002
//
//   int accpt_x = 0, rej = 0;
//   while (accpt_x < 1) {
//     x = runivm_single_onepar(k1, mu1);
//     sin_x_minus_mu1 = sin(x - mu1);
//     a = sqrt(k2 * k2 + k3 * k3 * sin_x_minus_mu1 * sin_x_minus_mu1);
//     U1 = R::unif_rand();
//     // I_upper_bd is numerically calculated in R
//     if (U1 <= BESSI0(a) / I_upper_bd)
//       accpt_x++;
//     else
//       rej++;
//     if(rej == 50) I_upper_bd = BESSI0(sqrt(k2 * k2 + k3 * k3));
//     // if there are lots of rejections, use the conservative upper bound.
//     if(rej == 5e4) warning("5e4 rejetions.. continuing..");
//   }
//
//   double beta = atan(k3 / k2 * sin_x_minus_mu1);
//
//   // Now generate y from vM(a, mu2 + beta) (Eqn. 2.4, Singh et al. 2002)
//   y = runivm_single_onepar(a, (mu2 + beta));
//
//   arma::rowvec final_sample(2);
//   final_sample[0] = x;
//   final_sample[1] = y;
//
//   return final_sample;
// }

// generates a single observation
arma::rowvec2 rsin_unimodal_single(double k1, double k2, double k3,
                                   double mu1, double mu2,
                                   double kappa_opt, double log_I0_kappa_opt,
                                   double logK,
                                   double log_const_vmsin) {
  double x, y, U1;


  double log_prop_den, log_target_den,
  log_2pi=log(2*M_PI), cos_x_mu1, sin_x_mu1, a_x;

  // do marginal x then conditional y|x generation


  // generating from marginal of x:
  // for marginal of x, proposal is vm(kappa_optim, mu1)
  int accpt_y = 0, ntry = 1; // will be 1 if accepted in rejection sampling from y marginal
  // First generate y from the marginal distribution given in Eqn. 5 in Mardia 2007
  while (accpt_y < 1) {
    x = runivm_single_onepar(kappa_opt, mu1);
    cos_x_mu1 = cos(x-mu1);
    sin_x_mu1 = sin(x-mu1);
    a_x = sqrt(k2*k2 + k3*k3*sin_x_mu1*sin_x_mu1);
    log_target_den = -log_const_vmsin + log_2pi + log(BESSI0(a_x)) + k2*cos_x_mu1;
    log_prop_den = -log_2pi - log_I0_kappa_opt + kappa_opt*cos_x_mu1;
    // Rcout << ntry << "\t";
    ntry++;
    U1 = R::unif_rand();
    if (log(U1) <= log_target_den-log_prop_den-logK)
      accpt_y++;
  }


  double mean_y = mu2+atan(k3/k2*sin_x_mu1);

  // generating x|y
  // generate x from vM(a_x, mean_y)
  y = runivm_single_onepar(a_x, mean_y);

  arma::rowvec2 final_sample = { x, y };

  return final_sample;
}



// [[Rcpp::export]]
arma::mat rsin_unimodal(int n, double k1, double k2, double k3,
                        double mu1, double mu2,
                        double kappa_opt, double log_I0_kappa_opt,
                        double logK,
                        double log_const_vmsin)
{
  arma::mat out(n, 2);
  for (int i = 0; i < n; i++)
    out.row(i) =
      rsin_unimodal_single(k1, k2, k3, mu1, mu2, kappa_opt,
                           log_I0_kappa_opt, logK, log_const_vmsin);

  return out;
}


arma::rowvec2 rsin_bimodal_single(double k1, double k2, double k3,
                                  double mu1, double mu2,
                                  double kappa_opt, double log_I0_kappa_opt,
                                  double logK,
                                  double log_const_vmsin,
                                  double mode_1, double mode_2,
                                  double vmpropn, double unifpropn)

{
  double x, y, U1, U0;



  double log_prop_den, log_target_den,
  log_2pi=log(2*M_PI), cos_x_mu1, sin_x_mu1, a_x;

  // do marginal x then conditional y|x generation


  // generating from marginal of x:


  // for marginal of x, proposal is vmpropn-vmpropn-unifpropn
  // mixture of vm(kappa_optim, mode_1), vm(kappa_optim, mode_2)
  // unif(0, 2*pi)
  int accpt_y = 0, ntry = 1; // will be 1 if accepted in rejection sampling from y marginal
  // First generate y from the marginal distribution given in Eqn. 5 in Mardia 2007
  while (accpt_y < 1) {

    U0 = R::unif_rand();

    if (U0 < vmpropn) {
      x = runivm_single_onepar(kappa_opt, mode_1);
    } else if (U0 < 2*vmpropn) {
      x = runivm_single_onepar(kappa_opt, mode_2);
    } else {
      x = R::runif(0, 2*M_PI);
    }
    cos_x_mu1 = cos(x-mu1);
    sin_x_mu1 = sin(x-mu1);
    a_x = sqrt(k2*k2 + k3*k3*sin_x_mu1*sin_x_mu1);
    log_target_den = -log_const_vmsin + log_2pi + log(BESSI0(a_x)) + k2*cos_x_mu1;
    log_prop_den =
      -log_2pi +
      log(vmpropn*exp(-log(BESSI0(kappa_opt))
                        + kappa_opt*cos(x-mode_1)
                        + log(1+exp(kappa_opt*(cos(x-mode_2) - cos(x-mode_1))))
      )
            + unifpropn
      );
    ntry++;
    if (ntry % 100 == 0) checkUserInterrupt();
    // Rcout << ntry << "\t";


    U1 = R::unif_rand();
    if (log(U1) <= log_target_den-log_prop_den-logK)
      accpt_y++;
  }



  double mean_y = mu2+atan(k3/k2*sin_x_mu1);

  // generating x|y


  // generate y from vM(a_x, mean_y)
  y = runivm_single_onepar(a_x, mean_y);

  arma::rowvec2 final_sample = { x, y };

  return final_sample;
}


//  [[Rcpp::export]]
arma::mat rsin_bimodal(int n, double k1, double k2, double k3,
                       double mu1, double mu2,
                       double kappa_opt, double log_I0_kappa_opt,
                       double logK,
                       double log_const_vmsin,
                       double mode_1, double mode_2,
                       double vmpropn, double unifpropn)

{
  arma::mat out(n, 2);
  for (int i = 0; i < n; i++)
    out.row(i) =
      rsin_bimodal_single(k1, k2, k3, mu1, mu2, kappa_opt,
                          log_I0_kappa_opt, logK, log_const_vmsin,
                          mode_1, mode_2, vmpropn, unifpropn);

  return out;
}






// generates a single observation
arma::rowvec rsin_single_onepar(double k1, double k2, double k3,
                                double mu1, double mu2, double I_upper_bd) {
  arma::vec par = { k1, k2, k3, mu1, mu2 };


  double x, y, U1;

  // brute force bivariate rejection sampling
  // with numerically computed upper bound I_upper_bd (passed from R)
  int accpt = 0;
  while (accpt < 1) {
    x = R::runif(0, 2*M_PI);
    y = R::runif(0, 2*M_PI);
    U1 = R::unif_rand();
    if (log(U1) <= (k1*cos(x-mu1)+k2*cos(y-mu2)+k3*sin(x-mu1)*sin(y-mu2))-I_upper_bd) {
      accpt++;
    }
  }

  arma::rowvec final_sample(2);
  final_sample[0] = x;
  final_sample[1] = y;

  return final_sample;
}






//  [[Rcpp::export]]
arma::mat rsin_onepar(int n, double k1, double k2, double k3,
                      double mu1, double mu2, double I_upper_bd) {
  if(n == 1){
    return(rsin_single_onepar(k1, k2, k3, mu1, mu2, I_upper_bd));
  } else {
    arma::mat all_sample(n, 2);
    for(int i = 0; i < n; i++)
      all_sample.row(i) = rsin_single_onepar(k1, k2, k3, mu1, mu2, I_upper_bd);
    return(all_sample);
  }
}


// // not needed
// //  [[Rcpp::export]]
// arma::mat rsin_manypar(arma::vec k1, arma::vec k2, arma::vec k3,
//                        arma::vec mu1, arma::vec mu2, arma::vec I_upper_bd) {
//   int n = k1.size();
//   arma::mat all_sample(n, 2);
//   for(int i = 0; i < n; i++)
//     all_sample.row(i) = rsin_single_onepar(k1[i], k2[i], k3[i], mu1[i], mu2[i], I_upper_bd[i]);
//   return(all_sample);
// }



// [[Rcpp::export]]
arma::mat mem_p_sin(arma::mat data, arma::mat par, arma::vec pi,
                    arma::vec log_c_von, int ncores = 1)
{
  int n = data.n_rows, K = par.n_cols,  j;
  double row_total;
  arma::mat den(n, K);
  for(int i = 0; i < n; i++){
    row_total = 0;
    for(j = 0; j < K; j++){
      den(i, j) = pi[j]*exp(ldsinnum(data(i, 0), data(i, 1), par.col(j)) - log_c_von[j]);;
      row_total += den(i, j);
    }
    row_total = maxi(row_total, 1e-50);
    for(j = 0; j < K; j++)
      den(i, j) /= row_total;
  }
  return(den);
}


// // not needed
// // [[Rcpp::export]]
// double llik_vmsin_full(arma::mat data, arma::mat par, arma::vec pi, arma::vec log_c, int ncores = 1)
// {
//
//   int n = data.n_rows, K = pi.size() ,j;
//   doubletemp, log_sum = 0.0;
//   arma::vec log_pi = log(pi);
//
//   if(K > 1) {
//     for(int i = 0; i < n; i++) {
//       temp = 0;
//       for(j = 0; j < K; j++)
//         temp += exp(ldsinnum(data(i,0), data(i,1), par.col(j)) - log_c[j] + log_pi[j] );
//       log_sum += log(temp);
//     }
//   }  else {
//     for(int i = 0; i < n; i++)
//       log_sum += ldsinnum(data(i,0), data(i,1), par);
//     log_sum -= n*log_c[0];
//   }
//   return(log_sum);
//
// }


// [[Rcpp::export]]
arma::vec llik_vmsin_contri_C(arma::mat data, arma::mat par, arma::vec pi, arma::vec log_c)
{

  int n = data.n_rows, K = pi.size() ,j;
  double temp;
  arma::vec llik_contri = arma::vec(n);
  arma::vec log_pi = log(pi);

  if(K > 1) {
    for(int i = 0; i < n; i++) {
      temp = 0;
      for(j = 0; j < K; j++)
        temp += exp(ldsinnum(data(i,0), data(i,1), par.col(j)) - log_c[j] + log_pi[j] );
      llik_contri[i] = log(temp);
    }
  }
  else {
    for(int i = 0; i < n; i++) {
      llik_contri[i] = ldsinnum(data(i,0), data(i,1), par) - log_c[0];
    }
  }
  return(llik_contri);

}




// [[Rcpp::export]]
double llik_vmsin_one_comp(arma::mat data, arma::vec par_vec,
                           double log_c)
{

  int n = data.n_rows;
  double log_sum = 0.0;

  for(int i = 0; i < n; i++)
    log_sum += ldsinnum(data(i,0), data(i,1), par_vec);
  log_sum -= n*log_c;
  return (log_sum);

}




// [[Rcpp::export]]
arma::vec grad_llik_vmsin_C(arma::mat data, arma::vec par) {
  int n = data.n_rows;
  double k1 = par[0], k2 = par[1], k3 = par[2], mu1 = par[3], mu2 = par[4];

  arma::vec grad_llik_wo_const = arma::zeros(6);
  // first 5 entries = grad, last entry = llik

  for(int i = 0; i < n; i++) {
    double x_mu1 = (data(i, 0) - mu1), y_mu2 = (data(i, 1) - mu2);
    double  cos_x_mu1 = cos(x_mu1), sin_x_mu1 = sin(x_mu1),
      cos_y_mu2 = cos(y_mu2), sin_y_mu2 = sin(y_mu2);

    grad_llik_wo_const[0] += cos_x_mu1; //(cos(x - mu1) - del_const_vmsin[0]/c_vmsin);
    grad_llik_wo_const[1] += cos_y_mu2; // - del_const_vmsin[1]/c_vmsin);
    grad_llik_wo_const[2] += sin_x_mu1 * sin_y_mu2; // - del_const_vmsin[2]/c_vmsin);
    grad_llik_wo_const[3] += k1*sin_x_mu1 - k3*cos_x_mu1*sin_y_mu2;
    grad_llik_wo_const[4] += k2*sin_y_mu2 - k3*sin_x_mu1*cos_y_mu2;
    grad_llik_wo_const[5] += k1*cos_x_mu1 + k2*cos_y_mu2 + k3*sin_x_mu1*sin_y_mu2;
  }

  double c_vmsin = const_vmsin(k1, k2, k3);
  arma::vec del_const_vmsin = d_const_vmsin(par);

  // adjust const terms in the first three pars
  for(int j = 0; j < 3; j++) {
    grad_llik_wo_const[j] -= n * del_const_vmsin[j]/c_vmsin;
  }

  // subtract constant term from llik
  grad_llik_wo_const[5] -= n * log(c_vmsin);

  return(grad_llik_wo_const);
}


// // not needed
// arma::vec grad_llik_vmsin_one_comp_i(double x, double y, arma::vec par,
//                                      double c_vmsin, arma::vec del_const_vmsin) {
//   double k1 = par[0], k2 = par[1], lambda = par[2], mu1 = par[3], mu2 = par[4];
//   arma::vec all_entries(6);
//   all_entries[0] = (cos(x - mu1) - del_const_vmsin[0]/c_vmsin);
//   all_entries[1] = (cos(y - mu2) - del_const_vmsin[1]/c_vmsin);
//   all_entries[2] = (sin(x - mu1)*sin(y - mu2) - del_const_vmsin[2]/c_vmsin);
//   all_entries[3] = (k1*sin(x - mu1) - lambda*cos(x - mu1)*sin(y - mu2));
//   all_entries[4] = (k2*sin(y - mu2) - lambda*sin(x - mu1)*cos(y - mu2));
//   all_entries[6] = ldsinnum(x, y, par) - log(c_vmsin);
//   return all_entries;
// }
//


//
// // not needed
// arma::vec grad_log_vmsin_one_comp_i(double x, double y, arma::vec par,
//                                     double c_vmsin, arma::vec del_const_vmsin) {
//   double k1 = par[0], k2 = par[1], lambda = par[2], mu1 = par[3], mu2 = par[4];
//   arma::vec all_entries(5);
//   all_entries[0] = (cos(x - mu1) - del_const_vmsin[0]/c_vmsin);
//   all_entries[1] = (cos(y - mu2) - del_const_vmsin[1]/c_vmsin);
//   all_entries[2] = (sin(x - mu1)*sin(y - mu2) - del_const_vmsin[2]/c_vmsin);
//   all_entries[3] = (k1*sin(x - mu1) - lambda*cos(x - mu1)*sin(y - mu2));
//   all_entries[4] = (k2*sin(y - mu2) - lambda*sin(x - mu1)*cos(y - mu2));
//
//   return all_entries;
// }
//
//
// // not needed
// // [[Rcpp::export]]
// arma::vec grad_vmsin_one_comp(arma::mat data, arma::vec par_vec,
//                               int ncores = 1)
// {
//   int n = data.n_rows;
//   double c_vmsin = const_vmsin(par_vec[0], par_vec[1], par_vec[2]);
//
//   arma::vec del_const_vmsin = d_const_vmsin(par_vec);
//
//   stdvec_vec grad_sum_part(ncores);
//   for(int g = 0; g < ncores; g++)
//     grad_sum_part[g] = arma::zeros(5);
//
//
//   for(int i = 0; i < n; i++) {
//     arma::vec grad_temp(5);
//     int g = get_thread_num_final();
//     // int g = 0;
//
//     // temp = exp(ldsinnum(data(i, 0), data(i, 1), par_mat_cols[j]) - l_c_vmsin[j] + log_pi[j]);
//     // denom += temp;
//     grad_temp = grad_log_vmsin_one_comp_i(data(i, 0), data(i, 1),
//                                           par_vec, c_vmsin, del_const_vmsin);
//
//     grad_sum_part[g] += grad_temp;
//   }
//
//   arma::vec grad_sum = arma::zeros(5);
//   for(int g = 0; g < ncores; g++)
//     grad_sum += grad_sum_part[g];
//
//   return (grad_sum);
// }
//
//
// // [[Rcpp::export]]
// arma::vec grad_llik_vmsin_C(arma::mat data, arma::vec par_vec)
// {
//
//     int n = data.n_rows;
//     double c_vmsin = const_vmsin(par_vec[0], par_vec[1], par_vec[2]);
//
//     arma::vec del_const_vmsin = d_const_vmsin(par_vec);
//
//     arma::vec grad_sum = arma::zeros(6);
//
//
//     for(int i = 0; i < n; i++) {
//
//       grad_sum += grad_llik_vmsin_one_comp_i(data(i, 0), data(i, 1),
//                                              par_vec, c_vmsin, del_const_vmsin);
//
//     }
//
//
//     return (grad_sum);
//   }


// // not needed
//
// // [[Rcpp::export]]
// arma::mat grad_vmsin_all_comp(arma::mat data, arma::mat par_mat, arma::vec pi, int ncores = 1)
// {
//   int n = data.n_rows, K = pi.size(), j;
//   arma::vec l_c_vmsin = log_const_vmsin_all(par_mat);
//   arma::vec c_vmsin = arma::exp(l_c_vmsin), log_pi = log(pi);
//
//   stdvec_vec del_const_vmsin_all(K);
//   stdvec_vec par_mat_cols(K);
//
//   for(j = 0; j < K; j++) {
//     del_const_vmsin_all[j] = d_const_vmsin(par_mat.col(j));
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
//       temp = exp(ldsinnum(data(i, 0), data(i, 1), par_mat_cols[j]) - l_c_vmsin[j] + log_pi[j]);
//       denom += temp;
//       grad_temp.col(j) = temp*grad_log_vmsin_one_comp_i(data(i, 0), data(i, 1), par_mat_cols[j], c_vmsin[j], del_const_vmsin_all[j]);
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
//
//
//
// // not needed
// // [[Rcpp::export]]
// double vmsinmix(double x, double y, arma::mat par, arma::vec pi, arma::vec log_c_von)
// {
//   double res = 0;
//   int K = par.n_cols;
//
//   for(int j = 0; j < K; j++)
//     res += exp(ldsinnum(x, y, par.col(j)) - log_c_von[j]) * pi[j];
//   return(res);
// } //likelihood contribution (density) of a single point in a mixture model
//
//
// // not needed
// // [[Rcpp::export]]
// arma::vec vmsinmix_manyx(arma::mat x, arma::mat par, arma::vec pi, arma::vec log_c_von)
// {
//   int n = x.n_rows;
//   arma::vec result(n);
//   for(int i = 0; i < n; i++)
//     result[i] = vmsinmix(x(i,0), x(i,1), par, pi, log_c_von);
//   return(result);
// }




// [[Rcpp::export]]
List vmsin_var_corr_anltc(double k1, double k2, double lambda)
{

  double c = 0, c_k1 = 0, c_k2 = 0, c_lambda = 0, c_k1k2 = 0,
    c_k1k1 = 0, c_k2k2 = 0, rat = lambda*lambda / (4*k1*k2);

  int p = 0;

  double choose = Rf_choose(2*p, p),
    rat_pow_p = pow(rat, p),
    I_p_k1 = BESSI(p, k1),
    I_p_k2 = BESSI(p, k2),
    I_pplus1_k1 = BESSI(p+1, k1),
    I_pplus1_k2 = BESSI(p+1, k2),
    I_pplus2_k1 = BESSI(p+2, k1),
    I_pplus2_k2 = BESSI(p+2, k2);

  double temp_c = choose *  rat_pow_p * I_p_k1 * I_p_k2,
    temp_c_k1 = choose * rat_pow_p * I_pplus1_k1 * I_p_k2,
    temp_c_k2 = choose * rat_pow_p * I_p_k1 * I_pplus1_k2,
    temp_c_lambda = p * choose *  rat_pow_p * I_p_k1 * I_p_k2,
    temp_c_k1k2 = choose * rat_pow_p * I_pplus1_k1 * I_pplus1_k2,
    temp_c_k1k1 =
      (k1 == 0) ?
      0 : choose * rat_pow_p * (I_pplus1_k1/k1 + I_pplus2_k1) * I_p_k2,
      temp_c_k2k2 =
        (k2 == 0) ?
        0 : choose * rat_pow_p * (I_pplus1_k2/k2 + I_pplus2_k2) * I_p_k1;

  c += temp_c;
  c_k1 += temp_c_k1;
  c_k2 += temp_c_k2;
  if(lambda != 0) c_lambda += temp_c_lambda;
  c_k1k2 += temp_c_k1k2;
  c_k1k1 += temp_c_k1k1;
  c_k2k2 += temp_c_k2k2;

  double den_min = std::min(c, std::min(k1, k2));

  while(temp_c/den_min > 1e-6) {
    p++;

    choose = Rf_choose(2*p, p);
    rat_pow_p = pow(rat, p);
    I_p_k1 = I_pplus1_k1;
    I_p_k2 = I_pplus1_k2;
    I_pplus1_k1 = I_pplus2_k1;
    I_pplus1_k2 = I_pplus2_k2;
    I_pplus2_k1 = BESSI(p+2, k1);
    I_pplus2_k2 = BESSI(p+2, k2);

    temp_c = choose *  rat_pow_p * I_p_k1 * I_p_k2;
    temp_c_k1 = choose * rat_pow_p * I_pplus1_k1 * I_p_k2;
    temp_c_k2 = choose * rat_pow_p * I_p_k1 * I_pplus1_k2;
    if (lambda != 0) temp_c_lambda = p * choose *  rat_pow_p * I_p_k1 * I_p_k2;
    temp_c_k1k2 = choose * rat_pow_p * I_pplus1_k1 * I_pplus1_k2;
    if (k1 != 0) temp_c_k1k1 = choose * rat_pow_p * (I_pplus1_k1/k1 + I_pplus2_k1) * I_p_k2;
    if (k1 != 0) temp_c_k2k2 = choose * rat_pow_p * (I_pplus1_k2/k2 + I_pplus2_k2) * I_p_k1;

    c += temp_c;
    c_k1 += temp_c_k1;
    c_k2 += temp_c_k2;
    if(lambda != 0) c_lambda += temp_c_lambda;
    c_k1k2 += temp_c_k1k2;
    c_k1k1 += temp_c_k1k1;
    c_k2k2 += temp_c_k2k2;

  }


  double pi_sq_times_4 = 4 * M_PI * M_PI;



  c *= pi_sq_times_4;
  c_k1 *= pi_sq_times_4;
  c_k2 *= pi_sq_times_4;
  if(lambda != 0) c_lambda *= pi_sq_times_4 * 2/lambda;
  c_k1k2 *= pi_sq_times_4;
  c_k1k1 *= pi_sq_times_4;
  c_k2k2 *= pi_sq_times_4;

  // double
  //   rho_js = c_lambda / sqrt((c-c_k1k1) * (c-c_k2k2)),
  //
  //   rho_fl = rho_js * c_k1k2 / sqrt(c_k1k1 * c_k2k2),
  //
  //     var1 = 1 - (c_k1/c),
  //     var2 = 1 - (c_k2/c);

  // rho_js = (c_k3) /
  //   sqrt((c-c_k1k1) * (c-c_k2k2)),
  double
    rho_js = (fabs(c_lambda) < 1e-10) ? 0 :
    sgn(c_lambda) *
      fmin(1, exp(log(plus(fabs(c_lambda)))
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
List vmsin_var_corr_mc(double k1, double k2, double k3,
                       arma::mat uni_rand, int ncores = 1)
{
  double x, y, expon, cos_x, cos_y, sin_x, sin_y,
  sum_exp, sum_exp1, sum_exp2, sum_exp3, sum_exp11,
  sum_exp12, sum_exp22,
  temp,
  // exp_temp,
  exp_expon_temp;

  int nsim = uni_rand.n_rows;

  x = uni_rand(0, 0) * 2 * M_PI;
  y = uni_rand(0, 1) * 2 * M_PI;
  cos_x = cos(x);
  cos_y = cos(y);
  sin_x = sin(x);
  sin_y = sin(y);
  sum_exp = 1.0;
  sum_exp1 = cos_x;
  sum_exp2 = cos_y;
  sum_exp3 = sin_x*sin_y;
  sum_exp11 = cos_x*cos_x;
  sum_exp22 = cos_y*cos_y;
  sum_exp12 = cos_x*cos_y;
  temp = k1 * cos_x + k2 * cos_y + k3 * sin_x * sin_y;
  // exp_temp = exp(temp);

  for(int i = 1; i < nsim; i++) {
    x = uni_rand(i, 0) * 2 * M_PI;
    y = uni_rand(i, 1) * 2 * M_PI;
    cos_x = cos(x);
    cos_y = cos(y);
    sin_x = sin(x);
    sin_y = sin(y);
    expon = k1 * cos_x + k2 * cos_y + k3 * sin_x * sin_y;
    exp_expon_temp = exp(expon - temp);

    sum_exp += exp_expon_temp;
    sum_exp1 += exp_expon_temp *  cos_x;
    sum_exp2 += exp_expon_temp *  cos_y;
    sum_exp3 += exp_expon_temp *   sin_x * sin_y;
    sum_exp11 += exp_expon_temp *  cos_x*cos_x;
    sum_exp22 += exp_expon_temp *  cos_y*cos_y;
    sum_exp12 += exp_expon_temp *  cos_x*cos_y;

  }



  // double pisq_4_exp_temp_over_nsim = 4 * M_PI * M_PI * exp_temp / nsim;
  double pisq_4_exp_temp_over_nsim = 4 * M_PI * M_PI / nsim;

  // all the following need to be multiplied by exp(temp)
  double c = sum_exp * pisq_4_exp_temp_over_nsim,
    c_k1 = sum_exp1 * pisq_4_exp_temp_over_nsim,
    c_k2 = sum_exp2 * pisq_4_exp_temp_over_nsim,
    c_k3 = sum_exp3 * pisq_4_exp_temp_over_nsim,
    c_k1k2 = sum_exp12 * pisq_4_exp_temp_over_nsim,
    c_k1k1 = sum_exp11 * pisq_4_exp_temp_over_nsim,
    c_k2k2 = sum_exp22 * pisq_4_exp_temp_over_nsim;


  double
    rho_js = (fabs(c_k3) < 1e-10) ? 0 :
    sgn(c_k3) *
      fmin(1, exp(log(plus(fabs(c_k3)))
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
List vmsin_var_cor_singlepar_cpp(double k1, double k2, double k3,
                                 arma::mat uni_rand, int ncores = 1)
{
  // if unimodal, then analytic formula for derivatives,
  // else numeric
  if (k3*k3 < k1*k2 && k1 <= 50 && k2 <= 50 && fabs(k3) <= 50) {
    return vmsin_var_corr_anltc(k1, k2, k3);
  } else {
    return vmsin_var_corr_mc(k1, k2, k3, uni_rand, ncores);
  }

}





// [[Rcpp::export]]
arma::vec ldsin_onex_manypar(arma::vec x, arma::vec k1, arma::vec k2, arma::vec k3,
                             arma::vec mu1, arma::vec mu2)
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

  arma::vec l_const_all = log_const_vmsin_all(all_par);
  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldsinnum(x[0], x[1], all_par.col(i));
  }

  return (ld_num - l_const_all);
}

// [[Rcpp::export]]
arma::vec ldsin_manyx_onepar(arma::mat x, double k1, double k2, double k3,
                             double mu1, double mu2)
{
  int n = x.n_rows;

  double l_const = log(const_vmsin(k1, k2, k3));
  arma::vec par(5);
  par[0] = k1; par[1] = k2; par[2] = k3; par[3] = mu1; par[4] = mu2;

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldsinnum(x(i,0), x(i,1), par);
  }

  return (ld_num - l_const);
}


// [[Rcpp::export]]
arma::vec ldsin_manyx_manypar(arma::mat x, arma::vec k1, arma::vec k2, arma::vec k3,
                              arma::vec mu1, arma::vec mu2)
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

  arma::vec l_const_all = log_const_vmsin_all(all_par);

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldsinnum(x(i,0), x(i,1), all_par.col(i));
  }

  return (ld_num - l_const_all);
}

