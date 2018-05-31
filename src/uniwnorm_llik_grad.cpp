#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

typedef std::vector<arma::mat> stdvec_mat;

// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
//
// #include "thread_num.h"
#include "bessel.h"




// [[Rcpp::export]]
double l_const_uniwnorm(double k)
{
  return 0.5*log(2 * M_PI / k);
}


// [[Rcpp::export]]
double const_uniwnorm(double k)
{
  return exp(l_const_uniwnorm(k));
}


// [[Rcpp::export]]
arma::vec log_const_uniwnorm_all(arma::mat par_mat)
{
  int K = par_mat.n_cols;
  arma::vec result(K);
  for(int j = 0; j < K; j++)
    result[j] = l_const_uniwnorm(par_mat(0, j));
  return(result);
}



// [[Rcpp::export]]
double lduniwnormnum(double x, arma::vec par, arma::vec omega_2pi_1d)
{
  int M = omega_2pi_1d.size();
  double expon_sum = 0.0;
  for(int j = 0; j < M; j++) {
    expon_sum += exp(-0.5 * par[0] * (omega_2pi_1d[j] - x + par[1]) * (omega_2pi_1d[j] - x + par[1]));
  }
  return (log(expon_sum));
}

// [[Rcpp::export]]
arma::vec duniwnorm_manyx_onepar(arma::vec x, double k, double mu, arma::vec omega_2pi_1d)
{
  int n = x.n_rows;

  double l_const = log(const_uniwnorm(k));
  arma::vec par(2);
  par[0] = k; par[1] = mu;

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = lduniwnormnum(x[i], par, omega_2pi_1d);
  }

  return arma::exp(ld_num - l_const);
}

// [[Rcpp::export]]
arma::vec duniwnorm_manyx_manypar(arma::vec x, arma::vec k, arma::vec mu, arma::vec omega_2pi_1d)
{
  int n = k.size();

  arma::mat all_par(2, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k[i];
    all_par(1,i) = mu[i];
  }

  arma::vec l_const_all = log_const_uniwnorm_all(all_par);

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = lduniwnormnum(x[i], all_par.col(i), omega_2pi_1d);
  }

  return arma::exp(ld_num - l_const_all);
}


// [[Rcpp::export]]
arma::vec duniwnorm_onex_manypar(double x, arma::vec k, arma::vec mu, arma::vec omega_2pi_1d)
{
  int n = k.size();

  arma::mat all_par(2, n);
  for(int i = 0; i < n; i++) {
    all_par(0,i) = k[i];
    all_par(1,i) = mu[i];
  }

  arma::vec l_const_all = log_const_uniwnorm_all(all_par);

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = lduniwnormnum(x, all_par.col(i), omega_2pi_1d);
  }

  return arma::exp(ld_num - l_const_all);
}


// // not needed
// // [[Rcpp::export]]
// double uniwnormmix(double x,  arma::mat par, arma::vec pi, arma::vec log_c_von, arma::vec omega_2pi_1d)
// {
//   double res = 0;
//   int K = par.n_cols;
//
//   for(int j = 0; j < K; j++)
//     res += exp(lduniwnormnum(x, par.col(j), omega_2pi_1d) - log_c_von[j]) * pi[j];
//   return(res);
// } //likelihood contribution (density) of a single point in a mixture model
//
//
// // not used
// // [[Rcpp::export]]
// arma::vec uniwnormmix_manyx(arma::vec x, arma::mat par, arma::vec pi, arma::vec log_c, arma::vec omega_2pi_1d)
// {
//   int n = x.n_rows;
//   arma::vec result(n);
//   for(int i = 0; i < n; i++)
//     result[i] = uniwnormmix(x[i], par, pi, log_c, omega_2pi_1d);
//   return(result);
// }


// [[Rcpp::export]]
arma::mat mem_p_uniwnorm(arma::vec data, arma::mat par, arma::vec pi,
                         arma::vec log_c_von, arma::vec omega_2pi_1d)
{
  int n = data.n_rows, K = par.n_cols, j;
  double row_total;
  arma::mat den(n, K);
  for(int i = 0; i < n; i++){
    row_total = 0;
    for(j = 0; j < K; j++){
      den(i, j) = pi[j]*exp(lduniwnormnum(data[i], par.col(j), omega_2pi_1d) - log_c_von[j]);;
      row_total += den(i, j);
    }
    row_total = maxi(row_total, 1e-50);
    for(j = 0; j < K; j++)
      den(i, j) /= row_total;
  }
  return(den);

}

// // unused
// // [[Rcpp::export]]
// double llik_uniwnorm_full(arma::vec data, arma::mat par, arma::vec pi,
//                           arma::vec log_c, arma::vec omega_2pi_1d,
//                           int ncores = 1)
// {
//   int n = data.n_rows, K = pi.size(), j;
//   double temp, log_sum = 0.0;
//   arma::vec log_pi = log(pi);
//
//   if(K > 1) {
//     for(int i = 0; i < n; i++) {
//       temp = 0;
//       for(j = 0; j < K; j++)
//         temp += exp(lduniwnormnum(data[i], par.col(j), omega_2pi_1d) - log_c[j] + log_pi[j] );
//       log_sum += log(maxi(temp, 1e-100));
//     }
//   } else {
//     for(int i = 0; i < n; i++)
//       log_sum += lduniwnormnum(data[i], par, omega_2pi_1d);
//     log_sum -= n*log_c[0];
//   }
//   return(log_sum);
// }



// [[Rcpp::export]]
arma::vec llik_uniwnorm_contri_C(arma::vec data, arma::mat par, arma::vec pi,
                                 arma::vec log_c, arma::vec omega_2pi_1d)
{
  int n = data.n_rows, K = pi.size(), j;
  double temp;
  arma::vec llik_contri(n);
  arma::vec log_pi = log(pi);

  if(K > 1) {
    for(int i = 0; i < n; i++) {
      temp = 0;
      for(j = 0; j < K; j++)
        temp += exp(lduniwnormnum(data[i], par.col(j), omega_2pi_1d) - log_c[j] + log_pi[j] );
      llik_contri[i] = log(maxi(temp, 1e-100));
    }
  } else {
    for(int i = 0; i < n; i++)
      llik_contri[i] = lduniwnormnum(data[i], par.col(0), omega_2pi_1d) - log_c[0];
  }
  return(llik_contri);
}




// [[Rcpp::export]]
double llik_uniwnorm_one_comp(arma::vec data, arma::vec par_vec,
                              double log_c, arma::vec omega_2pi_1d)
{

  int n = data.n_rows;
  double log_sum = 0.0;

  for(int i = 0; i < n; i++) {
    log_sum += lduniwnormnum(data[i], par_vec, omega_2pi_1d);
  }

  log_sum -= n*log_c;

  return (log_sum);
}


arma::vec grad_log_den_uniwnorm_one_comp_i(double x,
                                                 arma::vec par,
                                                 arma::vec omega_2pi_1d)
{
  double k = par[0], mu = par[1];
  int M = omega_2pi_1d.n_rows;

  double expon_j;
  arma::vec all_entries = arma::zeros(3);
  // first two entries = grad, last entry = llik

  for(int j = 0; j < M; j++) {
    expon_j = exp(-0.5 * k * (x - mu - omega_2pi_1d[j]) * (x - mu - omega_2pi_1d[j]));
    all_entries[0] += expon_j * 0.5 / k *
      (1 - k * (x - mu - omega_2pi_1d[j])*(x - mu - omega_2pi_1d[j]));
    all_entries[1] += k * expon_j * (x - mu - omega_2pi_1d[j]);
    all_entries[2] += expon_j;
  }


  // divide 1st 2 entries by the third to get deriv(log density)
  all_entries[0] /= all_entries[2];
  all_entries[1] /= all_entries[2];

  // Finally, divide the last entry by normalizing constant, and take log
  all_entries[2] = log(all_entries[2]) - l_const_uniwnorm(k);

  return all_entries;
}


// [[Rcpp::export]]
arma::vec grad_llik_uniwnorm_C(arma::vec data, arma::vec par,
                               arma::vec omega_2pi_1d) {
  int n = data.n_elem;

  arma::vec grad_llik = arma::zeros(3);
  // first 2 entries = grad, last entry = llik

  for(int i = 0; i < n; i++) {
    grad_llik +=
      grad_log_den_uniwnorm_one_comp_i(data[i], par, omega_2pi_1d);
  }

  return(grad_llik);
}



// // unused
// arma::vec grad_den_uniwnorm_one_comp_i_unadj(double x, arma::vec par, arma::vec omega_2pi_1d)
// {
//   double k = par[0], mu = par[1];
//   int M = omega_2pi_1d.n_rows;
//   double expon_j;
//   arma::mat all_entries(3, M);  ;
//   for(int j = 0; j < M; j++) {
//     expon_j = exp(-0.5 * k * (x - mu - omega_2pi_1d[j]) * (x - mu - omega_2pi_1d[j]));
//     all_entries(0, j) = expon_j * (1 - k*(x - mu - omega_2pi_1d[j])*(x - mu - omega_2pi_1d[j]));// / (2.0*sqrt_2pi*sqrt_k);
//     all_entries(1, j) = expon_j * (x - mu - omega_2pi_1d[j]); // * (1.0*k*sqrt_k/sqrt_2pi);
//     all_entries(2, j) = expon_j;
//   }
//
//   arma::vec sum_entries = arma::zeros(3);
//   for(int j = 0; j < M; j++)
//     sum_entries += all_entries.col(j);
//   return sum_entries;
// }

// // unused
// // [[Rcpp::export]]
// arma::mat grad_uniwnorm_all_comp(arma::vec data, arma::mat par_mat, arma::vec pi,
//                                  arma::vec omega_2pi_1d, int ncores)
// {
//   int n = data.size(), K = pi.size(), j;
//   double denom;
//   arma::mat  grad_sum = arma::zeros(2, K);
//   arma::vec uniwnorm_const(K); // log_uniwnorm_const(K) = log_const_uniwnorm_all(par_mat), log_pi = log(pi);
//   for(j = 0; j < K; j++)
//     uniwnorm_const[j] = sqrt(par_mat(0, j)); // / (2 * M_PI));// 0.5*log(2.0 * M_PI / k)
//
//   stdvec_mat grad_sum_part(ncores);
//   for(int g = 0; g < ncores; g++)
//     grad_sum_part[g] = grad_sum;
//
//
//   for(int i = 0; i < n; i++) {
//     arma::mat grad_den_temp(3, K);
//     int g = get_thread_num_final();
//     denom = 0;
//     for(j = 0; j < K; j++){
//       // denom += exp(lduniwnormnum(data[i], par_mat.col(j), omega_2pi_1d) - log_uniwnorm_const[j] + log_pi[j] );
//       grad_den_temp.col(j) = grad_den_uniwnorm_one_comp_i_unadj(data[i], par_mat.col(j), omega_2pi_1d);
//       denom += grad_den_temp(2, j) * uniwnorm_const[j] * pi[j];
//     }
//     grad_sum_part[g] += grad_den_temp.head_rows(2)/denom;
//   }  // the pi_h term in numerator hasn't been adjusted yet!
//
//   for(int g = 0; g < ncores; g++)
//     grad_sum += grad_sum_part[g];
//
//
// //  double sqrt_2pi = sqrt((double) 2*M_PI);
//   arma::vec kappa = arma::trans(par_mat.row(0)), sqrt_kappa = arma::sqrt(kappa);
//   arma::mat adj_factor(2, K);
//
//   for(j = 0; j < K; j++){
//     adj_factor(0,j) = pi[j] / (2 * sqrt_kappa[j]); // /sqrt_2pi
//     adj_factor(1,j) = pi[j] * kappa[j] * sqrt_kappa[j]; // /sqrt_2pi);
//   }
//
//   return (grad_sum % adj_factor);
// }
