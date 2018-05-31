#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// #ifdef _OPENMP
// #include <omp.h>
// #endif
//
//
// #include "thread_num.h"

using namespace Rcpp;

// typedef std::vector<arma::mat> stdvec_mat;



// [[Rcpp::export]]
double ldwnorm2_num( arma::vec x, arma::vec par,  arma::mat omega_2pi) {
  int M = omega_2pi.n_rows;
  arma::vec allexpons(M);
  double expon_sum = 0.0, term1, term2;
  for(int j = 0; j < M; j++) {
    term1 = (omega_2pi(j,0) - x[0] + par[3]);
    term2 = (omega_2pi(j,1) - x[1] + par[4]);
    allexpons[j] = -0.5 *(par[0] * term1 * term1 + par[1] * term2 * term2 +
      2 * par[2] * term1 * term2);
    expon_sum += exp(allexpons[j]);
  }
  return (log(expon_sum));
}




// [[Rcpp::export]]
double l_const_wnorm2(arma::vec par) {
  double det_sig_inv = par[0] * par[1] - par[2]*par[2];
  return (log(2*M_PI) - 0.5*log(det_sig_inv));
}

// [[Rcpp::export]]
double const_wnorm2(arma::vec par) {

  return exp(l_const_wnorm2(par));
}



// [[Rcpp::export]]
arma::vec log_const_wnorm2_all(arma::mat par_mat) {
  int K = par_mat.n_cols;
  arma::vec result(K);
  for(int j = 0; j < K; j++)
    result[j] = l_const_wnorm2(par_mat.col(j));
  return result;
}

// [[Rcpp::export]]
arma::mat mem_p_wnorm2(arma::mat data, arma::mat par_mat, arma::vec pi,
                       arma::vec log_c_wnorm, arma::mat omega_2pi, int ncores = 1)
{
  int n = data.n_rows, K = par_mat.n_cols, j;
  double row_total;
  arma::mat den(n, K);
  arma::vec log_pi = log(pi);
  for(int i = 0; i < n; i++) {
    row_total = 0;
    for(j = 0; j < K; j++){
      den(i, j) = exp(ldwnorm2_num(arma::trans(data.row(i)), par_mat.col(j), omega_2pi) - log_c_wnorm[j] + log_pi[j]);
      row_total += den(i, j);
    }
    if(row_total <= 1e-50) row_total = 1e-50;
    for(j = 0; j < K; j++)
      den(i, j) /= row_total;
  }
  return(den);
}


// // not used
// // [[Rcpp::export]]
// double llik_wnorm2_full(arma::mat data, arma::mat par, arma::vec pi,
//                              arma::vec log_c, arma::mat omega_2pi, int ncores = 1)
// {
//   int n = data.n_rows, K = pi.size(), j;
//   double temp, log_sum = 0.0;
//   arma::vec log_pi = log(pi);
//
//   if(K > 1) {
//     for(int i = 0; i < n; i++) {
//       temp = 0;
//       for(j = 0; j < K; j++)
//         temp += exp(ldwnorm2_num(arma::trans(data.row(i)), par.col(j), omega_2pi) - log_c[j] + log_pi[j] );
//       log_sum += log(temp);
//     }
//   } else {
//     for(int i = 0; i < n; i++) {
//       log_sum += ldwnorm2_num(arma::trans(data.row(i)), par, omega_2pi);
//     }
//     log_sum -= n*log_c[0];
//   }
//   return(log_sum);
// }
//


// [[Rcpp::export]]
arma::vec llik_wnorm2_contri_C(arma::mat data, arma::mat par, arma::vec pi,
                               arma::vec log_c, arma::mat omega_2pi)
{
  int n = data.n_rows, K = pi.size(), j;
  double temp;
  arma::vec log_pi = log(pi);
  arma::vec llik_contri = arma::vec(n);


  if(K > 1) {
    for(int i = 0; i < n; i++) {
      temp = 0;
      for(j = 0; j < K; j++)
        temp += exp(ldwnorm2_num(arma::trans(data.row(i)), par.col(j), omega_2pi) - log_c[j] + log_pi[j] );
      llik_contri[i] = log(temp);
    }
  } else {
    for(int i = 0; i < n; i++) {
      llik_contri[i] = ldwnorm2_num(arma::trans(data.row(i)),
                                    par, omega_2pi) - log_c[0];
    }
  }

  return(llik_contri);
}




// [[Rcpp::export]]
double llik_wnorm2_one_comp(arma::mat data, arma::vec par_vec,
                            double log_c, arma::mat omega_2pi)
{
  int n = data.n_rows;
  double log_sum = 0.0;

  for(int i = 0; i < n; i++) {
    log_sum += ldwnorm2_num(arma::trans(data.row(i)), par_vec, omega_2pi);
  }
  log_sum -= n*log_c;

  return log_sum;
}



// // [[Rcpp::export]]
arma::vec grad_log_den_wnorm2_1_comp_1_point(double x, double y,
                                             double s11, double s22, double s12,
                                             double mu1, double mu2,
                                             double det_prec,
                                             arma::mat omega_2pi)
{
  int M = omega_2pi.n_rows;
  double expon_j, term1, term2;
  arma::vec sum_entries = arma::zeros(6);
  for(int j = 0; j < M; j++) {
    term1 = (omega_2pi(j,0) - x + mu1);
    term2 = (omega_2pi(j,1) - y + mu2);
    expon_j = exp(-0.5 *( s11 * term1 * term1 + s22 * term2 * term2 +  2 * s12 * term1 * term2 ));

    sum_entries[0] += expon_j * (s22 - det_prec * term1 * term1 ) * 0.5 / det_prec;
    sum_entries[1] += expon_j * (s11 - det_prec * term2 * term2) * 0.5 / det_prec;
    sum_entries[2] += expon_j * (-s12 - det_prec * term1 * term2 ) / det_prec;
    sum_entries[3] += expon_j * (- s11 * term1 - s12 * term2 );
    sum_entries[4] += expon_j * (- s12 * term1 - s22 * term2 );
    sum_entries[5] += expon_j ;
  }

  // det_prec = det(precision matrix) = (k1k2 -k3^2)

  // Now divide all the 1st 5 terms by the last sum to get grad(log density)
  for(int j = 0; j < 5; j++) {
    sum_entries[j] /= sum_entries[5];
  }

  // finally return log of the last term (log density numerator)
  sum_entries[5] = log(sum_entries[5]);

  return sum_entries;
}






// [[Rcpp::export]]
arma::vec grad_llik_wnorm2_C(arma::mat data, arma::vec par, arma::mat omega_2pi) {
  int n = data.n_rows;
  arma::vec grad_llik_wo_const = arma::zeros(6);
  // first 5 entries = grad, last entry = llik
  double  s11 = par[0], s22 = par[1], s12 = par[2], mu1 = par[3], mu2 = par[4];

  double det_prec = s11*s22 - s12*s12;


  for(int i = 0; i < n; i++) {
    grad_llik_wo_const +=
      grad_log_den_wnorm2_1_comp_1_point(data(i, 0), data(i, 1),
                                         s11, s22, s12, mu1, mu2,
                                         det_prec,
                                         omega_2pi);
  }

  // Now adjust the constant in llik
  grad_llik_wo_const[5] -= n*(log(2*M_PI) - 0.5*log(det_prec));

  return grad_llik_wo_const;
}




// [[Rcpp::export]]
arma::vec grad_den_wnorm2_one_comp_i_unadj(double x, double y, arma::vec par, double det_sig_inv,
                                           double det_sig_inv_sqrt, arma::mat omega_2pi)
{
  double  s11 = par[0], s22 = par[1], s12 = par[2], mu1 = par[3], mu2 = par[4];
  int M = omega_2pi.n_rows;
  double expon_j, term1, term2;
  arma::mat all_entries(6, M);
  for(int j = 0; j < M; j++) {
    term1 = (omega_2pi(j,0) - x + mu1);
    term2 = (omega_2pi(j,1) - y + mu2);
    expon_j = exp(-0.5 *( s11 * term1 * term1 + s22 * term2 * term2 +  2 * s12 * term1 * term2 ));

    all_entries(0, j) = expon_j * (s22 - det_sig_inv * term1 * term1 );// / (4 * M_PI * det_sig_inv_sqrt);
    all_entries(1, j) = expon_j * (s11 - det_sig_inv * term2 * term2);//  / (4 * M_PI * det_sig_inv_sqrt);
    all_entries(2, j) = expon_j * (s12 - det_sig_inv * term1 * term2 );// / (2 * M_PI * det_sig_inv_sqrt);
    all_entries(3, j) = expon_j * (- s11 * term1 - s12 * term2 );// * det_sig_inv_sqrt / (2 * M_PI) ;
    all_entries(4, j) = expon_j * (- s12 * term1 - s22 * term2 );// * det_sig_inv_sqrt / (2 * M_PI) ;
    all_entries(5, j) = expon_j;
  } // multiply by the consts later

  arma::vec sum_entries = arma::zeros(6);
  for(int j = 0; j < M; j++)
    sum_entries += all_entries.col(j);

  return sum_entries;
}


// // not used
// // [[Rcpp::export]]
// arma::mat grad_wnorm2_all_comp(arma::mat data, arma::mat par_mat, arma::vec pi,
//                                arma::mat omega_2pi, int ncores = 1)
// {
//   int n = data.n_rows, K = pi.size(), j;
//   double denom;
//   arma::mat grad_sum = arma::zeros(5, K);
//
//   stdvec_mat grad_sum_part(ncores);
//   for(int g = 0; g < ncores; g++)
//     grad_sum_part[g] = grad_sum;
//
//   arma::vec det_sig_inv(K);
//   for(j = 0; j < K; j++)
//     det_sig_inv[j] = par_mat(0, j) * par_mat(1, j) - par_mat(2, j) * par_mat(2, j);
//
//   arma::vec det_sig_inv_sqrt = arma::sqrt(det_sig_inv);
//
//   for(int i = 0; i < n; i++) {
//     arma::mat grad_den_temp(6, K);
//     int g = get_thread_num_final();
//     // int g = 0;
//     denom = 0;
//     for(j = 0; j < K; j++){
//       // denom += exp(ldwnorm2_num(trans(data.row(i)), par_mat.col(j), omega_2pi) - log_wnorm2_const[j] + log_pi[j] );
//       grad_den_temp.col(j) = grad_den_wnorm2_one_comp_i_unadj(data(i, 0), data(i, 1), par_mat.col(j),
//                         det_sig_inv[j], det_sig_inv_sqrt[j], omega_2pi);
//       denom += grad_den_temp(5, j) * det_sig_inv_sqrt[j] * pi[j];
//     }
//     grad_sum_part[g] += grad_den_temp.head_rows(5) / denom;
//   }  // the pi_h term in numerator hasn't been adjusted yet!
//
//   for(int g = 0; g < ncores; g++)
//     grad_sum += grad_sum_part[g];
//
//   arma::mat adj_factor(5, K);
//   for(j = 0; j < K; j++){
//     adj_factor(0,j) = pi[j] / (2 * det_sig_inv_sqrt[j]);
//     adj_factor(1,j) = pi[j] / (2 * det_sig_inv_sqrt[j]);
//     adj_factor(2,j) = pi[j] / det_sig_inv_sqrt[j];
//     adj_factor(3,j) = pi[j] * det_sig_inv_sqrt[j];
//     adj_factor(4,j) = pi[j] * det_sig_inv_sqrt[j];
//   }
//
//   return (grad_sum % adj_factor);
// }



// // not used
// // [[Rcpp::export]]
// double wnorm2mix(arma::vec x, arma::mat par, arma::vec pi,
//                       arma::vec log_wnorm_const, arma::mat omega_2pi) {
//   double res = 0;
//   int K = pi.size();
//
//   for(int j = 0; j < K; j++)
//     res += exp(ldwnorm2_num(x, par.col(j), omega_2pi) - log_wnorm_const[j]) * pi[j];
//   return(res);
//
// }
//
// // [[Rcpp::export]]
// arma::vec wnorm2mix_manyx(arma::mat x, arma::mat par, arma::vec pi,
//                           arma::vec log_c_von, arma::mat omega_2pi)
// {
//   int n = x.n_rows;
//   arma::vec result(n);
//   for(int i = 0; i < n; i++)
//     result[i] = wnorm2mix(arma::trans(x.row(i)), par, pi, log_c_von, omega_2pi);
//   return(result);
// }


// [[Rcpp::export]]
arma::vec dwnorm2_onex_manypar(arma::vec x, arma::vec k1, arma::vec k2, arma::vec k3,
                               arma::vec mu1, arma::vec mu2, arma::mat omega_2pi)
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

  arma::vec l_const_all = log_const_wnorm2_all(all_par);
  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldwnorm2_num(x, all_par.col(i), omega_2pi);
  }

  return arma::exp(ld_num - l_const_all);
}

// [[Rcpp::export]]
arma::vec dwnorm2_manyx_onepar(arma::mat x, double k1, double k2, double k3,
                               double mu1, double mu2, arma::mat omega_2pi)
{
  int n = x.n_rows;


  arma::vec par(5);
  par << k1 << k2 << k3 << mu1 << mu2 << arma::endr;
  double l_const = l_const_wnorm2(par);
  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldwnorm2_num(arma::trans(x.row(i)), par, omega_2pi);
  }

  return arma::exp(ld_num - l_const);
}

// [[Rcpp::export]]
arma::vec dwnorm2_manyx_manypar(arma::mat x, arma::vec k1, arma::vec k2, arma::vec k3,
                                arma::vec mu1, arma::vec mu2, arma::mat omega_2pi)
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

  arma::vec l_const_all = log_const_wnorm2_all(all_par);

  arma::vec ld_num(n);
  for(int i = 0; i < n; i++) {
    ld_num[i] = ldwnorm2_num(arma::trans(x.row(i)), all_par.col(i), omega_2pi);
  }

  return arma::exp(ld_num - l_const_all);
}


//
// // [[Rcpp::export]]
// List wnorm2_var_corr(double k1, double k2, double k3)
// {
//   double den = k1*k2 - k3*k3,
//     sig1_sq = k2/den,
//     sig2_sq = k1/den,
//     rho = -k3/den;
//
//   double rho_fl = sinh(2*rho*sqrt(sig1_sq*sig2_sq))
//     /sqrt(sinh(2*sig1_sq)* sinh(2*sig2_sq));
//
//   return List::create(Rcpp::Named("var1") = 1-exp(-0.5*sig1_sq),
//                       Rcpp::Named("var2") = 1-exp(-0.5*sig2_sq),
//                       Rcpp::Named("rho") = rho);
//
// }

