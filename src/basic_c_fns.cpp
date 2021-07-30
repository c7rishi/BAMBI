#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "bessel.h"

// [[Rcpp::export]]
arma::vec rowVars(arma::mat mat_in)
{
  int nrow = mat_in.n_rows;
  arma::vec res(nrow);
  for(int i = 0; i < nrow; i++) {
    res[i] = arma::var(mat_in.row(i));
  }
  return res;
}



// perm_mat: n.pars x n.comp x n.iter
// perm.lab: n.iter x n.comp
// [[Rcpp::export]]
arma::cube par_mat_permute(arma::cube par_mat,
                           arma::umat perm_lab)
{
  int n_iter = par_mat.n_slices, n_row = par_mat.n_rows, n_col = par_mat.n_cols;
  arma::cube result(n_row, n_col, n_iter);
  for(int iter = 0; iter < n_iter; iter++) {
    for(int row = 0; row < n_row; row++) {
      for(int col = 0; col < n_col; col++) {
        result(row, col, iter) = par_mat(row, perm_lab(iter, col) - 1, iter);
      }
    }
  }

  return result;
}



// // [[Rcpp::export]]
// arma::cube comp_ind_permute(arma::cube clus_ind,
//                            arma::umat perm_lab)
// {
//   int n_iter = par_mat.n_slices, n_row = par_mat.n_rows, n_col = par_mat.n_cols;
//   arma::cube result(n_row, n_col, n_iter);
//   for(int iter = 0; iter < n_iter; iter++) {
//     for(int row = 0; row < n_row; row++) {
//       for(int col = 0; col < n_col; col++) {
//         result(row, col, iter) = par_mat(row, perm_lab(iter, col) - 1, iter);
//       }
//     }
//   }
//
//   return result;
// }



// [[Rcpp::export]]
Rcpp::NumericVector cID(arma::mat probs, int ncomp, arma::vec Uv) {
  double U;
  double* p = new double[ncomp];
  Rcpp::NumericVector clID(probs.n_rows);

  for (int i = 0; i < (int)probs.n_rows; i++) {
    U = Uv[i];
    //Rcout << U;
    p[0] = probs(i,0);
    if (U < p[0]) {
      clID[i] = 1;
      continue;
    }

    for (int j = 1; j < ncomp; j++) {
      p[j] = probs(i,j) + p[j-1];
      if (U < p[j]) {
        clID[i] = j+1;
        break;
      }
    }

  }

  delete[] p;
  return(clID);
}



// [[Rcpp::export]]
arma::uvec change_labs(arma::uvec orig, arma::uvec rand_perm) {
 unsigned int n = orig.n_elem;
 arma::uvec out(n);
  for(unsigned int i = 0; i < n; i++) {
    out[i] = rand_perm[orig[i]-1];
  }

  return out;
}


// int sgn(double val) {
//   if (val > 0) {
//     return 1;
//   } else if (val < 0) {
//     return -1;
//   } else {
//     return 0;
//   }
// }


// [[Rcpp::export]]
double calc_corr_tau_2(arma::mat samp_mat)
{
  int N = samp_mat.n_rows;
  double theta_ij, phi_ij;
  double num = 0;
  for(int i = 0; i < N-1; i++) {
    for(int j = i+1; j < N; j++) {
      theta_ij = samp_mat(i, 0) - samp_mat(j, 0);
      phi_ij = samp_mat(i, 1) - samp_mat(j, 1);
      if (theta_ij >= 0) {
        if (phi_ij >= 0) {
          num += sgn((theta_ij - M_PI)*(phi_ij - M_PI));
        } else {
          num += sgn((theta_ij - M_PI)*(phi_ij + M_PI));
        }
      } else {
        if (phi_ij >= 0) {
          num += sgn((theta_ij + M_PI)*(phi_ij - M_PI));
        } else {
          num += sgn((theta_ij + M_PI)*(phi_ij + M_PI));
        }
      }


    }
  }
  return 2 * num/(N * (N-1));
}


// [[Rcpp::export]]
double calc_corr_tau_1(arma::mat samp_mat)
{
  int N = samp_mat.n_rows;
  double sum_delta_ij = 0;
  for(int i = 0; i < N-2; i++) {
    for(int j = i+1; j < N-1; j++) {
      for(int k = j+1; k < N; k++) {
        sum_delta_ij +=
          sgn(samp_mat(i, 0)-samp_mat(j, 0)) *
          sgn(samp_mat(j, 0)-samp_mat(k, 0)) *
          sgn(samp_mat(k, 0)-samp_mat(i, 0))
        *
          sgn(samp_mat(i, 1)-samp_mat(j, 1)) *
          sgn(samp_mat(j, 1)-samp_mat(k, 1)) *
          sgn(samp_mat(k, 1)-samp_mat(i, 1));
      }
    }
  }
  return sum_delta_ij * 6.0 / (N*(N-1)*(N-2));
}


// [[Rcpp::export]]
double calc_corr_fl(arma::mat samp_mat)
{
  int N = samp_mat.n_rows;
  double sum_1 = 0, sum_2 = 0, sum_prod = 0, temp_1, temp_2;
  for(int i = 0; i < N-1; i++) {
    for(int j = i+1; j < N; j++) {
      temp_1 = sin(samp_mat(i, 0) - samp_mat(j, 0));
      temp_2 = sin(samp_mat(i, 1) - samp_mat(j, 1));
      sum_1 += temp_1*temp_1;
      sum_2 += temp_2*temp_2;
      sum_prod += temp_1*temp_2;
    }
  }
  return sum_prod/sqrt(sum_1*sum_2);
}



// // [[Rcpp::export]]
// double calc_corr_js(arma::mat samp_mat)
// {
//   int N = samp_mat.n_rows;
//   double sum_1 = 0, sum_2 = 0, sum_prod = 0, temp_1, temp_2;
//   for(int i = 0; i < N; i++) {
//       temp_1 = sin(samp_mat(i, 0));
//       temp_2 = sin(samp_mat(i, 1));
//       sum_1 += temp_1*temp_1;
//       sum_2 += temp_2*temp_2;
//       sum_prod += temp_1*temp_2;
//     }
//   return sum_prod/sqrt(sum_1*sum_2);
// }
