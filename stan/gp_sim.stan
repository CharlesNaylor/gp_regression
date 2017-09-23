// Full gp model, in which beta is matrix-variate
// wound up implementing the beta variance via 
// hierarchy, as in Trangucci's work.
functions {
  matrix eta_rng(int N,int B,vector alpha, vector length_scale,
                vector sigma, matrix L_omega) {
    matrix[N,B] beta_mu; //Instead of 0s, I'm going to 
    //put the beta-wise covariance on this term to avoid
    //having a massive kronecker product matrix.
    matrix[N,N] cov[B];
    matrix[B,B] L_Sigma_b;
    matrix[N,B] beta;
    real delta = 1e-9;
    real t[N];
    for(n in 1:N)
      t[n] = n;
    L_Sigma_b = diag_pre_multiply(sigma, L_omega);
    for(b in 1:B) 
      cov[b] = cov_exp_quad(t, alpha[b], length_scale[b]);
    for(n in 1:N) {
      beta_mu[n] = to_row_vector(multi_normal_cholesky_rng(rep_vector(0,B), 
                                                        L_Sigma_b));
      for(b in 1:B)
        cov[b, n, n] = cov[b, n, n] + delta;
    }
    for(b in 1:B)
      beta[,b] = multi_normal_cholesky_rng(beta_mu[,b], 
                            cholesky_decompose(cov[b]));
    return(beta);
  }
  matrix rho_rng(int B, real scale) {
    return(lkj_corr_cholesky_rng(B, scale));
  }
}
data {
  int<lower=1> N; // # of dts
  int<lower=1> A; // # of endos
  int<lower=1> B; // # of betas
  matrix[B,A] x[N]; // Exog
  cov_matrix[A] y_cov[N]; //contemporaneous covariance of ys
  matrix[N,B] beta;
  real<lower=0> alpha[B];
  vector<lower=0>[B] sigma;
  real<lower=0> length_scale[B];
  cholesky_factor_corr[B] L_omega; //eta correlation
}
transformed data {
  matrix[A,A] L_y_cov[N];
  for(n in 1:N)
    L_y_cov[n] = cholesky_decompose(y_cov[n]);
}
parameters {}
model {}
generated quantities {
  vector[A] y[N];
  for(n in 1:N)
    y[n] = multi_normal_cholesky_rng(to_vector(beta[n] * x[n]),
                                     L_y_cov[n]);
}

