//latent variable GP with single x and multiple ys
functions {
  vector eta_rng(int N, real alpha, real rho) {
    matrix[N,N] K;
    matrix[N,N] L_K;
    real time[N];

    for(n in 1:N)
      time[n] = n;
    K = cov_exp_quad(time, alpha, rho);
    for(n in 1:N)
      K[n,n] = K[n,n] + 1e-9;
    L_K = cholesky_decompose(K);
    return multi_normal_cholesky_rng(rep_vector(0,N), L_K);
  }
  matrix L_Omega_rng(int A, real eta) {
    return lkj_corr_cholesky_rng(A, eta);
  }
}
data {
  int<lower=1> N;
  int<lower=1> A; // # of assets
  matrix[N,A] x; //Exogs
  real<lower=0> rho;
  real<lower=0> alpha;
  vector[N] eta;
  matrix[A,A] Sigma_y; // covariance matrix for endos
}
transformed data {
  real delta = 1e-9;
  vector[N] mu = rep_vector(0, N);
  real time[N]; // Time sequence
  matrix[A,A] L_y = cholesky_decompose(Sigma_y);
  row_vector[A] mu_y[N]; // I'd like to do diag_pre_multiply(eta, x), but
  // that makes a matrix, which, AFAIK, I would still have to loop through
  // to convert to an array of row_vectors to provide to multi_normal_chol
  for(n in 1:N) {
    time[n] = n;
    mu_y[n] = eta[n] * x[n];
  }
}
parameters {
  row_vector[A] y[N]; 
}
model {
  y ~ multi_normal_cholesky(mu_y, L_y);
}
