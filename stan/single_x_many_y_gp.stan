//latent variable GP with single x and multiple ys
data {
  int<lower=1> N;
  int<lower=1> A; // # of assets
  matrix[N,A] x; //Exogs
  row_vector[A] y[N];
}
transformed data {
  real delta = 1e-9;
  vector[N] mu = rep_vector(0, N);
  real time[N]; // Time sequence
  for(n in 1:N)
    time[n] = n;
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  vector[N] eta;
  cholesky_factor_corr[A] L_omega; //y correlation
  vector<lower=0>[A] sigma_y; //y scale
}
transformed parameters {
  row_vector[A] mu_y[N];
  for(n in 1:N) {
    mu_y[n] = eta[n] * x[n];
  }
}
model {
  matrix[A, A] L_Sigma_y;
  matrix[N, N] K = cov_exp_quad(time, alpha, rho);
  matrix[N, N] L_K;
  //perturb diagonal elements
  for (n in 1:N) 
    K[n, n] = K[n, n] + delta;
  L_K = cholesky_decompose(K);
  L_Sigma_y = diag_pre_multiply(sigma_y, L_omega);

  rho ~ inv_gamma(5, 5);
  alpha ~ normal(0.5, 1);
  L_omega ~ lkj_corr_cholesky(5);
  sigma_y ~ cauchy(0, 2.5);
  eta ~ multi_normal_cholesky(mu, L_K);
  y ~ multi_normal_cholesky(mu_y, L_Sigma_y);
}
