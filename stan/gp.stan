// Full gp model, in which beta is matrix-variate
// wound up implementing the beta variance via 
// hierarchy, as in Trangucci's work.
data {
  int<lower=1> N; // # of dts
  int<lower=1> A; // # of endos
  int<lower=1> B; // # of betas
  matrix[N,A] y;
  matrix[B,A] x[N]; // Exog
  cov_matrix[A] y_cov[N]; //contemporaneous covariance of ys
}
transformed data {
  real t[N];
  matrix[A,A] L_y_cov[N];
  real delta = 1e-9;
  vector[B] beta_mu_mu = rep_vector(0, B);
  for(n in 1:N) {
    t[n] = n;
    L_y_cov[n] = cholesky_decompose(y_cov[n]);
  }
}
parameters {
  matrix[N,B] beta;
  real<lower=0> alpha[B];
  vector<lower=0>[B] sigma;
  real<lower=0> length_scale[B];
  cholesky_factor_corr[B] L_omega; //eta correlation
  matrix[N,B] beta_mu; //Instead of 0s, I'm going to 
  //put the beta-wise covariance on this term to avoid
  //having a massive kronecker product matrix.
}
model {
  matrix[N,N] cov[B];
  matrix[B,B] L_Sigma_b;

  to_vector(alpha) ~ normal(0.5, 1);
  sigma ~ cauchy(0, 1);
  L_omega ~ lkj_corr_cholesky(4);

  L_Sigma_b = diag_pre_multiply(sigma, L_omega);
  for(b in 1:B) 
    cov[b] = cov_exp_quad(t, alpha[b], length_scale[b]);
  for(n in 1:N) {
    beta_mu[n] ~ multi_normal_cholesky(beta_mu_mu, L_Sigma_b);
    for(b in 1:B)
      cov[b, n, n] = cov[b, n, n] + delta;
  }
  for(b in 1:B)
    beta[,b] ~ multi_normal_cholesky(beta_mu[,b], cholesky_decompose(cov[b]));
  for(n in 1:N)
    y[n] ~ multi_normal_cholesky(beta[n] * x[n], L_y_cov[n]);
}

