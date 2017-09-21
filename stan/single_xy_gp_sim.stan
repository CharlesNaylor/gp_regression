// multi-dimensional input latent variable GP with single x and y
// cf. manual 2.17.0, pp 253-254, 
// Gaussian Processes for Machine Learning, pp 119-122
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
}
data {
  int<lower=1> N;
  vector[N] x; //Exogs
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] eta;
}
transformed data {
  real delta = 1e-9;
  vector[N] mu = rep_vector(0, N);
  real time[N]; // Time sequence
  for(n in 1:N)
    time[n] = n;
}
parameters {
  vector[N] y;
}
model {
  y ~ normal(eta .* x, sigma);
}
