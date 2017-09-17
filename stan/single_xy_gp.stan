// multi-dimensional input latent variable GP with single y
// cf. manual 2.17.0, pp 253-254, 
// Gaussian Processes for Machine Learning, pp 119-122
data {
  int<lower=1> N;
  vector[N] x; //Exogs
  vector[N] y;
}
transformed data {
  real delta = 1e-9;
  real time[N]; // Time sequence
  vector[N] mu = rep_vector(0, N);
  for(n in 1:N)
    time[n] = n;
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] eta;
}
model {
  matrix[N, N] K = cov_exp_quad(time, alpha, rho);
  matrix[N, N] L_K;
  //perturb diagonal elements
  for (n in 1:N) 
    K[n, n] = K[n, n] + delta;
  L_K = cholesky_decompose(K);
    
  rho ~ inv_gamma(5, 5);
  alpha ~  normal(0.5, 1);
  sigma ~ normal(0, 1);
  eta ~ multi_normal_cholesky(mu, L_K);

  y ~ normal(eta .* x, sigma);
}

