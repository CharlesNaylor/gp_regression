// multi-dimensional input latent variable GP with single y
// cf. manual 2.17.0, pp 253-254
data {
  int<lower=1> N;
  int<lower=1> D;
  vector[D] x[N]; //Exogs for single currency
  vector[N] y;
}
transformed data {
  real delta = 1e-9;
  vector[D] ones = rep_vector(1,D);
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  matrix[N, D] eta;
}
model {
  vector[N] f;
  
  matrix[N, N] K = cov_exp_quad(x, alpha, rho);
  matrix[N, N] L_K;
  //perturb diagonal elements
  for (n in 1:N) 
    K[n, n] = K[n, n] + delta;
  L_K = cholesky_decompose(K);
  f = (L_K * eta) * ones; //post-multiply by a row of ones to sum by row
    
  rho ~ inv_gamma(5, 5);
  alpha ~  normal(0, 1);
  sigma ~ normal(0, 1);
  to_vector(eta) ~ normal(0, 1);

  y ~ normal(f, sigma);
}
