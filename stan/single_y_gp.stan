// multi-dimensional input latent variable GP with single y
// cf. manual 2.17.0, pp 253-254, 
// Gaussian Processes for Machine Learning, pp 119-122
data {
  int<lower=1> N;
  int<lower=1> D;
  matrix[N, D] x; //Exogs
  vector[N] y;
}
transformed data {
  real delta = 1e-9;
  vector[D] ones = rep_vector(1,D);
  real time[N]; // Time sequence
  for(n in 1:N)
    time[n] = n;
}
parameters {
  real<lower=0> rho;
  vector<lower=0>[D] alpha;
  real<lower=0> sigma;
  matrix[N, D] eta;
  cholesky_factor_corr[D] L_Omega;
}
model {
  vector[N] f;
  
  matrix[N, N] K = cov_exp_quad(time, 1.0, rho);
  matrix[N, N] L_K;
  //perturb diagonal elements
  for (n in 1:N) 
    K[n, n] = K[n, n] + delta;
  L_K = cholesky_decompose(K);
  f = L_K * ((eta * diag_pre_multiply(alpha, L_Omega)' .* x) * ones); //post-multiply by a vector of ones to sum by row
    
  rho ~ inv_gamma(5, 5);
  alpha ~  normal(0.5, 1);
  sigma ~ normal(0, 1);
  to_vector(eta) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(3);

  y ~ normal(f, sigma);
}

