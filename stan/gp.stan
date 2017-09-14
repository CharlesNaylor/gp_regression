data {
  int<lower=1> N; // # of dts
  int<lower=1> B; // # of betas
  matrix[N,K] y;
  matrix[B,K] x[N]; // Exog
  cov_matrix[K] y_cov; // contemporaneous covariance of ys
  int sim; // If this is flagged, dont fit, just generate ys.
}
transformed data {
  real t[N];
  matrix[K,K] L_y_cov;
  t[1] = 1;
  for(n in 2:N)
    t[n] = t[n-1]+1;
  L_y_cov = cholesky_decompose(y_cov);
}
parameters {
  matrix[N,B] beta;
  vector<lower=0>[B] alpha;
  real<lower=0> sigma[B];
  real<lower=0> rho[B];
  cholesky_factor_corr[B] L_Omega;
}
model {
  matrix[N,B] f; // factor effects
  {
      matrix[N,N] K = cov_exp_quad(x, 1.0, rho);
      matrix[N,N] L_K;
      
      // diagonal elements
      for(n in 1:N) { 
        K[n,n] = K[n,n] + delta;
     
        L_K = cholesky_decompose(K);
        f = L_K * beta * diag_pre_multiply(alpha, L_Omega)';
      }
  }
  alpha ~ normal(0.25, 1);
  sigma ~ normal(0, 1);
  to_vector(beta) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(3);
  rho ~ inv_gamma(5, 5);
  if(sim==0) {
      to_vector(y) ~ normal(to_vector(f), sigma);
  }
}
generated quantities {
  matrix[N,K] y_sim;
  for(n in 1:N)
    y_sim[n] = to_row_vector(multi_normal_cholesky_rng(to_vector(beta[n]*x[n]),
                          L_y_cov));
}
