data {
  int<lower=1> N; // # of dts
  int<lower=1> A; // # of endos
  int<lower=1> B; // # of betas
  matrix[N,A] y;
  matrix[B,A] x[N]; // Exog
  cov_matrix[A] y_cov[N]; //contemporaneous covariance of ys
  int sim; //If this is flagged, dont fit, just generate ys.
}
transformed data {
  real t[N];
  matrix[A,A] L_y_cov[N];
  vector[N] mu = rep_vector(0, N);
  for(n in 1:N) {
    t[n] = n;
    L_y_cov[n] = cholesky_decompose(y_cov[n]);
  }
}
parameters {
  matrix[N,B] beta;
  real<lower=0> alpha[B];
  real<lower=0> sigma[B];
  real<lower=0> length_scale[B];
}
model {
  to_vector(alpha) ~ normal(0.5, 1);
  to_vector(sigma) ~ cauchy(0, 1);
  {
      matrix[N,N] cov[B];
      for(b in 1:B) 
        cov[b] = cov_exp_quad(t, alpha[b], length_scale[b]);
      for(n in 1:N) {
        for(b in 1:B)
          cov[b, n, n] = cov[b, n, n] + square(sigma[b]);
      }
      for(b in 1:B)
        beta[,b] ~ multi_normal_cholesky(mu, cholesky_decompose(cov[b]));
      if(sim==0) {
        for(n in 1:N)
          y[n] ~ multi_normal_cholesky(beta[n] * x[n], L_y_cov[n]);
      }
  }
}

