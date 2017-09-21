# Simulate fake data then try to recover it with
# a latent GP regression with a single x and y.

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(8675309) #Let's keep the randomness consistent
fake_data = list(N=100, D=3, alpha=1.0, rho=5, sigma=0.05,
                 x=rnorm(N))

sim_model <- stan_model(file="single_xy_gp_sim.stan", 
                        model_name="sim")
rstan::expose_stan_functions(sim_model)
fake_data$eta <- with(fake_data, eta_rng(N, alpha, rho))

sim_samples <- sampling(sim_model, data=fake_data,
                        chains=2)
sim_y <- apply(rstan::extract(sim_samples)$y,2,median)

fit <- stan(file="single_xy_gp.stan", model_name="fit_y",
            data=list(y=sim_y,N=fake_data$N,x=fake_data$x))
print(density(extract(fit)$alpha))
