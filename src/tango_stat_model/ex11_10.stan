data {
  int<lower=1> N;
  int<lower=1> J;
  int mean_day;
  array[N] int<lower=1, upper=J> RatNo;
  vector[N] day;
  vector[N] weight;
}

parameters {
  real mu_a;
  real mu_b;

  real<lower=0> var_a;
  real<lower=0> var_b;
  real<lower=0> var_e;

  vector[J] a_raw;
  vector[J] b_raw;
}

transformed parameters {
  vector[J] alpha = mu_a + sqrt(var_a) * a_raw;
  vector[J] beta = mu_b + sqrt(var_b) * b_raw;
  real born_weight = mu_a - mu_b * mean_day;
  real<lower=0>  sigma_e = sqrt(var_e);
}

model {
  mu_a ~ normal(0, 100);
  mu_b ~ normal(0, 100);
  var_a ~ inv_gamma(0.001, 0.001);
  var_b ~ inv_gamma(0.001, 0.001);
  var_e ~ inv_gamma(0.001, 0.001);

  a_raw ~ normal(0, 1);
  b_raw ~ normal(0, 1);

  for (i in 1:N) {
    real mu = alpha[RatNo[i]] + beta[RatNo[i]] * day[i];
    weight[i] ~ normal(mu, sqrt(var_e));
  }
}
