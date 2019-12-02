data {
  int<lower=0> N;
  int<lower=1, upper=5> y[N];
}
parameters {
  real mu;
  real<lower=0> spread;
}
transformed parameters {
  vector[4] cutoffs;
  cutoffs[1] = mu - 3*spread;
  cutoffs[2] = mu - spread;
  cutoffs[3] = mu + spread;
  cutoffs[4] = mu + 3*spread; 

}
model {
  for (i in 1:N) y[i] ~ ordered_logistic(0.0, cutoffs);
}
generated quantities {
  vector[5] probs;
  for (i in 1:5) probs[i] = exp(ordered_logistic_lpmf(i | 0.0, cutoffs));
}

