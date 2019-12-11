library(dplyr)
library(rstan)
library(brms)
# yksinkertainen esimerkki, normaali tai t-jakauma tai molemmat
# schools tms hierarkinen
# aikasarja, esim. bal

d <- readr::read_csv("kilpisjarvi_raw.csv") %>% 
  setNames(c("vuosi", "kk", "pv", "klo", "tzone", "temp")) %>%
  mutate(t = vuosi + (kk-1)/12 - 1979, kkf = as.factor(kk), kymmen=(vuosi-1999)/10)

d %>% ggplot(aes(x=t, y=temp, color=kkf)) + geom_line()
d %>% ggplot(aes(x=t, y=temp, color=kkf)) + geom_smooth(method="gam")



m <- brm(bf(temp ~  kymmen + (1+kymmen|kkf), sigma ~ (1|kkf)), autocor=cor_ar(~ t, p=1), data=d, chains=2)

plot(m)
plot(ranef(m)$kkf[,,"kymmen"][,"Estimate"])
plot(ranef(m)$kkf[,,"Intercept"][,"Estimate"])
plot(ranef(m)$kkf[,,"sigma_Intercept"][,"Estimate"])
posterior_samples(m, pars="r_kkf\\[.*,kymmen\\]")

m <- stan_model("stars.stan")
fit <- sampling(m, data=list(y=6-c(1, 2, 2, 2, 2, 1, 2, 3, 2, 2, 3), N=11), chains=1)
plot(fit, par="probs")

## here we go

d <- readr::read_csv("kilpisjarvi_raw.csv") %>% 
  setNames(c("year", "month", "day", "_clock", "tzone", "temp")) %>%
  mutate(t = ISOdate(year, month, day), 
         f_month = as.factor(month),
         decade = as.numeric(t - ISOdate(2000, 1, 1), units="days")/365.25) %>%
  select(year, f_month, t, decade, temp)

ggplot(d, aes(x=t, y=temp)) + geom_line()
ggplot(d, aes(x=t, y=temp, color=f_month)) + geom_line()
ggplot(d, aes(x=t, y=temp, color=f_month)) + geom_smooth(method="lm")

d_year <- d %>% group_by(year) %>% summarise(n=n(), temp=mean(temp), decade=mean(decade)) %>% filter(n==12)

ggplot(d_year, aes(x=year, y=temp)) + geom_line()

model_str1 <- "
data {
  int N;
  real decade[N];
  real temp[N];
}
parameters {
  real<lower=0> sigma;
  real b;
  real k;
}
model {
  for (i in 1:N) temp[i] ~ normal(k*decade[i] + b, sigma);
}
"
m <- stan_model(model_code = model_str1)
fit <- sampling(m, data = with(d_year, list(N=length(temp), decade=decade, temp=temp)))

traceplot(fit)
plot(fit)
fit

model_str2 <- "
data {
  int N;
  real decade[N];
  real temp[N];
  int month[N];
}
parameters {
  real<lower=0> sigma;
  real b[12];
  real k[12];
}
model {
  for (i in 1:N) temp[i] ~ normal(k[month[i]]*decade[i] + b[month[i]], sigma);
  sigma ~ normal(0, 5);
  b ~ normal(0, 5);
  k ~ normal(0, 1);
}
"
m2 <- stan_model(model_code = model_str2)
fit2 <- sampling(m2, data = with(d, list(N=length(temp), decade=decade, month=as.integer(f_month), temp=temp)))

model_str3 <- "
data {
  int N;
  real decade[N];
  real temp[N];
  int month[N];
}
parameters {
  real<lower=0> sigma[12];
  real b[12];
  real k[12];
}
model {
  for (i in 1:N) {
     int m = month[i];
     temp[i] ~ normal(k[m] * decade[i] + b[m], sigma[m]); }
  sigma ~ normal(0, 5);
  b ~ normal(0, 5);
  k ~ normal(0, 1);
}
"
m3 <- stan_model(model_code = model_str3)
fit3 <- sampling(m3, data = with(d, list(N=length(temp), decade=decade, month=as.integer(f_month), temp=temp)))
plot(fit3, par="sigma")


model_str4 <- "
data {
  int N;
  real decade[N];
  real temp[N];
  int month[N];
}
parameters {
  real lsigma[12];
  real b[12];
  real k[12];
  real<lower=0> sigma_lsigma;
  real<lower=0> sigma_b;
  real<lower=0> sigma_k;
}
model {
  for (i in 1:N) {
     int m = month[i];
     temp[i] ~ normal(k[m] * decade[i] + b[m], exp(lsigma[m])); }
  lsigma ~ normal(0, sigma_lsigma);
  b ~ normal(0, sigma_b);
  k ~ normal(0, sigma_k);
  sigma_lsigma ~ normal(0, 2);
  sigma_b ~ normal(0, 5);
  sigma_k ~ normal(0, 5);
}
generated quantities {
  real sigma[12];
  for (m in 1:12) sigma[m] = exp(lsigma[m]);
}
"
m4 <- stan_model(model_code = model_str4)
fit4 <- sampling(m4, data = with(d, list(N=length(temp), decade=decade, month=as.integer(f_month), temp=temp)))

plot(fit4, par="k")
plot(fit3, par="k")

plot(fit4, par="sigma")
plot(fit3, par="sigma")
plot(fit4, par="sigma_lsigma")

hist(extract(fit4, "k[1]")[[1]] - extract(fit4, "k[6]")[[1]], n=100)
hist(extract(fit4, "b[1]")[[1]] - extract(fit4, "b[12]")[[1]], n=100)

model_str5 <- "
data {
   int N;
   real temp[N];
}
parameters {
   real tlevel[N];
   real<lower=0> obs_sigma;
   real<lower=9> lsigma;
}
model {
   for (i in 2:N) 
      tlevel[i] ~ normal(tlevel[i-1], lsigma);
   tlevel[1] ~ normal(-2, 5);
   for (i in 1:N)
       temp[i] ~ normal(tlevel[i], obs_sigma);
   obs_sigma ~ normal(0, 5);
   lsigma ~ normal(0, .5);
}"

m5 <- stan_model(model_code = model_str5)
fit5 <- sampling(m5, data = with(d_year, list(N=length(temp), temp=temp)))

