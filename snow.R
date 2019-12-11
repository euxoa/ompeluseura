library(dplyr)
library(rstan)
library(readr)

d <-  read_csv("https://raw.githubusercontent.com/dins/snow-depth/master/kaisaniemi.csv", 
               col_types = cols(`Lumensyvyys (cm)` = col_character())) %>%
      setNames(c("year", "month", "day", "clock", "tzone", "snow", "temp")) %>%
      mutate(date = ISOdate(year, month, day), 
             snow = ifelse(snow==FALSE, 0, as.numeric(snow)), 
             snow = ifelse(snow < 0, 0, snow)) %>%
      select(date, year, month, day, snow, temp)
saveRDS(d, file="kaisaniemi_daily.rds")
                  
christmas <- d %>% filter(day==24 & month==12)
readr::write_csv(christmas, "kaisaniemi_christmas.csv")

stan_data <- with(christmas %>% filter(!is.na(snow)), list(decade=(year-2000)/10, snow=snow, N=length(snow)))


first_model_code <- "
data {
   int N;
   int<lower=0, upper=1> is_snow[N]; }
parameters {
   real b; }
model {
   for (i in 1:N) is_snow[i] ~ bernoulli_logit(b); }
"


model_str <- "
data {
   int N;
   real<lower=0> snow[N];
   real decade[N]; }
transformed data {
   int<lower=0, upper=1> snowp[N];
   for (i in 1:N) snowp[i] = (snow[i] > 0); }
parameters {
   real b;
   real k; }
model {
   for (i in 1:N) snowp[i] ~ bernoulli_logit(k * decade[i] + b); }
generated quantities {
   real prob[N];
   for (i in 1:N) prob[i] = inv_logit(k* decade[i] + b);
}
"

m <- stan_model(model_code = model_str)
fit <- sampling(m, data=stan_data)
fit
hist(posterior_samples(fit, "k")$k, n=100)
plot(fit, pars=c("k", "b"))
plot(fit, pars="prob")

# Vectorized state-space model

model_str_statespace <- "
data {
   int N;
   real<lower=0> snow[N]; }
transformed data {
   int<lower=0, upper=1> snowp[N];
   for (i in 1:N) snowp[i] = (snow[i] > 0); }
parameters {
   vector[N] innovations;
   real base;
   real<lower=0> sigma_change;
   real ytrend;
   }
model {
   innovations ~ normal(ytrend, 1);
   snowp ~ bernoulli_logit(base + sigma * cumulative_sum(innovations));
   base ~ normal(0, 5);
   sigma ~ normal(0, .1);
   ytrend ~ normal(0, .1);
}
generated quantities {
   vector[N] prob;
   prob = inv_logit(base + sigma * cumulative_sum(innovations));
}
"
mss <- stan_model(model_code = model_str_statespace)
fit_ss <- sampling(m2, data=stan_data, control=list(adapt_delta=.95))
plot(fit_ss, par="prob")

christmas %>% filter(!is.na(snow)) %>% 
  mutate(p=apply(extract(fit_ss, "prob")[[1]], 2, mean)) %>% 
  ggplot(aes(x=year, y=p)) + geom_line() 

apply(extract(fit_ss, "prob")[[1]], 2, function (x) quantile(x, c(.25, .5, .75))) %>% 
  t() %>% as.data.frame() %>% setNames(c("pmin", "p", "pmax")) %>% 
  mutate(year=stan_data$decade) %>%
  ggplot(aes(x=year, ymin=pmin, y=p, ymax=pmax)) + geom_ribbon(alpha=.5) + geom_line()


plot(fit_ss, pars=c("ytrend", "sigma"))
ys <- data.frame(extract(fit_ss, pars=c("ytrend", "sigma")))
ys %>% ggplot(aes(x=ytrend, y=sigma)) + geom_point(alpha=.2)
cor(ys)
change_samples <- extract(fit_ss, "prob[1]")[[1]] - extract(fit_ss, "prob[60]")[[1]]
hist(change_samples, n=100)
mean(change_samples<0)
