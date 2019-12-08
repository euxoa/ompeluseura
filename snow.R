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

model_str2 <- "
data {
   int N;
   real<lower=0> snow[N]; }
transformed data {
   int<lower=0, upper=1> snowp[N];
   for (i in 1:N) snowp[i] = (snow[i] > 0); }
parameters {
   vector[N] innovations;
   real base;
   real<lower=0> sigma;
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
m2 <- stan_model(model_code = model_str2)
fit2 <- sampling(m2, data=stan_data, control=list(adapt_delta=.95))
plot(fit2, par="prob")

christmas %>% filter(!is.na(snow)) %>% 
  mutate(p=apply(extract(fit2, "prob")[[1]], 2, mean)) %>% 
  ggplot(aes(x=year, y=p)) + geom_line() 

plot(fit2, pars=c("ytrend", "sigma"))
change_samples <- extract(fit2, "prob[1]")[[1]] - extract(fit2, "prob[60]")[[1]]
hist(change_samples, n=100)
mean(change_samples<0)
