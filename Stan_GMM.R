model = "
data {
  int<lower=1> N;
  int<lower=1> k;
  real X[N];
}

parameters {
  simplex[k] theta[N];
  real<lower=0> mu[k];
  real<lower=0> sdc[k];
}

model {
  real ps[k];
  real ps2[k];
  real t[N];
  for(i in 1:N){
    for(j in 1:k){
      ps[j] = log(theta[i][j]) + normal_log(X[i], mu[j], sdc[j]);
    }
    target += log_sum_exp(ps);
  }
}
"
X = c(rnorm(100,15,2),rnorm(200,2,0.5),rnorm(300,5,2))
dat <- list(N = 600, k=3, X=X)
fit <- stan(model_code = model, model_name = "MM", 
            data = dat, iter = 1000, chains = 1, sample_file = 'GMM.csv',
            verbose = TRUE)
gmm = extract(fit)
mixtures_map = cbind(X,apply(cbind(apply(gmm$theta[,,1],2,median),apply(gmm$theta[,,2],2,median),apply(gmm$theta[,,3],2,median)),1,function(x){which(x==max(x))}))
plot(mixtures_map[,1])
points(mixtures_map[,1], col=c("red","blue","green")[mixtures_map[,2]])