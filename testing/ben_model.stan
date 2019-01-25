data {
    int          N;    //number of observations training and validation set
    int          P;    //number of predictors
    real         y[N]; //response vector
    matrix[N, P] X;    //model matrix
}

parameters {
    real            mu;       //intercept
    real<lower = 0> sigma2;   //error variance
    vector[P]       beta_raw; // regression parameters

    //hyperparameters prior
    vector<lower = 1>[P] tau;
    real<lower = 0>      lambda1;
    real<lower = 0>      lambda2;
}

transformed parameters {
    vector[P]       beta;
    real<lower = 0> sigma;   //error sd
    vector[N]       linpred; //mean normal model

    for(j in 1 : P){
    	beta[j] = sqrt(((sigma2 * (tau[j] - 1))/(lambda2 * tau[j]))) * beta_raw[j];
    }
    
    sigma   = sqrt(sigma2);
    linpred = mu + X * beta;
}

model {
    //prior regression coefficients
    beta_raw ~ normal(0, 1); //implies beta ~ normal(0, sqrt(1/(lambda2/sigma2 * (tau[j]/(tau[j]-1)))))
    tau      ~ gamma(0.5, (8 * lambda2 * sigma2) / lambda1^2);

    lambda1 ~ cauchy(0, 1);
    lambda2 ~ cauchy(0, 1);
	
   //priors nuisance parameters: uniform on log(sigma^2) & mu
   target += -2 * log(sigma); 
	
   //likelihood
   y ~ normal(linpred, sigma);
}
