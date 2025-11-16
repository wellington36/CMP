functions {
  real log_Z_terms(int j, real lambda, real nu){
    return(j * log(lambda) - nu * lgamma(j + 1));
  }
  real ddl_Z_terms(int j, real lambda, real nu){
    return(log(j) + (j - 1) * log(lambda) - nu * lgamma(j + 1));
  }
  real ddv_Z_terms(int j, real lambda, real nu){
    return(j * log(lambda) - log(lgamma(j + 1)) - nu * lgamma(j + 1));
  }
  real d2dv2_Z_terms(int j, real lambda, real nu){
    return(2*log(lgamma(j + 1)) + j * log(lambda) - nu * lgamma(j + 1));
  }
  real d2dl2_Z_terms(int j, real lambda, real nu){
    return(log(j) + log(j - 1) + (j - 2)*log(lambda) - nu * lgamma(j + 1));
  }
  real d2dldv_Z_terms(int j, real lambda, real nu){
    return(log(j) + (j + 1) * log(lambda) - log(lgamma(j + 1)) - nu * lgamma(j + 1));
  }
  // functions for information
}
data {
  int<lower=0> n;   // number of observations 
  int<lower=0> S1;  // sum of X_i's 
  array[n] int<lower=0> X; // data vectors for log_lik construction
  real<lower=0> S2; // sum of log(X_i!) 
}
parameters {
  real<lower=0> lambda;
  real<lower=0> nu;
}
model {
  array[101] real logZ;
  array[100] real ddl_Z;
  array[99] real ddv_Z;
  array[99] real d2dv2_Z;
  array[99] real d2dl2_Z;
  array[99] real d2dldv_Z;
  for (j in 0:100)
    logZ[j+1] = log_Z_terms(j, lambda, nu);
  for (j in 1:100)
    ddl_Z[j] = ddl_Z_terms(j, lambda, nu);
  for (j in 2:100)
    ddv_Z[j-1] = ddv_Z_terms(j, lambda, nu);
  for (j in 2:100)
    d2dv2_Z[j-1] = d2dv2_Z_terms(j, lambda, nu);
  for (j in 2:100)
    d2dl2_Z[j-1] = d2dl2_Z_terms(j, lambda, nu);
  for (j in 2:100)
    d2dldv_Z[j-1] = d2dldv_Z_terms(j, lambda, nu);
  target += S1*log(lambda) - nu*S2 - n*log_sum_exp(logZ) + (0.5)*log(
        ( (n/lambda)*( sum(exp(ddl_Z)) )/( sum(exp(logZ)) ) 
          + n*( ( sum(exp(d2dl2_Z)) * sum(exp(logZ)) 
          - ( sum(exp(ddl_Z)) )^2 )/( sum(exp(logZ))^2 ) )  ) * (
            n*( ( sum(exp(d2dv2_Z))*sum(exp(logZ)) - ( sum(exp(ddv_Z)) )^2 )/( sum(exp(logZ))^2 ) )
          ) - (
            ( n*( sum(exp(d2dldv_Z))*sum(exp(logZ)) - sum(exp(ddl_Z))*sum(exp(ddv_Z)) )/( sum(exp(logZ))^2 ) )^2
          )
     );
}
generated quantities {
  array[101] real logZ;
  array[n] real log_lik;
  for (j in 0:100)
    logZ[j+1] = log_Z_terms(j, lambda, nu);
  for (i in 1:n)
    log_lik[i] = X[i]*log(lambda) - nu*lgamma(X[n] + 1) - n*log_sum_exp(logZ);
}
