data {
  int<lower=1> D; // Number of documents
  int<lower=1> V; // Vocabulary size
  int<lower=1> K; // Number of topics
  int<lower=1> N; // Total word instances
  int<lower=1,upper=V> words[N]; // Word index for each word instance
  int<lower=1,upper=D> docs[N]; // Document index for each word instance
}

parameters {
  simplex[K] theta[D]; // Topic distribution for document D, influenced by spatial proximity
  simplex[V] phi[K]; // Word distribution for topic K
}

model {
  
  for (k in 1:K) {
    phi[k] ~ dirichlet(rep_vector(0.01, V)); // Prior on word distributions in each topic
  }
  
  
  for (d in 1:D) {
    theta[d] ~ dirichlet(rep_vector(0.01, K)); 
  }
  
  for (n in 1:N) {
    real log_prob[K];
    for(k in 1:K){
      log_prob[k] = log(theta[docs[n],k]) + log(phi[k,words[n]]);
    }
    target += log_sum_exp(log_prob);    
  }
      
  
}