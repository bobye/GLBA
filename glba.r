############################################################################
## Gated Latent Beta Allocation
## Jianbo Ye <jxy198 at ist.psu.edu>

############################################################################

## compute ratio statistics
ratio = function(I, rate, gamma) {
  J = (1-I);
  diag(J) = 0;
  logratio = log(rate / gamma) %*% I + 
             log((1-rate) / (1-gamma)) %*%  J;
  return (exp(logratio));
}
  
## compute alpha and beta statistics
alpha_and_beta = function(I, ab, tau) {
  omega = as.vector(I %*% tau);
  psi   = rowSums(I);
  return(list(alpha = ab$alpha + omega, 
              beta  = ab$beta + psi - omega));
}

## solve digamma equations
gamma_solve = function(rhs, x0, step=10, rate=0.67, delta=1.) {
  m = nrow(x0);
  x = x0;
  for (i in 1:step) {
    f = digamma(x) - matrix(digamma(rowSums(x)), m, 2) - rhs + 
        delta * matrix(rate * rowSums(x) -  log(rowSums(x)), m, 2);
    df1 =  trigamma(rowSums(x)) + delta / rowSums(x);
    df = cbind(trigamma(x[,1]) - df1 + rate * delta + 1, 
               -df1 + rate * delta + 1, 
               -df1 + rate * delta + 1, 
               trigamma(x[,2]) -df1 + rate * delta + 1);
    invdf = cbind(df[,4], -df[,3], -df[,2], df[,1]) / matrix(df[,1]*df[,4] - df[,2]*df[,3], m, 4);
    x = x - cbind(invdf[,1]*f[,1] + invdf[,2]*f[,2], invdf[,3]*f[,1] + invdf[,4]*f[,2]);
  }
  return(x)
}

## Variational EM update function
varEM_update = function(theta, hypergraph) {
  n = length(hypergraph$I); 
  m = length(theta$tau); 
  Psi = matrix(0, m, 3);
  Tau = as.vector(matrix(0, 1, m));
  Delta = as.vector(matrix(0, 1, m));
  gamma = list(omega = 0, psi = 0);
  for (i in 1:n) {
    I = hypergraph$I[[i]];
    U = hypergraph$inv_oracles[hypergraph$U[[i]]];
    ab = list(alpha = theta$ab[U,1], 
                   beta =  theta$ab[U, 2]);
    tau = theta$tau[U];
    ab_data = alpha_and_beta(I, ab, tau);
    Psi[U,1]=Psi[U,1] + digamma(ab_data$alpha);
    Psi[U,2]=Psi[U,2] + digamma(ab_data$beta);
    Psi[U,3]=Psi[U,3] + digamma(ab_data$alpha + ab_data$beta);
    r_data = ratio(I, ab$alpha / (ab$alpha + ab$beta), theta$gamma);
    Tau[U] = Tau[U] + r_data * tau / (r_data * tau + 1 - tau);
    gamma$omega = gamma$omega + sum(I %*% (1-tau));
    gamma$psi   = gamma$psi   + sum(ncol(I) - rowSums(I) - 1);
    Delta[U] = Delta[U] + 1;
  }
#  aPsi = exp((Psi[,1]-Psi[,3])/Delta);
#  bPsi = exp((Psi[,2]-Psi[,3])/Delta);
#  sPsi = ((aPsi) + (bPsi));
#  theta$ab[,1] = ((0.5*aPsi)/(1-sPsi)+0.5);
#  theta$ab[,2] = ((0.5*bPsi)/(1-sPsi)+0.5);
#  theta$ab = 2*theta$ab / matrix(rep(theta$ab[,1] + theta$ab[,2],2), m, 2);
  theta$ab = gamma_solve(cbind((Psi[,1]-Psi[,3])/Delta, (Psi[,2]-Psi[,3])/Delta), theta$ab, delta = 1/Delta);
  theta$tau = (Tau + 0.7) / (Delta+1);
  theta$gamma = gamma$omega / gamma$psi;
#  plot(theta$tau, ylim=c(0, 1))
#  plot(theta$ab[,1]/(theta$ab[,1]+theta$ab[,2]), ylim=c(0,1))
#  abline(h = theta$gamma, col = "blue", lwd = 2)
  plot((theta$ab[,1]+theta$ab[,2]))
  return(theta)
}


n = length(hypergraph$I);
x = as.data.frame(table(unlist(hypergraph$U)))
hypergraph$oracles = sort(unique (unlist(hypergraph$U)));
hypergraph$oracles = hypergraph$oracles[order(x$Freq, decreasing = TRUE)]
m = length(hypergraph$oracles);
hypergraph$inv_oracles = as.vector(matrix(0, 1, max(hypergraph$oracles)));
hypergraph$inv_oracles[hypergraph$oracles] = 1:m;

theta={};
theta$ab = cbind(matrix(2., m, 1), matrix(1., m, 1));
theta$tau = as.vector(matrix(0.5, m, 1));
theta$gamma = 0.37

for (i in 1:50) {
  theta = varEM_update(theta, hypergraph);
}

