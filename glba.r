# Copyright (C) 2016 The Pennsylvania State University, USA - All Rights Reserved
#
# Unauthorized copying of this file, via any medium is strictly prohibited
#
#  A non-exclusive, nontransferable, perpetual license is granted to
#  you to install and use the Software for academic, non-profit, or
#  government-sponsored research purposes. Use of the software under
#  this license is restricted to non-commercial purposes. COMMERCIAL
#  USE OF THE SOFTWARE REQUIRES A SEPARATELY EXECUTED WRITTEN LICENSE
#  AGREEMENT.
#
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
#  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#  PARTICULAR PURPOSE ARE DISCLAIMED. THE COPYRIGHT OWNER MAKES NO
#  REPRESENTATION OR WARRANTY THAT THE SOFTWARE WILL NOT INFRINGE ANY
#  PATENT OR OTHER PROPRIETARY RIGHT. IN NO EVENT SHALL THE COPYRIGHT OWNER BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
#  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
#  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
#  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
#  DAMAGE.
#
# NOTICE: All information contained herein is, and remains the property of The
# Pennsylvania State University. The intellectual and technical concepts
# contained herein are proprietary to The Pennsylvania State University and may
# be covered by U.S. and Foreign Patents, patents in process, and are protected
# by trade secret or copyright law. Dissemination of this information or
# reproduction of this material is strictly forbidden unless prior written
# permission is obtained from The Pennsylvania State University. If
# you obtained this code from other sources, please write to Jianbo Ye or James Z. Wang.
#
# Jianbo Ye, Jia Li, Michelle G. Newman, Reginald B. Adams, Jr. and James Z. Wang,
# ``Probabilistic Multigraph Modeling for Improving the Quality of Crowdsourced Affective Data,''
# IEEE Transactions on Affective Computing, vol. , no. , 14 pages, 2017.
#
# Written by Jianbo Ye <jxy198@psu.edu>, 2016


############################################################################
## Gated Latent Beta Allocation
## Jianbo Ye <jxy198 at ist.psu.edu>

############################################################################
## define prior
tau0 = 0.5
rate0 = 0.7
gamma0 = 0.3 # fixed prior, not critical
max_gamma = 0.5

## compute ratio statistics
ratio = function(I, alpha, beta, gamma, weight =1.) {
  J = (1-I) * weight;
  I = I * weight;
  diag(J) = 0;
  logratio = log(alpha / (alpha+beta) / gamma) %*% I + 
             log(beta / (alpha+beta) / (1-gamma)) %*%  J;
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
gamma_solve = function(rhs, x0, step=10, rate=rate0, delta=1.) {
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
    x[x<=0.01] = 0.01 # projected 
  }
  return(x)
}

draw_oracles = function(theta) {
  #plot(rowSums(theta$ab), theta$tau, ylim=c(0, 1), xlab = "Predictability", ylab = "Reliability")
  rbPal <- colorRampPalette(c('red','green'))
  colors <- rbPal(10)[as.numeric(cut(theta$tau,breaks = 10))]
  plot(theta$ab[,2], theta$ab[,1], xlab = "beta", ylab = "alpha", pch = 20, asp = 1., col = colors);
  abline(a=0, b=theta$gamma / (1-theta$gamma), col = "red", lwd = 2)  
}

draw_top = function(thetas, gammas = seq(0.3, 0.48, 0.02), top = 16) {
  plot(gammas, sapply(thetas, function (x) x$tau[1]), type = 'o', ylim = c(0,1), xlab = "gamma", ylab = "tau");
  for (i in 2:top) {
    lines(gammas, sapply(thetas, function (x) x$tau[i]), type = 'o', ylim = c(0,1));
  }
  taus=sapply(thetas, function(x,y) x$tau[y], y=1:length(thetas[[1]]$tau));
  return(taus);
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
    w = theta$w[U];
    ab = list(alpha = theta$ab[U,1], 
              beta =  theta$ab[U, 2]);
    tau = theta$tau[U];
    ab_data = alpha_and_beta(I, ab, tau * w);
    Psi[U,1]=Psi[U,1] + digamma(ab_data$alpha);
    Psi[U,2]=Psi[U,2] + digamma(ab_data$beta);
    Psi[U,3]=Psi[U,3] + digamma(ab_data$alpha + ab_data$beta);
    r_data = ratio(I, ab$alpha, ab$beta, theta$gamma, w);
    Tau[U] = Tau[U] + r_data * tau / (r_data * tau + 1 - tau);
    gamma$omega = gamma$omega + sum(I %*% ((1-tau) * w));
    gamma$psi   = gamma$psi   + sum(sum(w) - I %*% w - w);
    Delta[U] = Delta[U] + 1;
  }

  theta$ab = gamma_solve(cbind((Psi[,1]-Psi[,3])/Delta, (Psi[,2]-Psi[,3])/Delta), theta$ab, delta = 1/Delta);
  theta$tau = (Tau + tau0) / (Delta+1);
  #theta$gamma = min(max_gamma, (gamma$omega + gamma0) / (gamma$psi + 1));
  draw_oracles(theta);
  return(theta)
}

# glba_proj = function(estimated, hypergraph, theta, ab=cbind(3,1)) {
#   n=length(hypergraph$I);
#   Psi = as.vector(matrix(0, 1, 3));
#   Delta = 0;
#   for (i in 1:n) {
#     U = hypergraph$inv_oracles[hypergraph$U[[i]]];
#     I = matrix(estimated$I[[i]], nrow=1, ncol=length(U));
#     w = theta$w[U];
#     ab_local = list(alpha = ab[1], beta =  ab[2]);
#     tau = theta$tau[U];
#     ab_data = alpha_and_beta(I, ab_local, tau * w);
#     #print(ab_data)
#     Psi[1] = Psi[1] + digamma(ab_data$alpha);
#     Psi[2] = Psi[2] + digamma(ab_data$beta);
#     Psi[3] = Psi[3] + digamma(ab_data$alpha + ab_data$beta);
#     Delta = Delta + 1;
#     ab = gamma_solve(cbind((Psi[1]-Psi[3])/Delta, (Psi[2]-Psi[3])/Delta), ab, delta = 1/Delta);
#   }
#   return(ab);
# }

glba = function(hypergraph, weight = 1., tau = 0.5, iter = 100, gamma = 0.3, theta0 = NULL) {
  if (!is.null(theta0)) {
    theta = theta0
  } else {
    m = length(hypergraph$oracles);
    theta={};
    theta$ab = cbind(matrix(1., m, 1), matrix(1., m, 1)) / (2*tau0);
    theta$tau = tau * as.vector(matrix(1, m, 1));
    theta$gamma = gamma;
    theta$w = weight * as.vector(matrix(1, m, 1));    
  }
  
  for (i in 1:iter) {
    theta = varEM_update(theta, hypergraph);
    if (i > 50) { # empirical Bayes
      #tau0 <<- median(theta$tau);
      #rate0 <<- 1/median(theta$ab);
    }
  }
  return(theta);
}


glba_curve = function(hypergraph, gammas = seq(0.3, 0.48, 0.02)) {
  thetas = list(length(gammas));
  
  cat(sprintf("gamma: %.2f ... ", gammas[1]))
  thetas[[1]] = glba(hypergraph, gamma = gammas[1]);
  cat(sprintf("[done]\n"))
  for (i in 2:length(gammas)) {
    thetas[[i]] = thetas[[i-1]];
    thetas[[i]]$gamma = gammas[i];
    cat(sprintf("gamma: %.2f ... ", gammas[i]))
    thetas[[i]] = glba(hypergraph, theta0 = thetas[[i]]);
    cat(sprintf("[done]\n"))
  }
  taus=draw_top(thetas, gammas);
  return(list("thetas" = thetas, "taus" = taus));
}
