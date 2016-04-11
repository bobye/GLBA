# precision-recall curve


num_of_oracles=size(hypergraph$oracles);
oracle_labels=as.vector(array(0,c(num_of_oracles,1)));
oracle_labels[hypergraph$inv_oracles[labeled_spammers]]=1;

library(ROCR)
theta=glba(hypergraph, gamma = 0.3);
pred=prediction(1-theta$tau, oracle_labels)
perf1=performance(pred, "prec", "rec");


theta=glba(hypergraph, gamma = 0.37);
pred=prediction(1-theta$tau, oracle_labels)
perf2=performance(pred, "prec", "rec");


theta=glba(hypergraph, gamma = 0.44);
pred=prediction(1-theta$tau, oracle_labels)
perf3=performance(pred, "prec", "rec");

pred=prediction(1-rowMeans(result$taus), oracle_labels)
perf=performance(pred, "prec", "rec");

plot(perf, col=1, lwd=2)
par(new=T)
plot(perf1, col=2, lwd=2)
par(new=T)
plot(perf2, col=3, lwd=2)
par(new=T)
plot(perf3, col=4, lwd=2)
legend(.8, .8, c("aver.","γ=0.3", "γ=0.37", "γ=0.44"), col=c(1,2,3,4), lty=1, lwd=2, box.lwd=0)


plot(perf, col=1, lwd=2, ylim=c(0, 1.0), xlim=c(0, 1.0))
par(new=T)
plot(perf_arousal, col=2, lwd=2, ylim=c(0, 1.0), xlim=c(0, 1.0))
par(new=T)
plot(perf_dominance, col=3, lwd=2, ylim=c(0, 1.0), xlim=c(0, 1.0))
par(new=T)
plot(perf_likeness, col=4, lwd=2, ylim=c(0, 1.0), xlim=c(0, 1.0))
legend(.7, .9, c("valence","arousal", "dominance", "likeness"), col=c(1,2,3,4), lty=1, lwd=2, box.lwd=0)
