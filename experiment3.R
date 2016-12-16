
par(mfrow=c(1,2))
#subject centric overheads

taus = rowMeans(result$taus);

x = seq(0, 1, 0.01)
y = c()
for (t in x) {
  s = hypergraph$oracles[taus <= t]
  y=c(y,nnzero(sdata[,s]))
}

plot(x, y, 'l', xlab = 'Threshold of Subject Reliability', ylab = 'Overheads', lwd = 3)
#legend(0, 40000, c("Subject Centric Overheads"), lwd = 3, box.lwd = 0)


#image centric overheads

x = seq(0.4, 1, 0.01)
y = c()
for (t in x) {
  s = hypergraph$I[which(image_conf <=t)]
  y = c(y, Reduce("+", lapply(s, nrow)))
}

plot(x, y, 'l', xlab = 'Threshold of Image Confidence', ylab = 'Overheads', lwd = 3)
#legend(0.4, 40000, c("Image Centric Overheads"), lwd = 3, box.lwd = 0)