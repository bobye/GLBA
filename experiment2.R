image_conf=c();
for (pid in unique(emodata$pid)){
uid=which(sdata[pid,]!=0);
scores=sdata[pid,uid];
tmp=rowMeans(result$taus)[hypergraph$inv_oracles[uid]]; 

#if ((sum(tmp*scores)/sum(tmp) - mean(scores)) > 1.5) {
#  print(pid);
#  print(c(sum(tmp*scores)/sum(tmp)-1, mean(scores)-1))
#}
image_conf=c(image_conf, 1-prod(1-rowMeans(result$taus)[hypergraph$inv_oracles[uid]]));
}