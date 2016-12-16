image_conf=c();
high=c();
low=c();
for (pid in unique(emodata$pid)){
uid=which(sdata[pid,]!=0);
scores=sdata[pid,uid];
tmp=rowMeans(result$taus)[hypergraph$inv_oracles[uid]]; 
conf=1-prod(1-rowMeans(result$taus)[hypergraph$inv_oracles[uid]]);
if (sum(tmp*scores)/sum(tmp) < 2 && conf > 0.95) {
  #print(pid);
  low=c(low, pid);
#  print(c(sum(tmp*scores)/sum(tmp)-1, mean(scores)-1))
}
if (sum(tmp*scores)/sum(tmp) > 8 && conf > 0.95) {
  #print(pid);
  high=c(high, pid);
  #  print(c(sum(tmp*scores)/sum(tmp)-1, mean(scores)-1))
}
image_conf=c(image_conf, conf);
}


for (pid in high) {
  for (uid in which(sdata[pid,]!=0)){
    if (sdata[pid,uid] < 5) {
      print(uid);
    }
  }
}

for (pid in low) {
  for (uid in which(sdata[pid,]!=0)){
    if (sdata[pid,uid] > 5) {
      print(uid);
    }
  }
}