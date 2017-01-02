# GLBA
This repository includes the R code for Gated Latent Beta Allocation (GLBA) and some other codes that implement baseline methods, which appeared in the following paper:

Probabilistic Multigraph Modeling for Improving the Quality of Crowdsourced Affective Data
<i>Jianbo Ye, Jia Li, Michelle Newman, Reginald B. Adams, Jr., James Z. Wang</i>, 
IEEE Transactions on Affective Computing (To appear)

## Dataset

This repository compiles a (ID,Score) dataset `filter_four.csv` for demonstration purpose, which is a derivative of dataset EmoSet described in X. Lu's dissertation work. The use of this demo set is restricted only to explore the functionality of code.  

X. Lu, “Visual characteristics for computational prediction of aesthetics and evoked emotions,” Ph.D. dissertation, The Pennsylvania State Uni- versity, 2015, chapter 5. [Online available](https://etda.libraries.psu.edu/catalog/28857)

## How to run

Load model
```
> GLBA
```

Load data (default: valence, `the_metric = 1`)
```
> preprocssing
```
where `hypergraph` stores the userid map `hypergraph$oracles` and inverse map `hypergraph$inv_oracles`. 

Train GLBA to obtain user reliability
```
> result=glba_curve(hypergraph)
```
where the user reliability can be calculated by `rowMeans(result$taus)`. 


Remember valence is the best dimension to obtain reliability score. Later, if you want to load arousal data, you don't need to retrain. Change `the_metric = 2` in `preprocssing.r` to load arousal data.


Obtain image confidences and scores
```
> experiment2
```


