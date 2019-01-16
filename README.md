### Guide to scripts and data for "Sequence characteristics distinguish transcribed enhancers from promoters and predict their breadth of activity"
#### Laura Colbran
#### 01/15/2018

### bin/
avg_curves.py\
&nbsp;&nbsp;&nbsp;&nbsp;averages ROC and PR curves from many classifiers, used when we split larger sets up


enh-prom_analyses.ipynb\
&nbsp;&nbsp;&nbsp;&nbsp;contains R code for relative ROC calculation (Fig. 3), TF motif analyses (Fig. 4), PCA, kmer weights (Fig 1)


### data/
all_fantom_enhancers.bed\
&nbsp;&nbsp;&nbsp;&nbsp;Broad enhancers = all with #tiss >45\
&nbsp;&nbsp;&nbsp;&nbsp;Context Specific = random subset of those with #tiss = 1\
&nbsp;&nbsp;&nbsp;&nbsp;regions were set to 600bp     


all_fantom_prom.bed\
&nbsp;&nbsp;&nbsp;&nbsp;Broad Promoters = random subset of those with mean_act >372\
&nbsp;&nbsp;&nbsp;&nbsp;Context-Specific = all with mean_act <9\
&nbsp;&nbsp;&nbsp;&nbsp;regions were set to 600bp


roadmap_enhancers_600bp.bed\
&nbsp;&nbsp;&nbsp;&nbsp;filtered, set to 600bp


prom_enh_rel_ROC.txt\
&nbsp;&nbsp;&nbsp;&nbsp;values for Fig. 3 relative ROCs
 

roadmap_promoters_600bp.bed\
&nbsp;&nbsp;&nbsp;&nbsp;filtered, set to 600bp


tf_motif_specificity.csv\
&nbsp;&nbsp;&nbsp;&nbsp;FANTOM TSPS scores, IDs

### classifiers/
output and scripts from all SVM classifiers\
N.B. classifier script requires Python 2.7.8 and Shogun Machine Learning Toolbox v4.0.0

fantom_enhVSprom/\
&nbsp;&nbsp;&nbsp;&nbsp;direct classifiers between enhancers and promoters (Fig. 1)

fantom_enhVsprom_cgiMatched/\
&nbsp;&nbsp;&nbsp;&nbsp;direct classifiers between enhancers and promoters, stratified by CGI overlap

broadVSspecific/\
&nbsp;&nbsp;&nbsp;&nbsp;classifiers between broad and specific regions (Fig. 2)


cgi_analyses/\
&nbsp;&nbsp;&nbsp;&nbsp;stratified by CGI status (Fig. 3)


roadmap_enhVSprom/\
&nbsp;&nbsp;&nbsp;&nbsp;direct classifiers between enhancers and promoters (Fig. 5)

### enhVsprom_tf_matching/
tomtom output for top 6-mers in direct classifiers between enhancers and promoters (Fig 3B)

### motif_sim/
tomtom output for top 6-mer in other enhancer and promoter classifiers\
&nbsp;&nbsp;&nbsp;&nbsp;hocomoco/ (Fig 4)\
&nbsp;&nbsp;&nbsp;&nbsp;jaspar/ (Figs S11 & S12)\

### tf_counts/
overall broad and narrow tf counts in regions
