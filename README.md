### Guide to scripts and data for "Sequence characteristics distinguish transcribed enhancers from promoters and predict their breadth of activity"
#### Laura Colbran
#### November 30, 2018

### bin/
avg_curves.py\
:   averages ROC and PR curves from many classifiers, used when we split larger sets up


enh-prom_analyses.ipynb\
   contains R code for relative ROC calculation (Fig. 3), TF motif analyses (Fig. 4), PCA, kmer weights (Fig 1)


### data/
all_fantom_enhancers.bed\
   Broad enhancers = all with #tiss >45\
   Context Specific = random subset of those with #tiss = 1\
   regions were set to 600bp     


all_fantom_prom.bed\
   Broad Promoters = random subset of those with #tiss >372\
   Context-Specific = all with #tiss <9\
   regions were set to 600bp


roadmap_enhancers_600bp.bed\
   filtered, set to 600bp
 

roadmap_promoters_600bp.bed\
   filtered, set to 600bp


tf_motif_specificity.csv\
   FANTOM TSPS scores, IDs

### classifiers/
output and scripts from all SVM classifiers

fantom_enhVSprom/\
   direct classifiers between enhancers and promoters (Fig. 1)


broadVSspecific/\
   classifiers between broad and specific regions (Fig. 2)


cgi_analyses/\
   stratified by CGI status (Fig. 3)


roadmap_enhVSprom/\
   direct classifiers between enhancers and promoters (Fig. 5)
