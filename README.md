# IPMiner
Predicting ncRNA-protein interaction using high-level features <br>

we proposed a computational method IPMiner to predict ncRNA-protein interactions from sequences, which made use of deep learning and further improve its performance using stacked 
ensembling. It automatically extracted high-level features from conjoint triad features of protein and RNA sequence using stacked autoencoder, then the high-level features are fed into random forest 
to predict ncRNA-protein 
interaction. Finally stacked ensembling is used to integrate different predictors to further improve prediction performance.

Dependency: <br>
python 2.7 <br>
deep learning lib keras: https://github.com/fchollet/keras/  (version Keras-0.1.2 and the backend theano v0.9)<br>
machine learning lib scikit-learn v0.17: https://github.com/scikit-learn/scikit-learn <br>

Usage:
python IPMiner.py -datatype=RPI488  <br>
where RPI488 is lncRNA-protein interaction dataset, and IPMiner will do 5-fold cross-validation for it. you can also choose other datasets, such as RPI1807, RPI369,, RPI13254 and NPInter. <br>

python IPMiner.py -r=RNA_fasta_file -p=protein_fasta_file   <br>
it will predict pairwise interaction score for RNAs and protiens in input fasta file.

<b> Reference </b> <br>
Xiaoyong Pan, Yong-Xian Fan, Junchi Yan and Hong-Bin Shen. <a href=https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2931-8> IPMiner: hidden ncRNA-protein interaction sequential pattern mining with stacked autoencoder for accurate computational prediction</a>. BMC Genomics. 2016, 17:582 DOI: 10.1186/s12864-016-2931-8.

Contact: xypan172436atgmail.com
