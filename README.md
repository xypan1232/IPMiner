# IPMiner
Predicting ncRNA-protein interaction using high-level features <br>

we proposed a computational method IPMiner to predict ncRNA-protein interactions from sequences, which made use of deep learning and further improve its performance using stacked 
ensembling. It automatically extracted high-level features from conjoint triad features of protein and RNA sequence using stacked autoencoder, then the high-level features are fed into random forest 
to predict ncRNA-protein 
interaction. Finally stacked ensembling is used to integrate different predictors to further improve prediction performance.

Dependency: <br>
deep learning lib keras: https://github.com/fchollet/keras/ <br>
machine learning lib scikit-learn: https://github.com/scikit-learn/scikit-learn <br>

Usage:
python IPMiner.py -datatype=RPI488  <br>
where RPI488 is lncRNA-protein interaction dataset, and IPMiner will do 5-fold cross-validation for it. you can also choose other datasets, such as RPI1807, RPI369, RPI2241, RPI13254 and NPInter. 
