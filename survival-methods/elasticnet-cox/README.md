# **ElasticNet-Cox** using [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)

The *glmnet* package in R was used to fit Lasso, Ridge and Elastic-net regularization for the Cox model. This package penalizes the negative logarithm of the partial likelihood with an elastic-net penalty.

For these results the seed "1997" was set. The program developed allows to set the number of iterations to perform, the cross-validation measure and the alpha value to be used.

The lambda value controls how much penalty to apply to the regularization. The **cross-validation measure** helps to choose the lambda value that gives the best model from the 100 random values that are fitted by default. There are two options for this measure in *glmnet* with regularized Cox regression: **deviance** and **C**. If the first is selected, the lambda value chosen is the one that gives the model with the lowest partial-likelihood. If **C** is picked, the lambda value chosen is the one that delivers the model with the highest C-index. This measure, already explained previously, only considers comparable pairs (one pair is not comparable if one of the observations is right censored at a time before the other observation's death time).

The **alpha** is the value between zero and one that is used to weight the contributions of the $$l_{1}$$-norm and with one minus alpha the $$l_{2}$$-norm.

After the first trials, it was possible to verify that some iterations do not return any genes. In order to guarantee that executing the program will return some selected features, a **number of iterations to perform** can be set. 25 was the number of iterations chosen for the **deviance** measure runs and 50 was the number of iterations set for the **C** measure, as the executions with the last were always successful (returning selected covariates) and faster.

Each iteration's information is saved in a *list* with all the runs executed. The information of the iteration also consists in a *list* with eight fields: "measure", to save the cross-validation measure, "alpha", for the alpha value, "nr.selected.genes" with the number of genes selected, "list.of.genes" a *list* with the genes that were selected, "pvalue.train" with the p-value in the train set, "pvalue.test" for the p-value in the test set, "cindex.train" and "cindex.test" with the C-indexes for the train and test sets.

The p-values calculated result from the log-rank test via the Kaplan-Meier estimator. To achieve this, the patients were divided into two groups by the median of the fitted relative risk. The two Kaplan-Meier curves generated from an iteration (train set and test set) are added to a pdf that in the end of the program execution has the plots from every iteration.

The *list* that has the set of *lists* with the iterations executed allows to posteriorly calculate averages and standard deviations of the runs executed. The standard deviations are in parenthesis after the values in the tables that are presented in the following pages. The program also returns the iteration with the lowest p-value train and the iteration with the lowest p-value test.

Another functionality is the computation of the genes that were selected in at least half of the runs and the exact number of times they were, in fact, picked. This calculation uses the *hash* R library that allows to work with *dictionaries*. The structure saves the gene name in the key and the number of iterations where it was selected in the value part.

In Tables 4.7 and 4.8 can be seen the seen the results from 25 runs for each alpha values ranging from 0 to 1 in 0.1 intervals, using partial-likelihood as cross-validation measure for the lambda tuning.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/elasticnet-cox/table1_dev.png" alt="table1_dev" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/elasticnet-cox/table2_dev.png" alt="table2_dev" width="750"/>
</p>

As expected, and is displayed on 4.7, as the alpha value increases, the number of features selected dramatically reduces. The alpha equal to zero corresponds to the Ridge-Cox model, at this value, there was no variable selection, the genes selected are genes in the preprocessed data that was used (19045). As the alpha penalty increases towards the alpha equal to one, that corresponds to the Lasso-Cox model, the sparsity of the model increases as some of the coefficients for the covariates are forced to zero.

When genes are correlated, Lasso will choose between them based on its performance in the particular data sample being used, meaning that with different train and test sets a different feature from the set of correlated features could have been picked. With this in mind, it is relevant to note that the eleven executions, selected, in common, at least half of the runs the *PGK1* gene.

The evaluation of the models obtained was done using the average of the p-values and C-indexes in both train and test sets and can be seen on Table 4.8. Concerning the average p-values, in the test set, the models that are closer to being considered statistically significant (p-value less than or equal to 0.05) had p-values around 0.06. These were obtained with alphas between 0.7 and 1 (highlighted).

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/elasticnet-cox/table3_dev.png" alt="table3_dev" width="750"/>
</p>

In the last table for this measure, Table 4.9, is presented information about the best iterations in the 25 runs performed. The best iteration is the one that yields the smallest p-value in the test set, in this case, the lowest p-value obtained was 0.057965 for alpha 0.9. The four genes that were selected were *FIBCD1*, *LOC100128977*, *PGK1* and *PPFIA3*. The survival curve for this model can be seen below.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/elasticnet-cox/encox_dev_bestmodel.png" alt="Best ElasticNet-Cox model with dev" width="500"/>
</p>

All these runs and executions were important to understand the behaviour of the models and inspect the number of genes that could be selected. Besides the executions presented, more executions were made in order to extract a similar number of genes from the different models. These are omitted for simplification.

In the end, the objective set was to extract between 44 and 50 genes because these numbers were the closest numbers that was possible to collect of **genes selected at least half of the iterations** across all models experimented. The rationale behind this is based on a consensus approach, in which the systematic selection of a gene may indeed indicate biological and clinical relevance.

The same procedure was executed using the C-index as cross-validation measure for the lambda tuning. The results from 50 runs for eleven different alphas can be seen on Tables 4.10 and 4.11.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/elasticnet-cox/table1_c.png" alt="table1_c" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/elasticnet-cox/table2_c.png" alt="table2_c" width="750"/>
</p>

In Table 4.10 one of the aspects that stand out is the fact that all runs were successfully returning a set of genes. The number of genes extracted bigger: the mean of the average number of genes selected using the partial likelihood was 1737.41 for all the alpha values and 2224.88 was the mean of the average number of genes selected using the C-index for all the alpha values. The execution times were also smaller: the mean of the execution times for all alpha values using the partial likelihood was 7.89 hours and for all alpha values using the C-index the mean was 2.92 hours.

Some of the iterations using this measure, resulted in C-indexes of one for the train sets, however this value is not returned by the C-index function being used. Therefore, in Table 4.11 the averages of the C-indexes for the train set (the only set where some iterations yielded NA) are only computed for the iterations that returned a value different from NA. The corresponding standard deviations are also affected.

Concerning the p-values in the test set, the alpha values from 0.1 to 0.7 resulted in significant (less than or equal to 0.05) p-value averages (highlighted).

Table 4.12 presents the best iterations, again with the lowest p-values in the test partitions for this measure. As explained, the NA values represent C-indexes of one. The best iteration of them all was obtained for alpha equal to 0.1, with a p-value in the test set of 0.000207 (highlighted). This model's survival curves for low and high risk groups can be seen below.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/elasticnet-cox/table3_c.png" alt="table3_c" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/elasticnet-cox/encox_cindex_bestmodel.png" alt="Best ElasticNet-Cox model with c" width="500"/>
</p>

These executions, using the C-index measure, extracted a minimum of 250 genes selected at least half of the runs and, therefore, although the results were really good (because of the many significant p-values obtained) since there was the aim to extract between 44 and 50 genes (a smaller set), this model was not taken further to the comparison with other models as a proof-of-concept choice. Nonetheless, the list of the 250 genes can be checked in the /results folder.
