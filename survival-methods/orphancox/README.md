# **OrphanCox** using [glmsparsenet](https://www.bioconductor.org/packages/release/bioc/html/glmSparseNet.html)

The *glmSparseNet* package in R used one more time, in this occasion to fit OrphanCox models, with the same seed.

This model receives as input the same parameters as HubCox: the **network** option used was again "correlation", the **cutoff** value used was 0.001 and the **min-degree** used was 0.6. All these parameters are again the ones used in the [TCox paper](https://www.mdpi.com/2227-9059/8/11/488) with code available [here](https://github.com/sysbiomed/TCox) to fit OrphanCox models on colorectal cancer survival data from TCGA. The min-degree value 0.8 (used for HubCox) was also experimented with but was not resulting in any successful iterations.

The results for 25 successful iterations for each alpha value, are expressed in the two tables below.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/orphancox/table1.png" alt="table1" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/orphancox/table2.png" alt="table2" width="750"/>
</p>

In the first, Table 4.16 it is possible to see that in general, more runs had to be executed for this model than for the HubCox. This model required on average for all the alphas, 141.64 runs to obtain the 25 considered successful. The first model, with Lasso-Cox, Ridge-Cox and ElasticNet-Cox executed on average 98.91 runs and HubCox only 38.82.

There are three genes: *FIBCD1*, *PGK1* and *WNT3A* that were selected at least half of the runs for every alpha value.

The iterations for alpha equal to 0.1 have an execution time of 20.92 hours which is inconsistent with the other execution times, this results from the network computation time (18 hours) that is also included in that number because this was the first value to be executed for this model.

The average p-values obtained in the test set in Table 4.17 are very close to significant for alphas 0.6, 0.8, 0.9 and 1 with values around 0.06 (highlighted).

Table 4.18 with the best iterations for each alpha value shows significant p-values in the test set for alphas 0.3, 0.6, 0.7, 0.8, 0.9 and 1. Nonetheless, the best iteration has a p-value for the test set of 0.017882 and was collected using alpha 0.6 (highlighted). The corresponding Kaplan-Meier curves of this iteration that show the separation in low and high risk groups is in the figure below.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/orphancox/table3.png" alt="table3" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/orphancox/orphancox_bestmodel.png" alt="Best OrphanCox model" width="500"/>
</p>