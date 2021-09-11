# **HubCox** using [glmsparsenet](https://www.bioconductor.org/packages/release/bioc/html/glmSparseNet.html)

The *glmSparseNet* package in R was used to fit HubCox models. This package uses *glmnet*.

The same seed, "1997" was set. The code is available in [hubcox.R](https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/hubcox/hubcox.R) and also allows to define the number of iterations to perform and the alpha value to be used as previously addressed for the last models. The new parameters that can be chosen are the cutoff, the min-degree and the network.

The **network** parameter is responsible for the generation of the network. There are three options: "correlation", "covariance" and "sparsebn" (which uses a sparse bayesian network), but it is also possible to use as network an adjacency matrix or a metric calculated for each node. The "correlation" was the option selected.

The **cutoff** value is a real number that represents a threshold to remove edges from the network, 0.005 was used.

The **min-degree** defines the minimum value that the weights for the features can have, this was set to 0.6.

Both the network, the cutoff and the min-degree were set to those values because these allowed to obtain a set of genes with TCGA data for colorectal cancer survival data ([TCox paper](https://www.mdpi.com/2227-9059/8/11/488) with code available [here](https://github.com/sysbiomed/TCox)). When they were tested for the breast data, they also manage to select genes.

The network calculation step that involves the network parameter and the cutoff value is the one that takes the most time to complete. This step took a total of 18 hours to complete. For this reason these values were fixed and the network was only calculated one time.

Regarding the information saved in the *list* about each run, this program instead of saving the "measure" (as for the previous models), attaches the "min.degree" value that was set. The network or the cutoff were not saved because no other values of these were tried due to the extensive time needed for the computation.

In the two Tables below can be seen the results from requesting 25 successful runs for the same alpha values as for the other models.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/hubcox/table1.png" alt="table1" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/hubcox/table2.png" alt="table2" width="750"/>
</p>

The standard deviations, in parenthesis, in the first table allow to see that the average number of genes selected did not vary much between the 25 iterations of every alpha, except for alpha equal to 0.1, where the standard deviation was 53.56 genes, and also excluding alpha equal to zero where all the variables are selected, the average of the standard deviations for the other alpha values (0.2 to 1) was 1.36 genes.

It is also relevant that there are four genes: *LOC100128977*, *PGK1*, *PPFIA3* and *SHCBP1* that were selected at least half of the runs for every alpha value.

In the second table can be seen that the best models, with the lowest average p-values in the test set were obtained for alpha 0.9 and 1 (highlighted). However, none of the averages has a value considered significant (less or equal to 0.05).

In the last table are the best iterations from the 25 for every alpha value and the best model was obtained for alpha 0.9, with a p-value in the test of 0.062519, very close to being considered significant (highlighted). The three genes selected in this iteration were *LOC100128977*, *PGK1* and *PPFIA3*. The corresponding survival curves for low and high risk groups for this model is presented in below.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/hubcox/table3.png" alt="table3" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/hubcox/hubcox_bestmodel.png" alt="Best HubCox model" width="500"/>
</p>
