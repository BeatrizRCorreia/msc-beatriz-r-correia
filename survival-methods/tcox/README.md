# **TCox** using [TCox](https://github.com/sysbiomed/TCox)

The *TCox* model was also used, the same seed as the one used for the previous values was set. The original idea that is introduced in this thesis for this model was to insert a weights penalty factor that would encourage the selection of genes with similar correlation patterns between tumour tissue data from two distinct cancers.

In [Twiner](https://github.com/sysbiomed/twiner), BRCA and PRAD data from both normal tissue and tumour tissue were already analysed using this penalization. However, this was enforced in the context of a classification model using sparse logistic regression, never in the survival context.

As previously addressed at the beginning of the chapter, breast and prostate cancer are two of the most prevalent hormone-dependent cancers. The recognition of common biomarkers to these two cancers has the potential to boost targeted treatments to hormone-dependent disease conditions, such as bone relapse, by using the common gene signatures as targets.

The code developed starts by only keeping in the BRCA data and the PRAD data the genes that are present in both datasets. This leaves both datasets with 18767 features.

Then, the correlation matrix for each of the two datasets is computed and the dissimilarity measure of each gene between the two correlations matrices is calculated from the angle of the corresponding vectors.

The values are normalized by the maximum value so that they become between zero and one and the final distribution of the weights for every gene can be seen on the histogram below.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/histogram_weights.png" alt="histogram_weights" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/histogram_weights.png" alt="histogram_weights" width="700"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/histogram_weights.png" alt="histogram_weights" width="650"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/histogram_weights.png" alt="histogram_weights" width="600"/>
</p>

The information about each iteration execution saves the same eight fields already described but instead of the "measure" with the cross-validation measure (for the initial models), and in substitution of the "min.degree" for the *glmSparseNet* models, this one saves the description of the weights being used and is named **weights.used**. The weights were applied directly, thus w is saved in that field (if the inverse of the weights was being tested, 1/w could be saved in this field).

The next step consists on the fitting of the model. There was the option to fit it to the BRCA data or the PRAD data and the two alternatives were explored. Tables 4.19, 4.20 and 4.21 are referent to the analysis fitted to the breast dataset. Tables 4.22, 4.23 and 4.24 are referent to the prostate dataset.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/table1_BRCA.png" alt="table1_BRCA" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/table2_BRCA.png" alt="table2_BRCA" width="750"/>
</p>

One of the aspects that immediately stands out from Table 4.19 is the number of iterations needed to perform and consequent running times that are considerably bigger for PRAD when compared to BRCA executions. The high degree of censorship in prostate data (98%, *versus* 86% in breast data) is the reason for these differing statistics and a common problem when dealing with survival data.

There are seventeen genes that were selected at least half of the runs and at all sparsity degrees (for every alpha): *ACAP2*, *APOOL*, *ARID1B*, *COMTD1*, *DCTPP1*, *FAM98B*, *JAK1*, *MIER1*, *PCGF5*, *PCYT1A*, *PPFIA3*, *PSENEN*, *PSME2*, *RPL14*, *SECISBP2L*, *VHL* and *ZFC3H1*.

From Table 4.20 with the averages of the evaluation metrics for the 25 successful iterations, the best average p-values for the test set were obtained for alpha 0.2, with an average p-value of 0.571677 (highlighted).

None of the best iterations' p-values for the test set presented in Table 4.21 were significant, but the lowest p-value obtained was 0.075468 alpha 0.1, which is close to being less or equal to 0.05 (highlighted). The survival curve for this iteration is in the figure below Table 4.21.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/table3_BRCA.png" alt="table3_BRCA" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/tcox_brca_w_bestmodel.png" alt="Best TCox model BRCA" width="500"/>
</p>

Regarding the results from applying the penalization weights to the PRAD data solely, the average p-values for the test set, in Table 4.23 were smaller than the ones obtained for the BRCA data, these p-values were closer to being significant and with the best average p-value in the test set being 0.142810 for alpha 0.1 (highlighted).

The demanding execution times forced an adjustment to the number of successful runs that were required for each alpha. The justification for this modification can be recognized from Table 4.22 where from alpha 0.6 to alpha 0.5 the iterations set were dropped from 25 to 10 and the program still took 1.25 more hours to complete.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/table1_PRAD.png" alt="table1_PRAD" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/table2_PRAD.png" alt="table2_PRAD" width="750"/>
</p>

Another remark is that there are four genes that were selected at least half of the runs for every alpha value: *C9orf80*, *QTRTD1*, *REXO4* and *SURF6*.

The best model attained had a p-value in the test set of 0.138457 and is highlighted in Table 4.24. It was for alpha 0.1 and the corresponding survival curve for this iteration is in the figure below the table mentioned.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/tcox/tcox_prad_w_bestmodel.png" alt="Best TCox model PRAD" width="500"/>
</p>
