# **Random Survival Forest** using [PySurvival](https://square.github.io/pysurvival/models/random_survival_forest.html)

Lastly another non-parametric method was explored. The Random Survival Forest model was tested using Python and a package named *PySurvival*. The seed "5" was set.

After executing the same preprocessing steps, the idea here was to analyse the effects of the tuning of three different parameters in order to obtain the best model for the breast data.

The two evaluation metrics chosen to scrutinize the performance of this model were the C-index and the IBS. The C-index will allow to compare this model with the previous models and the IBS (that compares the predicted survival rate with the actual status of patients) was picked because it is another statistic usually reported in survival analysis and that had not been experimented for the former models.

The best models are the ones that obtain the highest C-indexes in the test sets and that obtain an IBS score in the test set below 0.25 as close to zero as possible.

The only parameter that was fixed for all the experiments was **max features**, this establishes the number of drawn candidate features (genes) to consider when looking for the best split and was set to the square root of the number of genes, the default value in several software packages. Nonetheless, as a starting point, the first experiment also used a value of 953 (5\%) for the max tree depth and 10 for the min node size.

In the first experiment, the aspect analysed was the impact of the **number of trees** in the forest. The idea was to compare different numbers of trees according to the C-index obtained in the test set to find the best value for this parameter. The limit number of trees that was experimented with was 1600 since beyond this value the program was killed multiple times due to extensive memory usage.

The results from two iterations for each number of trees can be seen in Table 4.25. The model that obtained the highest C-index in the test set (highlighted) built 1600 trees. This number makes sense since more trees allow to build more tree arrangements for the genes and, it is still a number very lower than the number of genes (19045). The number of trees was fixed as 1600 for the next tuning experiments.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/rsf/table1.png" alt="table1" width="750"/>
</p>

Each one of the executions for a number of trees produces a table with the VIMP for every gene. This information from the table was converted into a dictionary with the name of the gene in the key part and the VIMP as the value part. Then, the dictionaries for each execution were combined by making an average of the VIMPs of every iteration for a specific number of trees. Finally, another average was performed, this time taking into account the averages for each number of trees. The averages were considered the most appropriate form of combining the VIMPs in order to establish more stable values.

The results can be seen in Table 4.25: there were 10481 genes with positive VIMPs, 8522 genes with negative VIMPs and 42 genes with VIMPs equal to zero that can result from not being used for any tree (never seen) or the average of the averages of VIMPs yielding zero (extremely unlikely).

The second aspect considered was the impact of the **max tree depth**. For this argument the idea was to examine the effects of the forest having trees with a maximum depth of 5% of the dataset (around 953 genes) to 40% of the dataset (7618 genes). The reason for this was based on the fact that the dataset has 19045 genes and a value smaller that 5% seemed to not include enough information and beyond 40% was not possible to execute due to the high degree of memory usage. Again two iterations for each value were carried out and the results are in Table 4.26. The value of 7618 allowed to achieve the best C-index in the test set (highlighted) and was fixed for the further experiments.

From the VIMPs and, again using the same procedure described for the number of trees, it was interesting to notice that there were still 29 genes with VIMP equal to zero, that were most likely never positioned in the trees.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/rsf/table2.png" alt="table2" width="750"/>
</p>

The last parameter inspected was the effect of **min node size** values. This argument takes the minimum number of samples needed to be at a leaf node. Here, in Table 4.27, values between 2 and 20 were tested but it started to be questioned if the number of iterations being performed (two) could not be enough to choose a parameter as "best". For this reason the results were repeated, this time performing ten iteration using min node size 2, 10 and 18. The results are in Table 4.28 where the optimal value was found to be 10, as it originated the highest average C-index (highlighted).

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/rsf/table3.png" alt="table3" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/rsf/table4.png" alt="table4" width="750"/>
</p>

The last experiment was to perform thirty iterations of the RSF model using the parameters tuned in the previous attempts (considered best). The results of this execution are in Table 4.29. The average C-index was not as high as the ones that was possible to produce in the past executions. In the figure below table 4.29 are depicted the train set and test set C-indexes from all the iterations executed here. The values oscilate a lot for the test set.

Another aspect to note from Table 4.29 is that nine genes still have VIMP equal to zero and were probably never used in the trees.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/rsf/table5.png" alt="table5" width="750"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/rsf/RSF_bestmodel_30iterations.png" alt="RSF_bestmodel_30iterations" width="750"/>
</p>

Using the dictionary with the VIMPs resulting from the averaging of the thirty iterations of the best model (procedure already explained), the figure below was created. By looking at the image, there was the idea of retrieving the genes above a VIMP threshold, these are considered to be the most important genes. After testing with different thresholds in order to retrieve a number of genes similar to the numbers obtained for the other models tested (between 44 and 50 genes), the threshold of 1.8 was chosen.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/survival-methods/rsf/RSF_vimps_bestmodel_30iterations.png" alt="RSF_vimps_bestmodel_30iterations" width="750"/>
</p>

All the genes that have VIMPs above 1.8 (50 genes) were gathered: *CEACAM5*, *STXBP5*, *ENC1*, *DNAJC22*, *GLT25D2*, *CEL*, *SCG5*, *HPDL*, *CALML5*, *XRCC4*, *FIBCD1*, *WNT3A*, *PHB*, *HSPA8*, *DUS1L*, *RBBP8*, *CCDC54*, *C6orf141*, *GPR37L1*, *ZNF654*, *DAB2*, *CEACAM1*, *MLLT6*, *ADAMTSL1*, *ADAMTS7*, *TNFSF4*, *LOC728323*, *DNAJC14*, *NACC1*, *ABHD10*, *GADL1*, *NLE1*, *MYO16*, *OPLAH*, *ZNF707*, *SPRY4*, *MGC21881*, *CA11*, *MTMR7*, *NXN*, *MRPL38*, *BCLAF1*, *ZCCHC9*, *CHERP*, *PLXNB2*, *PUF60*, *TSGA10*, *LHX5*, *BAZ2A* and *C21orf57*.
