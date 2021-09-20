> # MSc code

## Survival analysis of transcriptomic high-dimensional oncological data for the identification of cancer biomarkers
___
> ### Repository content

The materials presented in this repository are organized in three folders:

* ``data`` folder with the two datasets that were used:
    * ``BRCA_primary_solid_tumor`` folder with the breast data files:
        * ``data2_01_BRCA_RNASeq2GeneNorm-20160128_transposed.csv``
        * ``data2_01_BRCA_RNASeq2GeneNorm-20160128.csv``
        * ``data2_colData.csv``
        * ``data2_META_0.csv``
        * ``data2_sampleMap.csv``
    * ``PRAD_primary_solid_tumor`` folder with the prostate data files:
        * ``data2_01_PRAD_RNASeq2GeneNorm-20160128_transposed.csv``
        * ``data2_01_PRAD_RNASeq2GeneNorm-20160128.csv``
        * ``data2_colData.csv``
        * ``data2_META_0.csv``
        * ``data2_sampleMap.csv``

* ``survival-methods`` folder with one folder for each of the Survival Analysis Statistical Methods used:
    * ``elasticnet-cox`` folder with the following content:
        * ``elasticnetcox.R``
        * ``results`` folder
        * ``README.md``
        * some README.md figures
    * ``hubcox`` folder with the following content:
        * ``hubcox.R``
        * ``results`` folder
        * ``README.md``
        * some README.md figures
    * ``orphancox folder`` with the following content:
        * ``orphancox.R``
        * ``results`` folder
        * ``README.md``
        * some README.md figures
    * ``rsf`` folder with the following content:
        * ``rsf.py``
        * ``results`` folder
        * ``README.md``
        * some README.md figures
    * ``tcox`` folder with the following content:
        * ``tcox.R``
        * ``results`` folder
        * ``README.md``
        * some README.md figures

* ``other-functions`` folder with the following files:
    * ``auxiliary_functions.R``
    * ``chat_tool_figures.R``
    * ``extract-dataset.R``
    * ``venndiagram.R``
___
> ### Model comparison in terms of biomarkers selected

In Table 4.31 can be seen the alpha values that allowed to obtain between 44 and 50 genes extracted at least half of the runs. These genes, also presented in that table were combined in the Venn diagram below.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/all_models_table.png" alt="all_models_table" width="650"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/venn_diagram.png" alt="venn_diagram" width="400"/>
</p>

The biomarkers obtained from the TCox fitted to the BRCA data and fitted to the PRAD data have no gene in common, for this reason and for simplification the genes from the two sets are combined in TCox area from the diagram (total of 94 genes).

There are **eight genes** that stand out because they were selected by **four different models**. Four of them were selected in common by the ElasticNet-Cox, HubCox, OrphanCox and RSF models. These genes are *FIBCD1*, *XRCC4*, *HSPA8* and *WNT3A*. The other four were selected in common by the ElasticNet-Cox, HubCox, OrphanCox and TCox (more specifically they result from the fitting to the BRCA data), and these biomarkers are *JAK1*, *PCYT1A*, *RPL3* and *PPFIA3*.

There are also **twenty four genes** that were selected by **three distinct models**. Twenty of them were selected in common by the ElasticNet-Cox, HubCox and OrphanCox models: *LSG1*, *TCP1*, *LOC220729*, *IMP5*, *SERPINA3*, *GLUL*, *IRF2*, *LOC100128977*, *AKR1E2*, *IYD*, *GCET2*, *PCMT1*, *DHX16*, *ABCG4*, *TPT1*, *MED17*, *SHCBP1*, *PELO*, *PGK1* and *IL18*. Two of them extracted from the three models: ElasticNet-Cox, HubCox and the TCox: *GPR172A* and *SPINT1*. Another two extracted from these three: ElasticNet-Cox, OrphanCox and TCox: *DCTPP1* and *C11orf20*.

Another important fact is that there is an overlap of thirteen genes between the ElasticNet-Cox and OrphanCox and, a smaller overlap of four genes between ElasticNet-Cox and HubCox, but there are no genes that were picked from both HubCox and OrphanCox. This makes sense since these network-based models promote opposing structures (hubs *versus* orphans), already addressed.

The [CHAT tool](https://chat.lionproject.net) can be used to verify if the genes selected were already described in the literature in the context of Cancer Hallmarks. The set of genes retrieved from each model (in Table 4.31) were looked up in this tool. In the six figures below are the number of hits found of a hallmark for each of the genes for ElasticNet-Cox, HubCox, OrphanCox, TCox fitted in the BRCA data, TCox fitted in the PRAD data and RSF models.

Apparently, the TCox model has more biomarkers that had no confirmation of hallmark hits than any of the other models (in total 42 of the 94 genes had no matches). Nevertheless, genes that do not appear to have hallmarks matches can only be an indication that those have not been deeply studied yet.

Counts of hallmarks association with the 50 genes selected by the **ElasticNet-Cox** model (the darker the cell, the higher number of hits were found of a hallmark for a particular gene in the CHAT database):

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/hallmarks_elasticnetcox.png" alt="hallmarks_elasticnetcox" width="700"/>
</p>

Counts of hallmarks association with the 44 genes selected by the **HubCox** model (the darker the cell, the higher number of hits were found of a hallmark for a particular gene in the CHAT database):

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/hallmarks_hubcox.png" alt="hallmarks_hubcox" width="700"/>
</p>

Counts of hallmarks association with the 45 genes selected by the **OrphanCox** model (the darker the cell, the higher number of hits were found of a hallmark for a particular gene in the CHAT database):

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/hallmarks_orphancox.png" alt="hallmarks_orphancox" width="700"/>
</p>

Counts of hallmarks association with the 47 genes selected by the **TCox** model fitted in the **BRCA** data (the darker the cell, the higher number of hits were found of a hallmark for a particular gene in the CHAT database):

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/hallmarks_tcox_brca.png" alt="hallmarks_tcox_brca" width="700"/>
</p>

Counts of hallmarks association with the 47 genes selected by the **TCox** model fitted in the **PRAD** data (the darker the cell, the higher number of hits were found of a hallmark for a particular gene in the CHAT database):

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/hallmarks_tcox_prad.png" alt="hallmarks_tcox_prad" width="700"/>
</p>

Counts of hallmarks association with the 50 genes selected by the **Random Survival Forest** model (the darker the cell, the higher number of hits were found of a hallmark for a particular gene in the CHAT database):

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/hallmarks_rsf.png" alt="hallmarks_rsf" width="700"/>
</p>

___
> ### How to run the code in this repository

**Linux environment instructions to run the files**

To execute the **ElasticNet-Cox model**:

1. Open a terminal.
2. Get to the folder "elasticnet-cox" where the file "elasticnetcox.R" is.
3. Run the command ``Rscript elasticnetcox.R | tee ~/Documents/msc-beatriz-r-correia/survival-methods/elasticnet-cox/results/dev_alpha0.1_test.txt``. This will execute the "elasticnetcox.R" file and both print the results from the execution in the terminal and to a file named "dev_alpha0.1_test.txt" (this file name is an example and can be adjusted in the command).

To execute the **HubCox model**:

1. Open a terminal.
2. Get to the folder "hubcox" where the file "hubcox.R" is.
3. Run the command ``Rscript hubcox.R | tee ~/Documents/msc-beatriz-r-correia/survival-methods/hubcox/results/alpha0.1_test.txt``. This will execute the "hubcox.R" file and both print the results from the execution in the terminal and to a file named "alpha0.1_test.txt" (this file name is an example and can be adjusted in the command).

To execute the **OrphanCox model**:

1. Open a terminal.
2. Get to the folder "orphancox" where the file "orphancox.R" is.
3. Run the command ``Rscript orphancox.R | tee ~/Documents/msc-beatriz-r-correia/survival-methods/orphancox/results/alpha0.1_test.txt``. This will execute the "orphancox.R" file and both print the results from the execution in the terminal and to a file named "alpha0.1_test.txt" (this file name is an example and can be adjusted in the command).

To execute the **Random Survival Forest model**:

1. Open a terminal.
2. Get to the folder "rsf" where the file "rsf_final.py" is.
3. Run the command ``python3 rsf_final.py``.

To execute the **TCox model**:

1. Open a terminal.
2. Get to the folder "tcox" where the file "tcox.R" is.
3. Run the command ``Rscript tcox.R | tee ~/Documents/msc-beatriz-r-correia/survival-methods/tcox/results/BRCA_alpha0.1_w_test.txt``. This will execute the "tcox.R" file and both print the results from the execution in the terminal and to a file named "BRCA_alpha0.1_w_test.txt" (this file name is an example and can be adjusted in the command).
