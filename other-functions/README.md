# **Other functions**

This folder contains 4 files:

* **extract-dataset.R:** A file with two functions for TCGA data extraction using the [curatedTCGAData](https://bioconductor.org/packages/release/data/experiment/html/curatedTCGAData.html) Bioconductor package:
    * `extract_data()`
    * `extract_primary_solid_tumor_data()`

* **auxiliary_functions.R:** A file with the following functions mainly for model analysis:
    * `get_genes_selected()`
    * `calculate_standard_deviation()`
    * `average_number_of_genes_selected()`
    * `average_pvalue_train()`
    * `average_pvalue_train_significant()`
    * `average_pvalue_train_non_significant()`
    * `average_pvalue_test()`
    * `average_pvalue_test_significant()`
    * `average_pvalue_test_non_significant()`
    * `run_with_lowest_pvalue_train()`
    * `run_with_lowest_pvalue_test()`
    * `c.index.calculator()`
    * `c.index.average.train()`
    * `c.index.average.test()`

* **chat_tool_figures.R:** A file that allows to create a figure with the counts of hallmarks association with the genes selected by the models. The darker the cell the higher the number of hits were found of a hallmark for a particular gene in the [CHAT database](https://chat.lionproject.net). The code presented in this file was adapted from [here](https://github.com/sysbiomed/glmSparseNet/blob/master/R/external_apis.R). The small adjustments are due to the website, at the time of this project, having an expired SSL certificate.

* **venn_diagram.R:** A file with the code used to make a Venn diagram with the genes retrieved from the different models.
