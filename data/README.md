# Data used
___
> ### BRCA data

One of the TCGA datasets downloaded is from the breast invasive carcinoma study (BRCA is the study abbreviation).

This study includes 14 experiments, from which the RNASeq2GeneNorm-20160128 experiment was selected. This experiment was chosen because it has gene mRNA abundance obtained using RNA-Seq data. The data is upper quartile normalized in RSEM TPM gene expression values. Then, the primary solid tumour assay (sample type code 01) was chosen. A fraction of the four csv files that were downloaded can be seen in Tables 4.1, 4.2, 4.3 and 4.4. There are 20501 genes.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/data/tables_BRCA_example.png" alt="tables_BRCA_example" width="650"/>
</p>

The data has information about 1093 patients. These patients are mostly females and, the majority 941 (86.1\%) are alive/censored (their vital status is a zero in the dataset), opposing to only 152 (13.9%) death patients (event of interest occurred, therefore their vital status is a one in the dataset). The patients were also diagnosed in very different age groups. The details can be seen in the infographics of Figure 4.1. Finally, these subjects had their pathological diagnosis between 1988 and 2013, most in the 2000s, as can be seen in Figure 4.2.

The BRCA data is the core data that was analysed through different models in this thesis. Breast cancer was chosen because it is the most commonly diagnosed cancer, with an estimated 2.3 million new cases (11.7% of the cancers diagnosed) in 2020, according to the Global Cancer Statistics of 2020. This [article](https://acsjournals.onlinelibrary.wiley.com/doi/epdf/10.3322/caac.21660) is based on the GLOBOCAN estimates for incidence and mortality worldwide for 36 cancers in 185 countries (an online database).

All the models described were fitted with the BRCA dataset.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/data/BRCA_infographics_1.png" alt="BRCA_infographics_1" width="650"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/data/BRCA_infographics_2.png" alt="BRCA_infographics_2" width="650"/>
</p>

___
> ### PRAD data

The other TCGA dataset downloaded is from the prostate adenocarcinoma study (PRAD is the study abbreviation).

This study is composed by 11 experiments, from which the RNASeq2GeneNorm-20160128 was the one selected. This experiment, as previously described for the BRCA data extracted, has gene mRNA abundance obtained from RNA-Seq data and is also upper quartile normalized in RSEM TPM gene expression values. Identically to what was done for the BRCA, the primary solid tumour assay (sample type code 01) was selected. The four files downloaded have the same structure as the ones presented in Tables 4.1, 4.2, 4.3 and 4.4 and therefore are omitted. There are also 20501 genes.

The dataset contains information on 497 patients, all males. Nearly all, specifically 487 (98%) are alive/censored and the remaining 10 (2%) are death. The range of ages at diagnosis present is smaller than in BRCA data. This analysis can be seen in Figure 4.3.

Concerning the pathological diagnosis, these were between 2000 and 2013, most in the later years, as can be seen in the Figure 4.4.

The PRAD data was extracted along with the BRCA data because [breast and prostate cancer are more similar than different](https://pubmed.ncbi.nlm.nih.gov/20147902/). They both are carcinomas that appear in hormonally regulated tissues, even though they arise in organs with very distinct anatomy and physiological function in women and men.

More specifically, in the genetics context, the [mutations of the BRCA2 gene have been linked to both breast and prostate cancer, and it has been proven that mutation carriers with prostate cancer have lower survival rates than patients with prostate cancer who do not carry the BRCA2 mutations](https://pubmed.ncbi.nlm.nih.gov/17565157/).

There are numerous other studies around the proximity between the two cancers, one of them, published in 2019, concluded that [men with a family history of female breast cancer in their first-degree relatives had a higher risk of developing prostate cancer](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-6055-9).

Lastly it is very important to note that breast cancer is the most common invasive cancer in women and prostate cancer is the second most common invasive cancer in men, according to the [Global Cancer Statistics of 2020](https://acsjournals.onlinelibrary.wiley.com/doi/epdf/10.3322/caac.21660).

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/data/PRAD_infographics_1.png" alt="PRAD_infographics_1" width="650"/>
</p>

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/data/PRAD_infographics_2.png" alt="PRAD_infographics_2" width="650"/>
</p>

___
> ### BRCA data preprocessing

There are two variables of interest for survival analysis: the vital status and the survival times of the subjects. The first one can be retrieved from one of the columns in the ColData.csv file as can be seen in Figure 4.1. The failure times for each patient can be calculated using the maximum of four ColData.csv columns: *days_to_death*, *days_to_last_followup*, *Days.to.date.of.Death* and *Days.to.Date.of.Last.Contact*. After doing this calculation, patients with survival times equal to zero or negative were eliminated. This left the dataset with 1080 patients (13 were dropped). In this 1080, there are 928 censored subjects (86%) and 152 death (14%).

Afterwards, the remaining dataset was randomly divided in 80% of patients for the train set and 20% of patients for the test set. The division had in consideration the 86-14, censored-death proportion, as can be seen in Figure 4.5. In this figure, the xtrain and the xtest mentioned refer to the part of the data that has the gene expression of the patients. The ytrain and ytest parts have the information about the survival times and the vital status of the patients. These separations of the data are implementation related.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/data/train_test_splits_BRCA.png" alt="train_test_splits_BRCA" width="650"/>
</p>

Next, the xtrain data is centered (the mean is subtracted) and scaled (divided by the standard deviation). Then, the xtest data is scaled based on the mean and standard deviation from the xtrain data.

After this scaling, some of the xtrain data predictor variables that have near zero variance are dropped (1456 genes) because they are considered to have less predictive power. The same variables are removed from the xtest data.

___
> ### PRAD data preprocessing

The same lineup of preprocessing actions executed for the breast data were followed for the prostate data.

Nonetheless, the name of the four columns in the ColData.csv file used for the failure times calculation, using the maximum were: *patient.days_to_death*, *patient.days_to_last_followup*, *days_to_death* and *days_to_last_followup* for this dataset.

There were no patients with survival times equal to zero or negative, therefore, the data continued with 497 patients. There are 487 censored subjects (98%) and 10 death (2%).

The next step was the random division of the data in 80% of patients for the train set and 20% for the test set. This division had awareness of the 98-2, censored-death proportion and can be seen in Figure 4.6.

The xtrain data, again with the gene expression (as explained for the BRCA), is also centered and scaled, and then used to scale the xtest data (also gene expression).

The next step was the removal from the xtrain data of the features that have near zero variance. Exactly 1482 genes were removed from both xtrain and xtest partitions due to their less predictive power.

<p align="center">
  <img src="https://github.com/BeatrizRCorreia/msc-beatriz-r-correia/blob/main/data/train_test_splits_PRAD.png" alt="train_test_splits_PRAD" width="650"/>
</p>
