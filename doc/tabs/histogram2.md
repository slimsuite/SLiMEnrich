The histogram shows the expected distribution of predicted DMIs from these randomised PPI data, and marks the observed number of [predicted DMIs](https://github.com/slimsuite/SLiMEnrich/wiki/Analysis-and-Outputs/#predicted) in the real data. The following values are calculated and displayed:

* **P-value:** this is an empirical p-value based on the proportion of randomised PPI datasets that equal or exceed the number of DMI observed in the real data.
* **Enrichment:** this is the ratio of observed non-redundant DMI to the mean of the random non-redundant DMI. 
* **FDR:** this is the estimated proportion of observed DMI that are false positives, based on the mean random DMI count, excluding any random datasets exceeding the observed predicted DMI count.

Clicking **SETTINGS** will open options for customising the histogram. The bin size (default 1) and x-axis extension can also be set for the commandline version (`-b` and `-x`, respectively). Note that the histogram x-axis will not truncate before the maximum real or random DMI count and any **Extend X-axis End** setting below this number will be ignored. The histogram can be saved as a PNG by clicking **DOWNLOAD**.

---
