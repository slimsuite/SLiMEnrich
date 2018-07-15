#### SLiMEnrich Histogram

This tab is the main SLiMEnrich output. Clicking on this tab will trigger the [PPI randomisation](https://github.com/slimsuite/SLiMEnrich/wiki/Randomisation) and generate the expected distribution of predicted DMIs based on purely chance associations between PPI `mProteins` and `dProteins`. In brief, the input PPI dataset is shuffled without replacement (i.e. keeping the original number of interacting partners per protein) and the number of predicted DMI calculated for the shuffled data. By default, this is done 1000 times, which can be set using the **Number of randomisations** (or `--random=INTEGER` on the commandline). 

