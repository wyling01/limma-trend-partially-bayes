# limma-trend-partially-bayes
R-scripts needed to reproduce the Limma Trend Partially Bayes paper.

The repository is structured as follows:

- **function**: Folder that contains functions required for reproducing the paper. 
- **real_datasets**: Folder that contains the four real data examples in the paper.
    - **rnaseq1_melanoma**: Contain the data, codes, and generated figures used for Section 7.1.
    - **rnaseq2_malaria**: Contain the data, codes, and generated figures used for Section 7.2. The data we used is from [Tonkin-Hill et al. (2018)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2004328), and the preprocessing part of `rnaseq2_malari.R` is adapted from the paper's [repository](https://github.com/gtonkinhill/falciparum_transcriptome_manuscript/blob/master/all_gene_analysis/TextS1_all_gene_analysis.Rmd). 
    - **chipseq**: Contain the data, codes, and generated figures used for Section 7.3. The data and preprocessing steps used in `chipseq.R` can be found in [MAnorm2 [Tu et al. (2021)]](https://github.com/tushiqi/MAnorm2/blob/master/utility/code-MAnorm2Paper/limmaTrend.r).
    - **protemics**: Contain the data, codes, and generated figures used for Section 7.4. The data and preprocessing steps used in can be found in [DEqMS [Zhu et al. (2020)]](https://bioconductor.statistik.tu-dortmund.de/packages/3.18/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html). 
- **sec1_plot**: Folder for introductory trend plots in Section 1.1.
- **simulation**: Folder for simulation studies in Section 6. In `simulation.R` we assess the theoretical results with our proposed methods, and in `sim_comparison.R`, we further compare with MAnorm2 and MAP methods. The generated simulation results are stored in folder `data`.