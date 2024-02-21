# Scyliorhinus_canicula_popgen

This repo contains code used to calculate differentiation between Atlantic and Mediterranean populations of small-spotted catshark (*Scyliorhinus canicula*).

The analysis steps (and where the code for each is found) are as follows:
  1. Download data: **code/main_analysis.sh**
  2. Per-sample mapping an variant calling: **code/main_analysis.sh**
  3. Joint genotyping and variant filtering: **code/main_analysis.sh**
  4. Population genetic statistic calculation: **code/main_analysis.sh**
  5. Population structure analysis: **code/main_analysis.sh** and **code/pop_structure_plots.R**
  6. Genome statistic calculation: **code/main_analysis.sh**
  7. Outlier analysis: **code/pcadapt.R** and **code/main_analysis.sh**
  8. Genome-wide plotting: **code/manhattan_plots.R** and **code/manhattan2.R**. *n.b. manhattan2.R is a modified version of the manhattan function from Stephen Turner's* [qqman](https://github.com/stephenturner/qqman) *package*

Additional metadata files needed for the analysis are in **data/**
