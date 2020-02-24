# bp205b-project

Assigns cell barcode to corresponding multiseq barcode
Applies Monocle Trajectory Analysis to cell barcoded expression data

## Jupyter notebooks

* data-analysis.ipynb - contains the functions necessary for mapping, as well as run_all function (takes in read1 and read2 gz files, cell barcode whitelist, multiseq barcode whitelist). Generates mapping csv and plots
* specific files - has files filled in, generates specific csv and plot when run
  * data-analysis-aa.ipynb
  * data-analysis-mk.ipynb
  * data-analysis-gbc.ipynb

## R Script

* TJ_Analysis_Clean.R - contains functions required for Monocle workflow. Applies clustering, trajectory analysis, DEx Gene branch analysis.

## Subfolders

* dataset folders - read1 and read2 gz files, whitelisted barcodes
  * aa_dataset
  * mk_dataset
  * gbc_dataset
* results - generated plots and mapping csvs

