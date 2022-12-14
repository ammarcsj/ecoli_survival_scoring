# Ecoli Survival Scoring

This repository contains the code as well as the necessary data to re-create the Bioinformatics analyses in the manuscript  <em>"Analysis of proteome adaptation reveals a key role of the bacterial envelope in starvation survival."</em> by [Schink et al.](https://www.biorxiv.org/content/10.1101/2022.05.18.492425v1.abstract) To this end, multiple proteomics dataset are analysed in order to determine proteins that show overall correlation with starvation survival.



<figure><img src="data/figures/method_box.jpg" width="70%"><figcaption><em>Summary of the steps performed in scoring the proteins. See manuscript for a detailed description.</em></figcaption></figure>

### Repository contents

* data folder: The datasets used in the analyses are stored in the data folder. Additional README files are provided in the respective folders. The general structure of the data folder comprises: 
    * The search engine results as well as differential analysis results for the three proteomics datasets by [Houser et al.](https://pubmed.ncbi.nlm.nih.gov/26275208/ ) (PRIDE: PXD002140), [Hui et al.](https://pubmed.ncbi.nlm.nih.gov/25678603/) (PRIDE: PXD001467) and [Schmidt et al.](https://pubmed.ncbi.nlm.nih.gov/26641532/) (PRIDE: PXD000498). 
    * Results of the Enrichment Analyses for several types of GO processes
    * Mapping files detailing gene name mappings and mass abundance estimates
    * Results tables of the analyses performed after running the code

* code folder: Contains Jupyter Notebooks (and some helper classes in Python files) that create the results tables of the study.

### Inspecting the data
If you only want to inspect the data, click on the green "Code" button on the upper right and select "Download Zip". After you have downloaded the zip file, you can unpack it and navigate through the data with your standard file browser.

### Re-creating the analyses
If you additionally want to run the code, install Python and the necessary packages as detailed below and then follow the instuctions given in the Jupyter Notebooks. The Notebooks are enumerated and you can follow the order indicated there.

### Code requirements
The code has been tested using Python 3.8 using the additional python packages numpy and pandas as well as jupyter. We recommend installing [the Conda package management system for Python](https://www.anaconda.com/products/distribution) and then creating an environment. This can be done by typing the following commands in the command prompt (after installation of Conda):

```bash
conda create -n ecoli_scoring python=3.8
conda activate ecoli_scoring
pip install pandas
pip install numpy
pip install jupyter
```

After this, the Jupyter notebooks can be executed. In case you are new to Jupyter Notebooks, see for example [this tutorial](https://realpython.com/jupyter-notebook-introduction/).

