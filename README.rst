
README
======

|

This folder contains the code used in `The evolution of hematopoietic cells under cancer therapy <LINK TO THE PAPER>`_.
If you use this code in a publication, please cite:

.. admonition:: Citation
   :class: note

   Oriol Pich, Albert Cortes-Bullich, Ferran Muiños, Marta Pratcorona, Abel Gonzalez-Perez,
   Nuria Lopez-Bigas, **The evolution of hematopoietic cells under cancer therapy**, 
   Nature Communications, 2021, `doi <https://doi.org/10.1038/s41467-021-24858-3>`_


Part of the data needed to reproduce the analyses is under restricted access. These datasets and their accession are described in specific sections of the methods.
The metastatic tumor cohort data (DR-024 version 2) is available from the `Hartwig Medical Foundation <(https://www.hartwigmedicalfoundation.nl/en>`_ for academic research upon request. Without the Hartwig Metastatic tumor cohort the results won't be reproducible.

The AML data from were obtained from:

- `Genomics of Acute Myeloid Leukemia <https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000159.v11.p5>`_ 
- `WTSI (EGAD00001005028) <https://ega-archive.org/datasets/EGAD00001005028/>`_
- `In-house cohort (EGAS00001005234)  <https://ega-archive.org/datasets/EGAS00001005234/>`_

|


Running this software
*********************

These analysis have been perform using software in Python, R, Julia, and GNU bash.

We have created a set of `Jupyter notebooks <http://jupyter.org/>`_
that you can run if you are interested in re-running partially or
totally our analysis.
In each notebook you will find further details for running them.

|

Requirements
************

To be able to run those notebooks you need to have the following
software installed:

```
conda create -n tAML python pandas numpy matplotlib seaborn tqdm scikit-learn holoviews bgreference ipykernel statsmodels
```

If there are other specific software requirements, they will be specified within each section. 

|


How to reproduce the analyses
*****************************

We divided the analyses in different figures. Please look at the different jupyter notebooks which will point
to the commands needed for reproducing the analyses. The supplementary figures are also included.

- `Figure 1 <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/evolution_hemato_therapy/raw/master/src/figures/Figure1.ipynb>`_

- `Figure 2 <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/evolution_hemato_therapy/raw/master/src/figures/Figure2.ipynb>`_

- `Figure 3 <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/evolution_hemato_therapy/raw/master/src/figures/Figure3.ipynb>`_

- `Figure 4 <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/evolution_hemato_therapy/raw/master/src/figures/Figure4.ipynb>`_

- `Supplementary Note <http://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/evolution_hemato_therapy/raw/master/src/figures/supplementary_note.ipynb>`_

Most of the code will try to find files specificed in the config file. Please change the path to your files accordingly. 

The mutations from our variant calling can be found here:

- `Mutations <https://bitbucket.org/bbglab/evolution_hemato_therapy/raw/master/data/samples/AML.all.muts.gz>`_

The code has been written by `Oriol Pich <https://github.com/oriolpich>`_ , except Figures 4e,d and most of the Logistic Regression analysis which was produced by `Ferran Muiños <https://github.com/koszulordie>`_. 
The code has been revised and tested by `Miguel Grau <https://github.com/migrau>`_ and `David Martinez <https://github.com/dmartmillan>`_.

 
