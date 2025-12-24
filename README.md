# Multichannel Correlation analysis applied to single-molecule FRET

## Overview
Python implementation of the multichannel correlation analysis framework described by M. Dhar and M. A. Berg [J. Chem. Phys. 163, 184113 (2025)] (DOI: [10.1063/5.0284658](https://doi.org/10.1063/5.0284658)) for smFRET application.

The goal of this project is to reproduce the results presented in the previously mentioned Dhar & Berg study (referred as "original study" in the following), within a Python framework. The practical implementation is largely based on the [Supplementary Information](https://doi.org/10.60893/figshare.jcp.30406921) associated with the original study, which provides a detailed section with pseudocodes implementation. Once the results of the original study will have been reproduced, the next step will be to apply the analysis framework to other smFRET datasets and compare the results to the existing litterature.

!!!

*Right now (2025-12-24), the implementation of the framework into Python is not fully complete; I cannot reproduce the figures of the original study. Debugging of the scripts is ongoing...*

!!!

## Original study

Mainak Dhar, Mark A. Berg; Nonparametric analysis of noisy, multivariable time series using high-order correlation functions: Single-molecule FRET as an example. J. Chem. Phys. 14 November 2025; 163 (18): 184113. DOI: [10.1063/5.0284658](https://doi.org/10.1063/5.0284658)

Supplementary information can be found from FigShare (contained guidelines for the pratical implementation of Dhar & Berg framework): [DOI: 10.60893/figshare.jcp.30406921](https://doi.org/10.60893/figshare.jcp.30406921)

## Resources
### Input dataset
The input dataset used in this study can be found from Zenodo [10.5281/zenodo.5701309](https://doi.org/10.5281/zenodo.5701309), file 'exp_dataset_Fig2_1ms.zip'.

These data have been provided by Ben Schuler (University of Zurich; [0000-0002-5970-4251](https://orcid.org/0000-0002-5970-4251)) in the frame of an smFRET tools benchmark study (Götz, M., Barth, A., Bohr, S.SR. *et al.* A blind benchmark of analysis tools to infer kinetic rate constants from single-molecule FRET trajectories. Nat Commun 13, 5402 (2022)) [10.1038/s41467-022-33023-3](https://doi.org/10.1038/s41467-022-33023-3).

The original data study by Zosel *et al.* has been published in PNAS: F. Zosel, A. Soranno, K.J. Buholzer, D. Nettels, & B. Schuler,  Depletion interactions modulate the binding between disordered proteins in crowded environments, Proc. Natl. Acad. Sci. U.S.A. 117 (24) 13480-13489, DOI: [10.1073/pnas.1921617117](https://doi.org/10.1073/pnas.1921617117) (2020).

The dataset probes the interactions between the nuclear-coactivator binding domain of CBP/p300 (NCBD) and the intrinsically disordered activation domain of the steroid receptor coactivator 3 (ACTR), measured using confocal single-photon detection; the data have been binned to 1 ms. The dataset comprises 19 molecules with average series length of 120 s.

### Google Colab notebook
The correlation analysis pipeline is available from Google Colab as an interactive Jupyter Notebook: [click here to access](https://colab.research.google.com/github/molcretb/smFRET_corr_analysis/blob/main/pipeline_Dhar_Berg_smFRET_correlation_framework.ipynb)

It contains the full instructions on how to conduct the analysis step by step.

### Content of the GitHub repository
This repository contains the following resources:
* correlation_module.py: Python functions to calculate the correlations and primed-correlations functions used in the analysis. It corresponds to the section "SII. PSEUDOCODE, A.Correlation Module" of the Supplementary Information associated with the original study.
* inversion_module.py: Python functions used to retrieve the FRET states space from the computed primed-correlation functions. It corresponds to the section "SII. PSEUDOCODE, B.Inversion Module" of the Supplementary Information associated with the original study.
* pipeline_Dhar_Berg_smFRET_correlation_framework.ipynb: Jupyter notebook used to run the analysis pipeline from the extraction of the smFRET traces from the input dataset to the recovering of the states space. This Jupyter Notebook can either be run locally or on Google Colab.
* selection_data_traces.csv: list of the smFRET data files used as input dataset. Each row corresponds to the path of one of the data files.

## License
The codes associated with this repository are licensed under the MIT license. The original dataset from the smFRET tools benchmark study is licenced under the CC-BY-4.0 license

## Contact

Bastien Molcrette [0000-0002-5995-5376](https://orcid.org/0000-0002-5995-5376) (Schmid group, Department of Chemistry, University of Basel) bastien.molcrette [@] unibas.ch
