# smFRET_corr_analysis
Python implementation of the multichannel correlation analysis framework described by M. Dhar and M. A. Berg [J. Chem. Phys. 163, 184113 (2025)] (DOI: [10.1063/5.0284658](https://doi.org/10.1063/5.0284658)) for smFRET application

## Original study

Mainak Dhar, Mark A. Berg; Nonparametric analysis of noisy, multivariable time series using high-order correlation functions: Single-molecule FRET as an example. J. Chem. Phys. 14 November 2025; 163 (18): 184113. DOI: [10.1063/5.0284658](https://doi.org/10.1063/5.0284658)

## Input dataset
The input dataset used in this study can be found from Zenodo [10.5281/zenodo.5701309](https://doi.org/10.5281/zenodo.5701309), file 'exp_dataset_Fig2_1ms.zip'.

These data have been provided by Ben Schuler (University of Zurich; [0000-0002-5970-4251](https://orcid.org/0000-0002-5970-4251)) in the frame of an smFRET tools benchmark study (Götz, M., Barth, A., Bohr, S.SR. *et al.* A blind benchmark of analysis tools to infer kinetic rate constants from single-molecule FRET trajectories. Nat Commun 13, 5402 (2022)) [10.1038/s41467-022-33023-3](https://doi.org/10.1038/s41467-022-33023-3).

The original data study by Zosel *et al.* has been published in PNAS: F. Zosel, A. Soranno, K.J. Buholzer, D. Nettels, & B. Schuler,  Depletion interactions modulate the binding between disordered proteins in crowded environments, Proc. Natl. Acad. Sci. U.S.A. 117 (24) 13480-13489, DOI: [10.1073/pnas.1921617117](https://doi.org/10.1073/pnas.1921617117) (2020).

The dataset probes the interactions between the nuclear-coactivator binding domain of CBP/p300 (NCBD) and the intrinsically disordered activation domain of the steroid receptor coactivator 3 (ACTR), measured using confocal single-photon detection; the data have been binned to 1 ms. The dataset comprises 19 molecules with average series length of 120 s.

## Google Colab notebook
The correlation analysis pipeline is available from Google Colab as an interactive Jupyter Notebook: [click here to access](https://colab.research.google.com/github/molcretb/smFRET_corr_analysis/blob/main/pipeline_Dhar_Berg_smFRET_correlation_framework.ipynb)

It contains the full instructions on how to conduct the analysis step by step.

## License
The codes associated with this repository are licensed under the MIT license. The original dataset from the smFRET tools benchmark study is licenced under the CC-BY-4.0 license

## Contact

Bastien Molcrette [0000-0002-5995-5376](https://orcid.org/0000-0002-5995-5376) (Schmid group, Department of Chemistry, University of Basel) bastien.molcrette [@] unibas.ch
