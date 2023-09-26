# Dose optimization of an adjuvanted peptide-based personalized neoantigen melanoma vaccine (2023)

**Authors**: Authors: Wencel Valega-Mackenzie, Marisabel Rodriguez Messan, Osman N. Yogurtcu, Ujwani Nukala, Zuben E. Sauna, and Hong Yang

## Overview

This repository provides the code to reproduce the numerical results presented in our paper [citation here]. It includes a Jupyter Notebook named [CanVaxDOpt](CanVaxDOpt.ipynb), which guides users through the process of simulating optimal vaccine doses, immune responses, and tumor responses for each patient.

### Key Features:

- Numerical solution of the optimization problem using the Forward Backward Sweep method (FBSM).
- Selection of the optimal peptide dose with the most clinical benefit than any other tested dose.
- Reproduction of figures from our paper: Optimal vaccine doses (Fig 5), Activated T cells (Fig 6), Tumor diameter (Fig 7).
- Comparison of clinical benefits using the predicted optimal vaccine dose against clinical trial doses or other tested doses (Tables 3 and 4).
- Tables presenting optimal vaccine dose concentrations of peptides and adjuvants per vaccination day (Supplemental Information).

### Usage:

To reproduce the results, execute each cell in the Jupyter Notebook from top to bottom. Note that running cells to find optimal solutions may take approximately 1-2 hours per patient. Other cells should run quickly.

### Dependencies:

Ensure that you have all the necessary Python libraries installed as described on the first cell.

### Excel File:

The Excel file, [Patients_HLA_BA](Patients_HLA_BA.xlsx), contains binding affinities of MHC-I/II to specific alleles associated with the immunogenic peptides of each patient. These affinities were used to compute the dissociation constant, $K_{D,j}$. The relevant information from this spreadsheet has been integrated into the Jupyter Notebook.

### Citation:

If you use this code in your research, please cite our paper [citation here].
